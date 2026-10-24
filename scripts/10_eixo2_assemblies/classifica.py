#!/usr/bin/env python3
"""Uso: python3 classifica.py refs.faa genoma1.fna genoma2.fna ...
Le <genoma>.gff (miniprot) na pasta atual. Saida: classificacao.tsv"""
import sys, re, os
from collections import defaultdict
def lens_fa(p):
    d, k = {}, None
    for l in open(p):
        if l[0] == '>': k = l[1:].split()[0]; d[k] = 0
        else: d[k] += len(l.strip())
    return d
refs = lens_fa(sys.argv[1]); END = 5000  # distancia (pb) p/ considerar "na borda do contig"
nome = {'P46548':'NMT','Q93971':'gon-2','G5ECM4':'gtl-1','Q20766':'gtl-2'}
def curto(r):
    for k, v in nome.items():
        if k in r: return v
    return r
out = open('classificacao.tsv', 'w')
out.write('\t'.join('genoma gene n_aln cobertura_uniao melhor_cob identidade contigs dist_borda_pb frameshift_stop estado'.split()) + '\n')
for fna in sys.argv[2:]:
    g = os.path.basename(fna).rsplit('.', 1)[0]; cl = lens_fa(fna); al = []
    if not os.path.exists(g + '.gff'): continue
    for l in open(g + '.gff'):
        f = l.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'mRNA': continue
        a = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
        t = a['Target'].split(); qs, qe = int(t[1]), int(t[2])
        al.append(dict(ref=t[0], c=f[0], s=int(f[3]), e=int(f[4]), sc=float(f[5]), qs=qs, qe=qe,
                       idn=float(a['Identity']), fs=('Frameshift' in a) or ('StopCodon' in a), rank=int(a['Rank'])))
    # resolve parálogos: cada locus pertence ao alinhamento de maior score
    al.sort(key=lambda x: -x['sc']); keep = []
    for x in al:
        if any(k['c'] == x['c'] and min(k['e'], x['e']) - max(k['s'], x['s']) > 0.5 * (x['e'] - x['s']) for k in keep): continue
        keep.append(x)
    for r in refs:
        L = refs[r]; h = sorted([x for x in keep if x['ref'] == r], key=lambda x: -x['sc']); sel = []; cv = set()
        for x in h:  # so fragmentos colineares: residuos novos >=50 e sobreposicao <=30%
            iv = set(range(x['qs'], x['qe'] + 1))
            if not sel or (len(iv - cv) >= 50 and len(iv & cv) <= 0.3 * len(iv)): sel.append(x); cv |= iv
        h = sel
        if not h: out.write(f'{g}\t{curto(r)}\t0\t0\t0\t-\t-\t-\t-\tausente (sem hit)\n'); continue
        cov = set()
        for x in h: cov.update(range(x['qs'], x['qe'] + 1))
        u = len(cov) / L; best = max((x['qe'] - x['qs'] + 1) / L for x in h)
        dist = min(min(x['s'] - 1, cl[x['c']] - x['e']) for x in h)
        fs = any(x['fs'] for x in h); ctg = sorted({x['c'] for x in h}); idn = sum(x['idn'] for x in h) / len(h)
        if best >= 0.9: est = 'integro' + (' (frameshift/stop: checar)' if fs else '')
        elif len(ctg) > 1 and u >= 0.8: est = f'fragmentado entre {len(ctg)} contigs (checar colinearidade)'
        elif len(h) > 1 and u >= 0.8: est = 'multiplos alinhamentos no mesmo contig (checar)'
        elif dist < END: est = 'parcial na borda do contig (provavel truncamento do assembly)'
        else: est = 'parcial interno (divergencia/modelo ruim: checar)'
        out.write(f'{g}\t{curto(r)}\t{len(h)}\t{u:.2f}\t{best:.2f}\t{idn:.2f}\t{",".join(ctg)}\t{dist}\t{fs}\t{est}\n')
out.close(); print(open('classificacao.tsv').read())
