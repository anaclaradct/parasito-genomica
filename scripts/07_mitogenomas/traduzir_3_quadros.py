#!/usr/bin/env python3
# Traduz uma sequencia de nucleotideos nos 3 quadros de leitura (+1,+2,+3) usando o codigo
# genetico mitocondrial de invertebrados (transl_table=5), sem depender de nenhum alinhamento.
# Conta stop codons internos por quadro: um quadro sem nenhum stop e candidato ao quadro real;
# se NENHUM dos 3 quadros vier limpo, e evidencia forte de pseudogene/NUMT (Crainey et al. 2018).
# Uso: traduzir_3_quadros.py entrada.fasta > saida.tsv
import sys
bases = "TCAG"
codons = [a + b + c for a in bases for b in bases for c in bases]
aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
TAB = dict(zip(codons, aas))
TAB["AGA"] = "S"; TAB["AGG"] = "S"; TAB["ATA"] = "M"; TAB["TGA"] = "W"
def traduzir(seq, frame):
    seq = seq.upper().replace("-", "")[frame:]
    prot = []
    for i in range(0, len(seq) - 2, 3):
        prot.append(TAB.get(seq[i:i+3], "X"))
    return "".join(prot)
def ler_fasta(c):
    d, n = {}, None
    for l in open(c):
        l = l.strip()
        if l.startswith(">"): n = l[1:].split()[0]; d[n] = []
        elif n: d[n].append(l)
    return {k: "".join(v) for k, v in d.items()}
seqs = ler_fasta(sys.argv[1])
print("sequencia\tquadro\tn_stop_interno\tn_aa\tmelhor_quadro\tproteina")
for nome, s in seqs.items():
    resultados = []
    for f in (0, 1, 2):
        p = traduzir(s, f)
        resultados.append((f, p.count("*"), len(p), p))
    melhor = min(resultados, key=lambda r: r[1])
    for f, nstop, naa, p in resultados:
        marca = "SIM" if f == melhor[0] and nstop == melhor[1] else "-"
        print("%s\t%d\t%d\t%d\t%s\t%s" % (nome, f + 1, nstop, naa, marca, p))
