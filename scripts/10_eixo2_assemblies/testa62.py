#!/usr/bin/env python3
"""Rode em ~/Mansonella_projeto/08_EIXO2:  python3 testa62.py
hmmsearch dos BUSCO ausentes em todos os Mansonella (62) e de 62 genes-controle (Single nos 11 genomas) nos proteomas anotados."""
import subprocess, tempfile, pandas as pd
L = 'lineages/nematoda_odb10'
cut = {l.split()[0]: float(l.split()[-1]) for l in open(f'{L}/scores_cutoff') if l.strip()}
c = pd.read_csv('resultados_busco/classes_genes.tsv', sep='\t', index_col=0)
M = pd.read_csv('resultados_busco/genes_status.tsv', sep='\t', index_col=0)
alvo = list(c.index[c.classe.str.contains('TODOS')])
ctl = list(M.drop(columns='len').eq('Single').all(axis=1).pipe(lambda s: s[s].index).to_series().sample(len(alvo), random_state=1))
P = {'Bmal(controle)': 'Bmal_prot/ncbi_dataset/data/GCF_000002995.4/protein.faa'}
for n in ['Moz_Brazil1', 'Moz_Venz1', 'Mpe_Cam1', 'Mpe_Cam2']: P[n] = f'proteomas_sinha/{n}.proteins.fasta'
tmp = tempfile.mkdtemp()
def rodar(ids):
    rows = []
    for g in ids:
        r = {'gene': g, 'corte': cut.get(g)}
        for n, f in P.items():
            t = f'{tmp}/{g}_{n}.tbl'
            subprocess.run(['hmmsearch', '--noali', '--cpu', '2', '-o', '/dev/null', '--tblout', t, f'{L}/hmms/{g}.hmm', f], check=True)
            sc = [float(x.split()[5]) for x in open(t) if not x.startswith('#')]
            r[n] = max(sc) if sc else 0.0
        rows.append(r)
    return pd.DataFrame(rows)
A = rodar(alvo); B = rodar(ctl)
A.assign(conjunto='ausentes_em_todos_Mansonella').to_csv('resultados_busco/teste_62_proteomas.tsv', sep='\t', index=False, float_format='%.1f')
B.assign(conjunto='controle_Single_em_todos').to_csv('resultados_busco/teste_controle_proteomas.tsv', sep='\t', index=False, float_format='%.1f')
print(f'Genes com hit >= corte BUSCO no proteoma (n={len(A)} por conjunto)\n{"proteoma":16s} {"62 ausentes":>12s} {"62 controle":>12s}')
for n in P: print(f'{n:16s} {(A[n] >= A.corte).sum():12d} {(B[n] >= B.corte).sum():12d}')
