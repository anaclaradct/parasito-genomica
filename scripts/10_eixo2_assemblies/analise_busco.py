#!/usr/bin/env python3
"""Rode em ~/Mansonella_projeto/08_EIXO2:  python3 analise_busco.py
Le busco_compleasm_odb10/*/nematoda_odb10/full_table.tsv + lineages/nematoda_odb10/lengths_cutoff.
Gera resultados_busco/ (tabelas TSV + figura PNG)."""
import glob, os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
BASE = 'busco_compleasm_odb10'; LIN = 'lineages/nematoda_odb10/lengths_cutoff'; OUT = 'resultados_busco'
N50 = dict(Avit_ref=14021571, Dimm_ref=15667105, Bpah_ref=10892846, Bmal_ref=14214749, Mpe_mperst1=45005,
           Mde_mdeux2=97646, Mde_mdeux3=173569, Moz_Brazil1=37748, Moz_Venz1=38906, Mpe_Cam1=285814, Mpe_Cam2=29133)
CTRL = ['Avit_ref', 'Dimm_ref', 'Bpah_ref', 'Bmal_ref']
os.makedirs(OUT, exist_ok=True)
ln = pd.read_csv(LIN, sep='\t', header=None, names=['Gene', 'a', 'b', 'len']).set_index('Gene')['len']
st = {}
for g in N50:
    f = f'{BASE}/{g}/nematoda_odb10/full_table.tsv'
    if os.path.exists(f): st[g] = pd.read_csv(f, sep='\t').drop_duplicates('Gene').set_index('Gene')['Status']
    else: print('FALTA', g, file=sys.stderr)
M = pd.DataFrame(st); M['len'] = ln.reindex(M.index)
M.to_csv(f'{OUT}/genes_status.tsv', sep='\t')
gs = list(st); grp = {g: ('controle' if g in CTRL else 'Mansonella') for g in gs}
rows = []
for g in gs:
    c = M[g].value_counts(); n = len(M)
    rows.append(dict(genoma=g, grupo=grp[g], N50=N50[g], **{k: int(c.get(k, 0)) for k in ['Single', 'Duplicated', 'Fragmented', 'Interspersed', 'Missing']}))
R = pd.DataFrame(rows); R['pct_completo'] = 100 * (R.Single + R.Duplicated) / n
R['pct_frag'] = 100 * R.Fragmented / n; R['pct_ausente'] = 100 * R.Missing / n; R['pct_frag_ausente'] = R.pct_frag + R.pct_ausente
R.sort_values('N50').to_csv(f'{OUT}/resumo_busco.tsv', sep='\t', index=False, float_format='%.2f')
print(R.sort_values('N50')[['genoma', 'N50', 'pct_completo', 'pct_frag', 'pct_ausente']].round(2).to_string(index=False))
rho = pd.Series(np.log10(R.N50)).rank().corr(R.pct_frag_ausente.rank()); print(f'\nSpearman(log N50, %frag+ausente) todos={rho:.2f} (n={len(R)}; exploratorio)')
m = R[R.grupo == 'Mansonella']; print(f'Spearman so Mansonella={m.N50.rank().corr(m.pct_frag_ausente.rank()):.2f} (n={len(m)})')
# fragmentacao por tamanho do gene (5 faixas por quantil)
edges = np.quantile(M['len'].dropna(), np.linspace(0, 1, 6)); M['faixa'] = np.clip(np.searchsorted(edges, M['len'], side='right') - 1, 0, 4)
lab = [f'{int(edges[i])}-{int(edges[i+1])}' for i in range(5)]; T = []
for g in gs:
    for b in range(5):
        s = M.loc[M.faixa == b, g]; T.append(dict(genoma=g, faixa=lab[b], n=len(s), frac_frag=(s == 'Fragmented').mean(), frac_frag_ausente=s.isin(['Fragmented', 'Missing']).mean()))
T = pd.DataFrame(T); T.to_csv(f'{OUT}/fragmentacao_por_tamanho.tsv', sep='\t', index=False, float_format='%.4f')
# classes de genes (ok = Single/Duplicated)
ok = M[gs].isin(['Single', 'Duplicated']); ctrl = [g for g in gs if g in CTRL]; man = [g for g in gs if g not in CTRL]
bad_man = (~ok[man]); k = bad_man.sum(axis=1)
cls = pd.Series('sem problema', index=M.index)
cls[(k > 0) & ok[ctrl].all(axis=1)] = 'problema em alguns Mansonella, ok nos controles (assembly?)'
cls[(k == len(man)) & ok[ctrl].all(axis=1)] = 'problema em TODOS os Mansonella, ok nos controles (divergencia/perda?)'
cls[(k > 0) & ~ok[ctrl].all(axis=1)] = 'problema tambem em controle (geral)'
out = M[['len']].copy(); out['classe'] = cls; out['genomas_com_problema'] = bad_man.apply(lambda r: ','.join(r.index[r]), axis=1)
out[cls != 'sem problema'].sort_values('classe').to_csv(f'{OUT}/classes_genes.tsv', sep='\t'); print('\n', cls.value_counts().to_string())
# figura
C1, C2 = '#1b6b73', '#8a8f98'; fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
OFF = {g: [(5, 6), (5, -12), (-8, 10), (6, 0)][i % 4] for i, g in enumerate(R.sort_values('N50').genoma)}
for g, r in R.set_index('genoma').iterrows():
    a = ax[0]; a.scatter(r.N50, r.pct_frag_ausente, s=55, color=C1 if r.grupo == 'Mansonella' else C2, edgecolor='white', zorder=3)
    a.annotate(g.replace('_ref', ''), (r.N50, r.pct_frag_ausente), xytext=OFF[g], textcoords='offset points', fontsize=8, color='#333')
ax[0].set_xscale('log'); ax[0].set_xlabel('N50 do assembly (pb, escala log)'); ax[0].set_ylabel('% dos 3131 genes BUSCO fragmentados ou ausentes'); ax[0].set_title('A  Contiguidade e completude', loc='left', fontsize=11)
x = np.arange(5)
for g in man: s = T[T.genoma == g]; ax[1].plot(x, 100 * s.frac_frag_ausente, color=C1, alpha=.35, lw=1.2)
ax[1].plot(x, 100 * T[T.genoma.isin(man)].groupby('faixa', sort=False).frac_frag_ausente.mean().reindex(lab).values, color=C1, lw=2.6, label='Mansonella (media; finas = cada assembly)')
ax[1].plot(x, 100 * T[T.genoma.isin(ctrl)].groupby('faixa', sort=False).frac_frag_ausente.mean().reindex(lab).values, color=C2, lw=2.6, label='Controles (media)')
for nome, L in [('NMT', 472), ('gon-2', 1991)]:  # linhas verticais dentro do grafico
    b = int(np.clip(np.searchsorted(edges, L, side='right') - 1, 0, 4)); ax[1].axvline(b + (0.12 if nome == 'gon-2' else -0.12), color='#111', lw=1, ls=(0, (3, 3)), zorder=1)
    ax[1].text(b + (0.12 if nome == 'gon-2' else -0.12), 0.55, f' {nome} ({L} aa)' if nome == 'NMT' else f'{nome} ({L} aa) ', transform=ax[1].get_xaxis_transform(), rotation=90, va='bottom', ha='right' if nome == 'gon-2' else 'left', fontsize=8)
ax[1].set_xticks(x); ax[1].set_xticklabels(lab, fontsize=8); ax[1].set_xlabel('Tamanho de referência do gene BUSCO (aa)'); ax[1].set_ylabel('% fragmentados ou ausentes'); ax[1].legend(frameon=False, fontsize=8); ax[1].set_title('B  Fragmentação por tamanho de gene', loc='left', fontsize=11)
for a in ax: a.spines[['top', 'right']].set_visible(False); a.grid(axis='y', color='#e6e6e6', lw=.6); a.set_axisbelow(True)
plt.tight_layout(); plt.savefig(f'{OUT}/figura_busco.png', dpi=200); print(f'\nOK: {OUT}/ (resumo_busco.tsv, fragmentacao_por_tamanho.tsv, classes_genes.tsv, genes_status.tsv, figura_busco.png)')
