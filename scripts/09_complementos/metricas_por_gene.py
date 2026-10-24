#!/usr/bin/env python3
"""
metricas_por_gene.py - variabilidade por gene e concordancia de cada arvore de gene
com a arvore concatenada (mitogenoma).

(1) Por alinhamento (um FASTA alinhado por gene): comprimento, sitios constantes,
    variaveis, parcimonia-informativos (>=2 estados, cada um em >=2 sequencias) e
    proporcoes. Gaps e ambiguos (N, R, Y..., X) contam como dado ausente.
(2) Opcional, se houver arvores por gene: RF de cada gene contra a concatenada
    (bruto e colapsado nos nos de baixo suporte) e uma matriz splits x genes
    (P = presente e suportado, p = presente e fraco, - = ausente, . = nao aplicavel, ? = sem rotulo).

Arvores por gene podem vir de:
  --arvores-dir DIR [--sufixo-arvore .treefile]   (um arquivo por gene: <gene><sufixo>)
  --arvores-multi ARQ --nomes-genes g1,g2,...      (uma arvore por linha, na MESMA ordem
                                                    do arquivo de particoes do IQ-TREE -S)
IMPORTANTE: confira a ordem dos genes contra o arquivo de particoes antes de interpretar.

Exemplo:
  python3 metricas_por_gene.py --alinhamentos-dir 04_alinhamentos_genes/ \
     --ref matriz_codificante.treefile --arvores-multi loci.treefile \
     --nomes-genes @genes_ordem.txt --excluir "Brugia|NC_004298" --saida resultados_por_gene
"""
import argparse
import glob
import os
import sys
import re
from filo_utils import (ler_fasta, ler_newick, ler_arvore_arquivo, renomear_folhas, folhas,
                        podar_arvore, obter_splits, filtrar_suportados, rf, restringir_split,
                        passa, fabrica_normalizador, carregar_mapa, ACC_RE)

VALIDOS_AA = set("ACDEFGHIKLMNPQRSTVWY")
VALIDOS_NT = set("ACGT")


def detectar_tipo(seqs):
    tot = sum(len(s) for s in seqs)
    nt = sum(sum(1 for c in s if c in "ACGTUNRYKMSWBDHV-?") for s in seqs)
    return "nt" if tot and nt / tot >= 0.9 else "aa"


def metricas_alinhamento(regs):
    seqs = [s for _, s in regs]
    L = {len(s) for s in seqs}
    if len(L) != 1:
        raise ValueError("sequencias com comprimentos diferentes (nao alinhado?)")
    tipo = detectar_tipo(seqs)
    validos = VALIDOS_NT if tipo == "nt" else VALIDOS_AA
    if tipo == "nt":
        seqs = [s.replace("U", "T") for s in seqs]
    cte = var = pi = com_dados = 0
    for col in zip(*seqs):
        cont = {}
        for c in col:
            if c in validos:
                cont[c] = cont.get(c, 0) + 1
        if not cont:
            continue
        com_dados += 1
        if len(cont) == 1:
            cte += 1
        else:
            var += 1
            if sum(1 for v in cont.values() if v >= 2) >= 2:
                pi += 1
    n_cel = len(seqs) * len(seqs[0])
    n_val = sum(1 for s in seqs for c in s if c in validos)
    return {"tipo": tipo, "n_seqs": len(seqs), "comprimento": len(seqs[0]),
            "sitios_com_dados": com_dados, "constantes": cte, "variaveis": var,
            "pars_informativos": pi,
            "prop_variaveis": var / com_dados if com_dados else float("nan"),
            "prop_pi": pi / com_dados if com_dados else float("nan"),
            "pct_gap_ou_ambiguo": 100.0 * (1 - n_val / n_cel) if n_cel else float("nan")}


def lista_arg(v):
    if v.startswith("@"):
        with open(v[1:]) as fh:
            return [l.strip() for l in fh if l.strip()]
    return [x.strip() for x in v.split(",") if x.strip()]


def preparar_arvore(raiz, norm, excluir):
    renomear_folhas(raiz, norm)
    if excluir:
        rx = re.compile(excluir)
        raiz = podar_arvore(raiz, lambda n: bool(rx.search(n)))
    return raiz


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alinhamentos-dir", required=True)
    ap.add_argument("--padroes", default="*.fasta,*.fa,*.fas,*.aln")
    ap.add_argument("--ref", help="arvore concatenada (Newick, com SH-aLRT/UFBoot se houver)")
    ap.add_argument("--arvores-dir")
    ap.add_argument("--sufixo-arvore", default=".treefile")
    ap.add_argument("--arvores-multi")
    ap.add_argument("--nomes-genes")
    ap.add_argument("--sh-min", type=float, default=80.0)
    ap.add_argument("--ufb-min", type=float, default=95.0)
    ap.add_argument("--excluir", help="regex de taxons a ignorar em alinhamentos E arvores (ex.: outgroup)")
    ap.add_argument("--min-ramo", type=float, default=1e-5,
                    help="ramos internos menores que isto sao colapsados (padrao 1e-5; 0 desliga)")
    ap.add_argument("--regex-id", default=ACC_RE)
    ap.add_argument("--mapa")
    ap.add_argument("--saida", default="resultados_por_gene")
    a = ap.parse_args()

    # ---------- alinhamentos
    arqs = []
    for pat in a.padroes.split(","):
        arqs += glob.glob(os.path.join(a.alinhamentos_dir, pat))
    arqs = sorted(set(arqs))
    if not arqs:
        sys.exit("nenhum alinhamento em " + a.alinhamentos_dir)
    tab = []
    for f in arqs:
        gene = os.path.splitext(os.path.basename(f))[0]
        try:
            regs_f = ler_fasta(f)
            if a.excluir:
                rx_al = re.compile(a.excluir)
                regs_f = [r for r in regs_f if not rx_al.search(r[0])]
            m = metricas_alinhamento(regs_f)
        except ValueError as e:
            print("AVISO: %s ignorado: %s" % (f, e), file=sys.stderr)
            continue
        m["gene"] = gene
        tab.append(m)
    cols = ["gene", "tipo", "n_seqs", "comprimento", "sitios_com_dados", "constantes",
            "variaveis", "pars_informativos", "prop_variaveis", "prop_pi", "pct_gap_ou_ambiguo"]
    with open(a.saida + "_variabilidade.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for m in tab:
            fh.write("\t".join(("%.4f" % m[c]) if isinstance(m[c], float) else str(m[c])
                               for c in cols) + "\n")
    print("variabilidade -> %s_variabilidade.tsv (%d genes)" % (a.saida, len(tab)))

    # ---------- arvores por gene
    if not a.ref:
        return
    mapa = carregar_mapa(a.mapa) if a.mapa else None
    norm = fabrica_normalizador(a.regex_id, mapa)
    ref = preparar_arvore(ler_arvore_arquivo(a.ref), norm, a.excluir)
    todos_ref = frozenset(folhas(ref))
    ref_splits, _ = obter_splits(ref, None, a.min_ramo)
    ref_rotulado = any(v is not None for v in ref_splits.values())
    if not ref_rotulado:
        print("AVISO: arvore de referencia sem rotulos de suporte; so RF bruto.", file=sys.stderr)

    genes_arv = {}
    if a.arvores_dir:
        for m in tab:
            p = os.path.join(a.arvores_dir, m["gene"] + a.sufixo_arvore)
            if os.path.exists(p):
                genes_arv[m["gene"]] = ler_arvore_arquivo(p)
    elif a.arvores_multi:
        if not a.nomes_genes:
            sys.exit("--arvores-multi exige --nomes-genes (ordem das particoes)")
        nomes = lista_arg(a.nomes_genes)
        with open(a.arvores_multi) as fh:
            linhas = [l for l in fh.read().split("\n") if l.strip()]
        if len(linhas) != len(nomes):
            sys.exit("%d arvores no arquivo mas %d nomes de genes" % (len(linhas), len(nomes)))
        for n, l in zip(nomes, linhas):
            genes_arv[n] = ler_newick(l)
    else:
        sys.exit("com --ref informe --arvores-dir ou --arvores-multi")

    faltam = [m["gene"] for m in tab if m["gene"] not in genes_arv]
    if faltam:
        print("AVISO: sem arvore para: " + ", ".join(faltam), file=sys.stderr)

    linhas_rf, colunas_matriz, presenca = [], [], {}
    ref_sup = filtrar_suportados(ref_splits, a.sh_min, a.ufb_min)
    for gene, raiz in genes_arv.items():
        raiz = preparar_arvore(raiz, norm, a.excluir)
        T = frozenset(folhas(raiz)) & todos_ref
        if len(T) < 4:
            print("AVISO: %s tem <4 taxons em comum com a referencia" % gene, file=sys.stderr)
            continue
        g_spl, _ = obter_splits(raiz, T, a.min_ramo)
        r_spl, _ = obter_splits(ref, T, a.min_ramo)
        N = len(T)
        bruto = rf(g_spl, r_spl, N)
        g_rot = any(v is not None for v in g_spl.values())
        if g_rot and ref_rotulado:
            col = rf(filtrar_suportados(g_spl, a.sh_min, a.ufb_min),
                     filtrar_suportados(r_spl, a.sh_min, a.ufb_min), N)
            col_rf, col_norm = col["rf"], "%.3f" % col["rf_norm"]
        else:
            col_rf, col_norm = "NA", "NA"
        linhas_rf.append([gene, N, len(g_spl), len(r_spl), bruto["compartilhados"],
                          bruto["rf"], bruto["rf_max"], "%.3f" % bruto["rf_norm"],
                          col_rf, col_norm, "sim" if g_rot else "nao"])
        # matriz
        colunas_matriz.append(gene)
        for i, (lado, vref) in enumerate(sorted(ref_splits.items(),
                                                key=lambda kv: (len(kv[0]), sorted(kv[0])))):
            sid = "S%02d" % (i + 1)
            sp = restringir_split(lado, todos_ref, T)
            if sp is None:
                cod = "."
            elif sp not in g_spl:
                cod = "-"
            else:
                pv = passa(g_spl[sp], a.sh_min, a.ufb_min) if g_rot else None
                cod = "?" if pv is None else ("P" if pv else "p")
            presenca.setdefault((sid, lado), {})[gene] = cod

    with open(a.saida + "_RF_por_gene.tsv", "w") as fh:
        fh.write("\t".join(["gene", "n_taxons", "splits_gene", "splits_ref", "compartilhados",
                            "RF", "RF_max", "RF_norm", "RF_colapsado", "RF_colapsado_norm",
                            "gene_tem_suporte"]) + "\n")
        for l in linhas_rf:
            fh.write("\t".join(map(str, l)) + "\n")
    with open(a.saida + "_matriz_splits.tsv", "w") as fh:
        fh.write("\t".join(["split", "n_taxons_lado", "ref_SH", "ref_UFB", "genes_com_split",
                            "genes_aplicaveis", "taxons_do_lado"] + colunas_matriz) + "\n")
        for (sid, lado), d in sorted(presenca.items(), key=lambda kv: kv[0][0]):
            v = ref_splits[lado]
            sh = "%g" % v[0] if v else "NA"
            ufb = "%g" % v[1] if v and len(v) > 1 else "NA"
            cods = [d.get(g, ".") for g in colunas_matriz]
            aplic = sum(1 for c in cods if c != ".")
            pres = sum(1 for c in cods if c in ("P", "p", "?"))
            fh.write("\t".join([sid, str(len(lado)), sh, ufb, str(pres), str(aplic),
                                ";".join(sorted(lado))] + cods) + "\n")
    print("RF por gene -> %s_RF_por_gene.tsv ; matriz -> %s_matriz_splits.tsv" % (a.saida, a.saida))
    print("Lembrete: RF de gene x concatenada e parcialmente circular (o gene faz parte da "
          "concatenacao). Use como descricao; a leitura principal e a matriz de splits.")


if __name__ == "__main__":
    main()
