#!/usr/bin/env python3
"""
resolucao_arvores.py - "proporcao de nos internos fortemente suportados" e RF entre
arvores, SOMENTE no conjunto de taxons comuns.

Uso:
  python3 resolucao_arvores.py \
     --arvore mito=mtDNA.treefile --arvore rrnS_mito=rrnS_mito.treefile \
     --arvore cox1_mito=cox1_mito.treefile --arvore 12S_publico=12S.treefile \
     --ref mito --excluir "Brugia|NC_004298" --saida resolucao_taxons_comuns

- IDs sao normalizados para o accession sem versao/coordenadas (KX822021.1:9215-9879 ->
  KX822021). Se os nomes nao casarem, use --mapa (TSV nome_original<TAB>id_comum).
- O conjunto comum e a intersecao dos taxons de todas as arvores (ou --lista-taxons).
- ATENCAO: o suporte NAO e recalculado ao restringir a arvore. Para numeros finais, o
  mais limpo e podar os alinhamentos (podar_para_dambe.py --manter) e refazer o IQ-TREE.
- O COI publico (68 seqs) nao compartilha accessions com os mitogenomas; use o cox1
  extraido dos mitogenomas para entrar nesta comparacao.
"""
import argparse
import re
import sys
from filo_utils import (ler_arvore_arquivo, renomear_folhas, folhas, podar_arvore, obter_splits,
                        filtrar_suportados, rf, passa, fabrica_normalizador, carregar_mapa, ACC_RE)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arvore", action="append", required=True, metavar="ROTULO=ARQUIVO")
    ap.add_argument("--ref", help="rotulo da arvore de referencia para o RF")
    ap.add_argument("--lista-taxons", help="arquivo com IDs a manter (um por linha)")
    ap.add_argument("--excluir", help="regex de taxons a ignorar (outgroup)")
    ap.add_argument("--sh-min", type=float, default=80.0)
    ap.add_argument("--ufb-min", type=float, default=95.0)
    ap.add_argument("--min-ramo", type=float, default=1e-5,
                    help="ramos internos menores que isto sao colapsados (padrao 1e-5; 0 desliga)")
    ap.add_argument("--regex-id", default=ACC_RE)
    ap.add_argument("--mapa")
    ap.add_argument("--saida", default="resolucao_taxons_comuns")
    a = ap.parse_args()

    mapa = carregar_mapa(a.mapa) if a.mapa else None
    norm = fabrica_normalizador(a.regex_id, mapa)
    rx = re.compile(a.excluir) if a.excluir else None

    arvores = {}
    for item in a.arvore:
        if "=" not in item:
            sys.exit("use ROTULO=ARQUIVO: " + item)
        rot, arq = item.split("=", 1)
        raiz = ler_arvore_arquivo(arq)
        if rx:  # exclui pelo nome ORIGINAL ou normalizado
            raiz = podar_arvore(raiz, lambda n: bool(rx.search(n)) or bool(rx.search(norm(n))))
        renomear_folhas(raiz, norm)
        arvores[rot] = raiz

    conjuntos = {r: set(folhas(t)) for r, t in arvores.items()}
    if a.lista_taxons:
        with open(a.lista_taxons) as fh:
            comum = {norm(l.strip()) for l in fh if l.strip()}
        comum &= set.intersection(*conjuntos.values())
    else:
        comum = set.intersection(*conjuntos.values())
    print("taxons por arvore: " + "; ".join("%s=%d" % (r, len(c)) for r, c in conjuntos.items()))
    print("taxons comuns: %d" % len(comum))
    if len(comum) < 4:
        sys.exit("menos de 4 taxons comuns - nada a comparar (confira IDs / use --mapa)")
    comum = frozenset(comum)
    N = len(comum)

    spl = {r: obter_splits(t, comum, a.min_ramo)[0] for r, t in arvores.items()}
    linhas = []
    for r, s in spl.items():
        rotulado = any(v is not None for v in s.values())
        sh_ok = sum(1 for v in s.values() if passa(v, a.sh_min, None) is True)
        ufb_ok = sum(1 for v in s.values() if passa(v, None, a.ufb_min) is True)
        ambos = sum(1 for v in s.values() if passa(v, a.sh_min, a.ufb_min) is True)
        sem_info = sum(1 for v in s.values() if v is None)
        lin = {"arvore": r, "N_taxons": N, "nos_internos": len(s), "max_possivel": N - 3,
               "SH_ok": sh_ok, "UFB_ok": ufb_ok, "ambos_ok": ambos, "sem_rotulo": sem_info,
               "prop_ambos_sobre_existentes": ("%.3f" % (ambos / len(s))) if s and rotulado else "NA",
               "prop_ambos_sobre_maximo": ("%.3f" % (ambos / (N - 3))) if rotulado else "NA"}
        if a.ref:
            if a.ref not in spl:
                sys.exit("--ref nao e um dos rotulos: " + a.ref)
            b = rf(s, spl[a.ref], N)
            lin.update({"RF_vs_ref": b["rf"], "RF_max": b["rf_max"], "RF_norm": "%.3f" % b["rf_norm"]})
            ref_rot = any(v is not None for v in spl[a.ref].values())
            if rotulado and ref_rot:
                c = rf(filtrar_suportados(s, a.sh_min, a.ufb_min),
                       filtrar_suportados(spl[a.ref], a.sh_min, a.ufb_min), N)
                lin.update({"RF_colapsado": c["rf"], "RF_colapsado_norm": "%.3f" % c["rf_norm"]})
            else:
                lin.update({"RF_colapsado": "NA", "RF_colapsado_norm": "NA"})
        linhas.append(lin)

    cols = list(linhas[0].keys())
    with open(a.saida + ".tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for l in linhas:
            fh.write("\t".join(str(l.get(c, "")) for c in cols) + "\n")
    for l in linhas:
        print("%-14s nos=%2d  suportados(SH>=%g & UFB>=%g)=%2d  prop/existentes=%s" %
              (l["arvore"], l["nos_internos"], a.sh_min, a.ufb_min, l["ambos_ok"],
               l["prop_ambos_sobre_existentes"]))
    print("-> %s.tsv" % a.saida)
    print("Suporte nao recalculado apos restringir aos taxons comuns; nos colapsados usam o menor suporte.")


if __name__ == "__main__":
    main()
