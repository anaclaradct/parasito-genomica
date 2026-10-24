#!/usr/bin/env python3
"""
podar_para_dambe.py - remove sequencias (por accession) de um alinhamento FASTA e,
opcionalmente, da arvore correspondente. Serve para a analise de sensibilidade do
teste de saturacao (Xia et al. 2003) no DAMBE sem as sequencias em quarentena
(ex.: NUMTs do Peru JF412338-JF412347).

O casamento e feito pelo accession SEM versao (JF412338.1 == JF412338).

Exemplos:
  python3 podar_para_dambe.py --fasta COI_alinhamento.fasta \
      --remover JF412338-JF412347 --saida COI_sem_NUMT.fasta \
      --arvore COI.treefile --saida-arvore COI_sem_NUMT.treefile

  --remover aceita: intervalos (JF412338-JF412347 ou JF412338-347), accessions soltos,
  listas separadas por virgula/espaco, ou @arquivo.txt (um por linha).
  --manter faz o inverso (fica so com os listados), util para restringir a taxons comuns.
"""
import argparse
import re
import sys
from filo_utils import (ler_fasta, escrever_fasta, acc_key, ler_arvore_arquivo,
                        podar_arvore, para_newick, folhas)


def expandir_ids(spec):
    itens = []
    for tok in re.split(r"[,\s]+", spec.strip()):
        if not tok:
            continue
        if tok.startswith("@"):
            with open(tok[1:]) as fh:
                for l in fh:
                    itens.extend(expandir_ids(l))
            continue
        itens.append(tok)
    ids = set()
    for tok in itens:
        tok = re.sub(r"\.\d+$", "", tok)
        m = re.fullmatch(r"([A-Za-z_]+)(\d+)-([A-Za-z_]*)(\d+)", tok)
        if m:
            p1, n1, p2, n2 = m.groups()
            p2 = p2 or p1
            if p1 != p2:
                sys.exit("intervalo com prefixos diferentes: " + tok)
            if len(n2) < len(n1):          # forma curta: JF412338-347 == JF412338-JF412347
                n2 = n1[:len(n1) - len(n2)] + n2
            a, b = int(n1), int(n2)
            if b < a:
                sys.exit("intervalo invertido: " + tok)
            larg = len(n1)
            for k in range(a, b + 1):
                ids.add("%s%0*d" % (p1, larg, k))
        else:
            ids.add(tok)
    return ids


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--saida", required=True)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--remover")
    g.add_argument("--manter")
    ap.add_argument("--arvore")
    ap.add_argument("--saida-arvore")
    ap.add_argument("--remover-colunas-vazias", action="store_true",
                    help="apaga colunas que ficam so com gaps depois da remocao")
    a = ap.parse_args()

    ids = expandir_ids(a.remover or a.manter)
    modo_remover = a.remover is not None
    regs = ler_fasta(a.fasta)
    if not regs:
        sys.exit("FASTA vazio: " + a.fasta)

    def alvo(nome):
        return acc_key(nome.split()[0]) in ids

    mantidos, removidos = [], []
    for cab, seq in regs:
        sai = alvo(cab) if modo_remover else not alvo(cab)
        (removidos if sai else mantidos).append((cab, seq))

    achados = {acc_key(c.split()[0]) for c, _ in regs}
    faltando = sorted(ids - achados)
    if faltando:
        print("AVISO: %d accession(s) pedido(s) nao encontrado(s) no FASTA (confira se o intervalo pedido tem lacunas legitimas): %s" %
              (len(faltando), ", ".join(faltando)), file=sys.stderr)

    if a.remover_colunas_vazias and mantidos:
        L = len(mantidos[0][1])
        vivas = [j for j in range(L) if any(s[j] not in "-?" for _, s in mantidos)]
        mantidos = [(c, "".join(s[j] for j in vivas)) for c, s in mantidos]

    escrever_fasta(mantidos, a.saida)
    print("FASTA: %d -> %d sequencias (removidas: %d)  -> %s" %
          (len(regs), len(mantidos), len(removidos), a.saida))
    if mantidos:
        print("comprimento do alinhamento de saida: %d" % len(mantidos[0][1]))

    if a.arvore:
        if not a.saida_arvore:
            sys.exit("--arvore exige --saida-arvore")
        raiz = ler_arvore_arquivo(a.arvore)
        antes = len(folhas(raiz))
        nova = podar_arvore(raiz, (lambda n: alvo(n)) if modo_remover
                            else (lambda n: not alvo(n)))
        with open(a.saida_arvore, "w") as fh:
            fh.write(para_newick(nova) + "\n")
        print("ARVORE: %d -> %d folhas -> %s" % (antes, len(folhas(nova)), a.saida_arvore))
        print("Obs.: os rotulos de suporte da arvore original sao mantidos, nao recalculados. "
              "Para a figura/relato final, o mais limpo e refazer o IQ-TREE no FASTA podado.")


if __name__ == "__main__":
    main()
