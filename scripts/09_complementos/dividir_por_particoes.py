#!/usr/bin/env python3
"""
dividir_por_particoes.py - corta uma matriz concatenada (FASTA) em um alinhamento por
gene, usando um arquivo de particoes no formato do IQ-TREE/RAxML:
    DNA, atp6 = 1-672
    DNA, cob = 673-1773
(coordenadas 1-based, inclusivas). Garante que os genes correspondem EXATAMENTE aos
locos usados nas arvores por locus (-S) dessa mesma matriz.

Exemplo:
  python3 dividir_por_particoes.py \
     --fasta 04_alinhamentos_genes/matriz_principal.fasta \
     --particoes 04_alinhamentos_genes/matriz_principal_particoes.txt \
     --saida-dir genes_da_matriz_principal
"""
import argparse
import os
import re
import sys
from filo_utils import ler_fasta, escrever_fasta

LINHA = re.compile(r"^\s*[A-Za-z0-9_]+\s*,\s*(\S+)\s*=\s*(\d+)\s*-\s*(\d+)\s*;?\s*$")


def ler_particoes(caminho):
    partes = []
    with open(caminho) as fh:
        for l in fh:
            if not l.strip() or l.lstrip().startswith("#"):
                continue
            m = LINHA.match(l)
            if not m:
                print("AVISO: linha de particao ignorada: " + l.strip(), file=sys.stderr)
                continue
            partes.append((m.group(1), int(m.group(2)), int(m.group(3))))
    return partes


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--particoes", required=True)
    ap.add_argument("--saida-dir", required=True)
    a = ap.parse_args()

    regs = ler_fasta(a.fasta)
    if not regs:
        sys.exit("FASTA vazio")
    L = {len(s) for _, s in regs}
    if len(L) != 1:
        sys.exit("sequencias com comprimentos diferentes - nao e uma matriz alinhada")
    L = L.pop()
    partes = ler_particoes(a.particoes)
    if not partes:
        sys.exit("nenhuma particao lida")
    os.makedirs(a.saida_dir, exist_ok=True)
    print("matriz: %d sequencias x %d colunas; %d particoes" % (len(regs), L, len(partes)))

    anterior_fim = 0
    for nome, ini, fim in partes:
        if fim > L or ini < 1 or fim < ini:
            sys.exit("particao %s (%d-%d) fora da matriz de %d colunas" % (nome, ini, fim, L))
        if ini != anterior_fim + 1:
            print("AVISO: %s comeca em %d mas a anterior terminou em %d (lacuna/sobreposicao)" %
                  (nome, ini, anterior_fim), file=sys.stderr)
        anterior_fim = fim
        saida = [(cab, seq[ini - 1:fim]) for cab, seq in regs]
        escrever_fasta(saida, os.path.join(a.saida_dir, nome + ".fasta"))
        print("  %-8s %5d-%-5d  %5d colunas" % (nome, ini, fim, fim - ini + 1))
    if anterior_fim != L:
        print("AVISO: particoes terminam em %d mas a matriz tem %d colunas" % (anterior_fim, L),
              file=sys.stderr)
    print("-> " + a.saida_dir)


if __name__ == "__main__":
    main()
