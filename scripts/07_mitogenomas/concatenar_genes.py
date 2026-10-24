#!/usr/bin/env python3
# Concatena varios alinhamentos de genes (mesmo comprimento dentro de cada arquivo) numa
# supermatriz, preenchendo com gap quem nao tem um gene, e gera arquivo de particao pro IQ-TREE.
# Uso: concatenar_genes.py saida_prefixo excluir(accA,accB ou vazio) gene1=arq1.fasta gene2=arq2.fasta ...
import sys
def ler_fasta(c):
    d, n = {}, None
    for l in open(c):
        l = l.strip()
        if l.startswith(">"): n = l[1:].split()[0]; d[n] = []
        elif n: d[n].append(l)
    return {k: "".join(v) for k, v in d.items()}
prefixo = sys.argv[1]
excluir = set(x for x in sys.argv[2].split(",") if x)
pares = [a.split("=", 1) for a in sys.argv[3:]]
genes = {}
for nome, arq in pares:
    d = ler_fasta(arq)
    for acc in excluir: d.pop(acc, None)
    comprimentos = set(len(s) for s in d.values())
    if len(comprimentos) != 1:
        sys.exit("ERRO: %s tem sequencias de tamanhos diferentes: %s" % (nome, comprimentos))
    genes[nome] = (d, comprimentos.pop())
todos_acc = sorted(set.union(*[set(d) for d, _ in genes.values()]))
if not todos_acc: sys.exit("nenhum accession sobrou apos exclusao")
with open(prefixo + ".fasta", "w") as fh, open(prefixo + "_particoes.txt", "w") as fp:
    pos = 1
    for nome, (d, L) in genes.items():
        fp.write("DNA, %s = %d-%d\n" % (nome, pos, pos + L - 1))
        pos += L
    for acc in todos_acc:
        fh.write(">%s\n" % acc)
        for nome, (d, L) in genes.items():
            fh.write(d.get(acc, "-" * L))
        fh.write("\n")
print("taxons: %d" % len(todos_acc))
print("genes: %d" % len(genes))
print("comprimento total: %d" % sum(L for _, L in genes.values()))
for nome, (d, L) in genes.items():
    print("  %s: %d pb, %d/%d taxons presentes" % (nome, L, len(d), len(todos_acc)))
