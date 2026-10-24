#!/usr/bin/env python3
# Retrotraduz um alinhamento de proteina para nucleotideo, usando a sequencia de nucleotideos
# original (nao alinhada) de cada registro para recuperar os codons certos.
# Uso: retrotraduzir.py proteina_alinhada.fasta nucleotideo_original.fasta > nt_alinhado.fasta
import sys
def ler_fasta(c):
    d, n = {}, None
    for l in open(c):
        l = l.strip()
        if l.startswith(">"): n = l[1:].split()[0]; d[n] = []
        elif n: d[n].append(l)
    return {k: "".join(v) for k, v in d.items()}
prot_aln = ler_fasta(sys.argv[1]); nt_orig = ler_fasta(sys.argv[2])
falt = [n for n in prot_aln if n not in nt_orig]
if falt: sys.exit("sem nucleotideo original para: " + ",".join(falt))
for nome, p in prot_aln.items():
    nt = nt_orig[nome]
    pos, saida = 0, []
    for aa in p:
        if aa == "-":
            saida.append("---")
        else:
            saida.append(nt[pos:pos+3])
            pos += 3
    print(">%s\n%s" % (nome, "".join(saida)))
