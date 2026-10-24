#!/usr/bin/env python3
# Traduz um FASTA de nucleotideos (assume quadro 0, codigo mitocondrial de invertebrados
# transl_table=5) para proteina, pronta para alinhar com MAFFT. Remove o stop final se houver.
# Se achar stop CODON INTERNO, avisa em stderr (sinal de possivel NUMT/erro de coordenada) e
# marca esse aminoacido como X, sem travar o script.
# Uso: traduzir_para_alinhar.py entrada.fasta > proteina.fasta
import sys
bases = "TCAG"
codons = [a + b + c for a in bases for b in bases for c in bases]
aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
TAB = dict(zip(codons, aas))
TAB["AGA"] = "S"; TAB["AGG"] = "S"; TAB["ATA"] = "M"; TAB["TGA"] = "W"
def ler_fasta(c):
    d, n = {}, None
    for l in open(c):
        l = l.strip()
        if l.startswith(">"): n = l[1:].split()[0]; d[n] = []
        elif n: d[n].append(l)
    return {k: "".join(v) for k, v in d.items()}
def traduzir(seq):
    seq = seq.upper()
    prot = [TAB.get(seq[i:i+3], "X") for i in range(0, len(seq) - 2, 3)]
    return "".join(prot)
seqs = ler_fasta(sys.argv[1])
for nome, s in seqs.items():
    p = traduzir(s)
    if p.endswith("*"): p = p[:-1]
    n_int = p.count("*")
    if n_int:
        sys.stderr.write("AVISO possivel NUMT: %s tem %d stop(s) codon interno(s)\n" % (nome, n_int))
        p = p.replace("*", "X")
    print(">%s\n%s" % (nome, p))
