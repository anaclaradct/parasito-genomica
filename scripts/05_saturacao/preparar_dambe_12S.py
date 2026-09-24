#!/usr/bin/env python3
"""Converte o alinhamento FASTA (686pb, sem outgroup) para formato PHYLIP relaxado,
compativel com DAMBE, para o teste de saturacao de Xia et al. (2003)."""
import sys

def ler_fasta(path):
    seqs = {}
    header = None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                header = line[1:].split()[0]
                seqs[header] = []
            else:
                seqs[header].append(line)
    return {h: "".join(s) for h, s in seqs.items()}

if __name__ == "__main__":
    infile, outfile = sys.argv[1], sys.argv[2]
    seqs = ler_fasta(infile)
    n = len(seqs)
    comp = len(next(iter(seqs.values())))

    with open(outfile, "w") as out:
        out.write(f" {n} {comp}\n")
        for h, s in seqs.items():
            out.write(f"{h}  {s}\n")

    print(f"PHYLIP relaxado salvo: {n} sequencias, {comp} pb -> {outfile}")
