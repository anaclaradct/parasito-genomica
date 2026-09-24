#!/usr/bin/env python3
"""Remove sequencias especificas do alinhamento (por ID), sem alterar as colunas."""
import sys

def ler_fasta(path):
    seqs = {}
    header = None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                header = line[1:]
                seqs[header] = []
            else:
                seqs[header].append(line)
    return {h: "".join(s) for h, s in seqs.items()}

if __name__ == "__main__":
    infile, outfile = sys.argv[1], sys.argv[2]
    ids_remover = sys.argv[3:]

    seqs = ler_fasta(infile)
    removidas = []
    mantidas = {}
    for h, s in seqs.items():
        if any(h.startswith(rid) for rid in ids_remover):
            removidas.append(h)
        else:
            mantidas[h] = s

    print(f"Sequencias no arquivo original: {len(seqs)}")
    print(f"Removidas ({len(removidas)}):")
    for r in removidas:
        print(f"  - {r}")
    print(f"Mantidas: {len(mantidas)}")

    with open(outfile, "w") as out:
        for h, s in mantidas.items():
            out.write(f">{h}\n{s}\n")

    print(f"Salvo em: {outfile}")
