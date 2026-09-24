#!/usr/bin/env python3
"""Recorta o alinhamento para a regiao coberta pelos fragmentos mais curtos
(ex: guaxinim, 96pb), preservando o alinhamento completo original intacto."""
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

def primeira_ultima_base(seq):
    primeira = next((i for i, c in enumerate(seq) if c != "-"), None)
    ultima = next((i for i in range(len(seq)-1, -1, -1) if seq[i] != "-"), None)
    return primeira, ultima

if __name__ == "__main__":
    infile = sys.argv[1]
    outfile = sys.argv[2]
    ref_ids = sys.argv[3:]  # ids das sequencias curtas que definem o nucleo

    seqs = ler_fasta(infile)

    inicio, fim = None, None
    for rid in ref_ids:
        matches = [h for h in seqs if h.startswith(rid)]
        if not matches:
            print(f"AVISO: {rid} nao encontrado no alinhamento", file=sys.stderr)
            continue
        h = matches[0]
        p, u = primeira_ultima_base(seqs[h])
        print(f"{h}: cobre colunas {p}-{u}")
        if inicio is None or p > inicio:
            inicio = p
        if fim is None or u < fim:
            fim = u

    print(f"\nRegiao de nucleo comum (intersecao): colunas {inicio}-{fim} ({fim-inicio+1} pb)")

    with open(outfile, "w") as out:
        for h, s in seqs.items():
            trecho = s[inicio:fim+1]
            out.write(f">{h}\n{trecho}\n")

    print(f"Nucleo comum salvo em: {outfile}")
