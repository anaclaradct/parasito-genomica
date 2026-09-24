import sys
import csv

def parse_fasta(path):
    records = {}
    header, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(seq).upper()
                header = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if header is not None:
            records[header] = "".join(seq).upper()
    return records

fasta_path = sys.argv[1]
acc1, acc2 = sys.argv[2], sys.argv[3]

recs = parse_fasta(fasta_path)
s1, s2 = recs[acc1], recs[acc2]
diffs = [(i+1, a, b) for i, (a, b) in enumerate(zip(s1, s2)) if a != b and a != '-' and b != '-']
print(f"{acc1} vs {acc2}: {len(diffs)} diferencas (de {min(len(s1),len(s2))} pb comparados)")
for pos, a, b in diffs[:20]:
    print(f"  posicao {pos}: {a} -> {b}")
if len(diffs) > 20:
    print(f"  ... e mais {len(diffs)-20}")
