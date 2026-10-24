import sys, csv
def fasta(p):
    r, h = {}, None
    for l in open(p):
        l = l.strip()
        if l.startswith(">"): h = l[1:].split()[0]; r[h] = ""
        elif h: r[h] += l.upper()
    return r
seqs = fasta(sys.argv[1])
sp = {r["accession"]: r["especie_pasta"] for r in csv.DictReader(open(sys.argv[2]), delimiter="\t")}
acc = [a for a in seqs if a in sp]
res = {}
for i, a in enumerate(acc):
    for b in acc[i+1:]:
        if sp[a] == sp[b]: continue
        n = d = 0
        for x, y in zip(seqs[a], seqs[b]):
            if x in "ACGT" and y in "ACGT":
                n += 1; d += x != y
        if n == 0: continue
        for s in (sp[a], sp[b]):
            res.setdefault(s, []).append((d / n, n, a, b))
for s in sorted(res):
    print("==", s)
    for p, n, a, b in sorted(res[s])[:3]:
        print(f"  p={p:.4f}  sitios={n:4d}  {a} ({sp[a]}) x {b} ({sp[b]})")
