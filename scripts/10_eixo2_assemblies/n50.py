import sys
for f in sys.argv[1:]:
    t, L = 0, []
    for l in open(f):
        if l.startswith(">"):
            L.append(0)
        else:
            L[-1] += len(l.strip())
    L.sort(reverse=True)
    t = sum(L); a = 0
    for i, x in enumerate(L, 1):
        a += x
        if a >= t / 2:
            print(f"{f}\t{t}\t{len(L)}\t{x}\t{i}\t{L[0]}")
            break
