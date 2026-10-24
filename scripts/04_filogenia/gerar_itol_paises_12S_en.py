import csv, re, sys
tree, meta, out = sys.argv[1:4]
CORES = {"Gabon": "#FF4500", "Cameroon": "#9400D3", "Peru": "#FFD700",
         "Colombia": "#00CED1", "Malaysia": "#1E90FF", "Japan": "#8B4513",
         "Brazil": "#32CD32", "Venezuela": "#FF69B4", "Costa Rica": "#DC143C",
         "Togo": "#808000", "Haiti": "#00008B", "Unknown": "#CCCCCC"}
pais = {r["accession"].split(".")[0]: r["pais"]
        for r in csv.DictReader(open(meta, encoding="utf-8"), delimiter="\t")}
folhas = re.findall(r"[(,]([^(),:;]+(?::\d+-\d+)?)", open(tree).read())
linhas, usados, faltando = [], [], []
for f in folhas:
    acc = re.split(r"[._:]", f.strip("'"))[0]
    if acc not in pais:
        faltando.append(f); continue
    p = pais[acc] if pais[acc] in CORES else "Unknown"
    linhas.append(f"{f}\t{CORES[p]}\t{p}")
    if p not in usados: usados.append(p)
ordem = [p for p in CORES if p in usados]
with open(out, "w") as o:
    o.write("DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\tCountry of origin\nCOLOR\t#000000\n\n")
    o.write("LEGEND_TITLE\tCountry of origin\n")
    o.write("LEGEND_SHAPES\t" + "\t".join("1" for _ in ordem) + "\n")
    o.write("LEGEND_COLORS\t" + "\t".join(CORES[p] for p in ordem) + "\n")
    o.write("LEGEND_LABELS\t" + "\t".join(ordem) + "\n\nDATA\n")
    o.write("\n".join(linhas) + "\n")
print(len(linhas), "folhas com pais;", "sem metadado:", faltando)
