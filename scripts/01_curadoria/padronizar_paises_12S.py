import csv, sys, re, os
arq = sys.argv[1]
PT = {"Brasil": "Brazil", "Camarões": "Cameroon", "Gabão": "Gabon", "Japão": "Japan",
      "Malásia": "Malaysia", "Colômbia": "Colombia", "Guiana Francesa": "French Guiana"}
INFERIDO = {"LT623914": "Brazil", "KP760340": "Haiti", "OR271611": "Colombia"}
linhas = list(csv.reader(open(arq, encoding="utf-8"), delimiter="\t"))
cab, dados = linhas[0], linhas[1:]
n, ip = len(cab), cab.index("pais")
tem_nota = "pais_nota" in cab
if not tem_nota:
    cab.append("pais_nota")
for r in dados:
    while len(r) > n and "" in r[ip + 1:]:
        del r[r.index("", ip + 1)]
    if len(r) != n:
        sys.exit(f"linha com colunas erradas, nada foi gravado: {r}")
    acc = r[0].split(".")[0]
    p = PT.get(re.split(r"[:(]", r[ip])[0].strip(), re.split(r"[:(]", r[ip])[0].strip())
    nota = ""
    if acc in INFERIDO:
        p, nota = INFERIDO[acc], "inferido pela publicacao/isolado associado"
    r[ip] = p if p and p not in ("NA", "Sem_dado") else "Unknown"
    if not tem_nota:
        r.append(nota)
tmp = arq + ".tmp"
csv.writer(open(tmp, "w", encoding="utf-8", newline=""), delimiter="\t", lineterminator="\n").writerows([cab] + dados)
os.replace(tmp, arq)
print("ok:", len(dados), "linhas")
