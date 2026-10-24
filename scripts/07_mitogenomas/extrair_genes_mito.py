#!/usr/bin/env python3
"""Extrai genes mitocondriais (12 CDS + 2 rRNA) de um arquivo GenBank multi-registro.
Uso: python3 extrair_genes_mito.py mitogenomas.gb pasta_saida
Saidas: <gene>.fasta (header >ACCESSION) e presenca_genes.tsv
Sem dependencias externas. Trata complement(), join(), <, >.
"""
import sys, os, re

GENES = ["atp6", "cob", "cox1", "cox2", "cox3", "nad1", "nad2", "nad3",
         "nad4", "nad4l", "nad5", "nad6", "rrnS", "rrnL"]
SINONIMOS = {
    "atp6": ["atp6", "atpase6", "atp synthase f0 subunit 6", "atp synthase subunit 6"],
    "cob": ["cob", "cytb", "cyt b", "cytochrome b"],
    "cox1": ["cox1", "coi", "co1", "coxi", "cytochrome c oxidase subunit i"],
    "cox2": ["cox2", "coii", "co2", "coxii", "cytochrome c oxidase subunit ii"],
    "cox3": ["cox3", "coiii", "co3", "coxiii", "cytochrome c oxidase subunit iii"],
    "nad1": ["nad1", "nd1", "nadh dehydrogenase subunit 1"],
    "nad2": ["nad2", "nd2", "nadh dehydrogenase subunit 2"],
    "nad3": ["nad3", "nd3", "nadh dehydrogenase subunit 3"],
    "nad4": ["nad4", "nd4", "nadh dehydrogenase subunit 4"],
    "nad4l": ["nad4l", "nd4l", "nadh dehydrogenase subunit 4l"],
    "nad5": ["nad5", "nd5", "nadh dehydrogenase subunit 5"],
    "nad6": ["nad6", "nd6", "nadh dehydrogenase subunit 6"],
    "rrnS": ["rrns", "12s ribosomal rna", "small subunit ribosomal rna", "s-rrna", "12s rrna"],
    "rrnL": ["rrnl", "16s ribosomal rna", "large subunit ribosomal rna", "l-rrna", "16s rrna"],
}
MAPA = {s: g for g, lista in SINONIMOS.items() for s in lista}
COMP = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")

def revcomp(s):
    return s.translate(COMP)[::-1]

def extrair_local(loc, seq):
    loc = loc.replace("<", "").replace(">", "")
    if loc.startswith("complement(") and loc.endswith(")"):
        return revcomp(extrair_local(loc[11:-1], seq))
    if loc.startswith("join(") and loc.endswith(")"):
        partes, nivel, atual = [], 0, ""
        for c in loc[5:-1]:
            if c == "," and nivel == 0:
                partes.append(atual); atual = ""
                continue
            nivel += (c == "(") - (c == ")")
            atual += c
        partes.append(atual)
        return "".join(extrair_local(p.strip(), seq) for p in partes)
    m = re.fullmatch(r"(\d+)\.\.(\d+)", loc)
    if m:
        a, b = int(m.group(1)), int(m.group(2))
        return seq[a - 1:b]
    m = re.fullmatch(r"(\d+)", loc)
    if m:
        return seq[int(m.group(1)) - 1]
    raise ValueError("local nao suportado: " + loc)

def ler_registros(caminho):
    with open(caminho) as f:
        texto = f.read()
    for bloco in texto.split("\n//"):
        if "LOCUS" not in bloco:
            continue
        yield bloco[bloco.index("LOCUS"):]

def parse(bloco):
    acc = re.search(r"^VERSION\s+(\S+)", bloco, re.M).group(1)
    org = re.search(r'/organism="([^"]+)"', bloco)
    org = org.group(1) if org else "?"
    orig = bloco.split("\nORIGIN", 1)[1] if "\nORIGIN" in bloco else ""
    seq = re.sub(r"[^A-Za-z]", "", orig).upper()
    feats_txt = bloco.split("\nFEATURES", 1)[1].split("\nORIGIN", 1)[0] if "\nFEATURES" in bloco else ""
    feats, atual = [], None
    for linha in feats_txt.split("\n")[1:]:
        if len(linha) > 5 and linha[5] != " ":
            if atual: feats.append(atual)
            tipo, loc = linha[5:21].strip(), linha[21:].strip()
            atual = {"tipo": tipo, "loc": loc, "quals": []}
        elif atual is not None:
            t = linha.strip()
            if t.startswith("/"):
                atual["quals"].append(t)
            elif atual["quals"]:
                atual["quals"][-1] += (" " if not atual["quals"][-1].startswith("/translation") else "") + t
            else:
                atual["loc"] += t  # local quebrado em varias linhas
    if atual: feats.append(atual)
    return acc, org, seq, feats

def qual(feat, nome):
    for q in feat["quals"]:
        if q.startswith("/" + nome + "="):
            return q.split("=", 1)[1].strip('"')
    return None

def identificar(feat):
    for campo in ("gene", "product"):
        v = qual(feat, campo)
        if v and v.lower().strip() in MAPA:
            return MAPA[v.lower().strip()]
    return None

def main(gb, saida):
    os.makedirs(saida, exist_ok=True)
    seqs = {g: [] for g in GENES}
    tabela = []
    for bloco in ler_registros(gb):
        acc, org, seq, feats = parse(bloco)
        achados, avisos = {}, []
        for ft in feats:
            if ft["tipo"] not in ("CDS", "rRNA"):
                continue
            g = identificar(ft)
            if g is None:
                continue
            if g in GENES[:12] and ft["tipo"] != "CDS":
                continue
            if g in achados:
                avisos.append(g + ":duplicado")
                continue
            s = extrair_local(ft["loc"], seq)
            achados[g] = (s, ft["loc"])
        for g, (s, loc) in achados.items():
            seqs[g].append((acc, org, loc, s))
        tabela.append([acc, org, str(len(seq)), str(len(achados))] +
                      [str(len(achados[g][0])) if g in achados else "0" for g in GENES] +
                      [";".join(avisos) or "-"])
    for g in GENES:
        with open(os.path.join(saida, g + ".fasta"), "w") as f:
            for acc, org, loc, s in seqs[g]:
                f.write(">%s\n%s\n" % (acc, s))
    with open(os.path.join(saida, "presenca_genes.tsv"), "w") as f:
        f.write("\t".join(["accession", "organismo", "tamanho", "n_genes"] + GENES + ["avisos"]) + "\n")
        for linha in tabela:
            f.write("\t".join(linha) + "\n")
    print("accession\torganismo\tn_genes(de 14)\tfaltando")
    for linha in tabela:
        falt = [g for g, v in zip(GENES, linha[4:18]) if v == "0"]
        print("%s\t%s\t%s\t%s" % (linha[0], linha[1], linha[3], ",".join(falt) or "-"))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
