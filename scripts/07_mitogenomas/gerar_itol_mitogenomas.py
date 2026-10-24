#!/usr/bin/env python3
# Gera datasets do iTOL (especie/cor + localidade) para a arvore de mitogenomas (13 taxons:
# 12 Mansonella + Brugia malayi como outgroup). Nomes cientificos em italico Unicode.
# Uso: gerar_itol_mitogenomas.py pasta_saida
import sys, os
def italico(s):
    out = []
    for c in s:
        if c == "h": out.append("\u210E")
        elif "A" <= c <= "Z": out.append(chr(0x1D434 + ord(c) - ord("A")))
        elif "a" <= c <= "z": out.append(chr(0x1D44E + ord(c) - ord("a")))
        else: out.append(c)
    return "".join(out)
TAB = {
    "KX822021.1":  ("Mansonella ozzardi", "Brazil"),
    "OR271611.1":  ("Mansonella ozzardi", "Colombia"),
    "OR288092.1":  ("Mansonella ozzardi", "Colombia"),
    "OR296456.1":  ("Mansonella ozzardi", "Colombia"),
    "OR296457.1":  ("Mansonella ozzardi", "Colombia"),
    "OR296458.1":  ("Mansonella ozzardi", "Brazil"),
    "OR296459.1":  ("Mansonella ozzardi", "Venezuela"),
    "MT361687.1":  ("Mansonella perstans", "Cameroon"),
    "OQ633017.1":  ("Mansonella perstans", "Gabon"),
    "OQ633019.1":  ("Mansonella perstans", "Gabon"),
    "OQ633020.1":  ("Mansonella sp. DEUX", "Gabon"),
    "OQ633021.1":  ("Mansonella sp. DEUX", "Gabon"),
    "NC_004298.1": ("Brugia malayi outgroup", "-"),
}
CORES = {
    "Mansonella perstans": "#e6194B",
    "Mansonella ozzardi": "#3cb44b",
    "Mansonella sp. DEUX": "#4363d8",
    "Brugia malayi outgroup": "#999999",
}
LEGENDA = {
    "Mansonella perstans": italico("Mansonella perstans"),
    "Mansonella ozzardi": italico("Mansonella ozzardi"),
    "Mansonella sp. DEUX": italico("Mansonella") + ' sp. "DEUX"',
    "Brugia malayi outgroup": italico("Brugia malayi") + " (outgroup)",
}
def main(out):
    os.makedirs(out, exist_ok=True)
    ordem = ["Mansonella perstans", "Mansonella ozzardi", "Mansonella sp. DEUX", "Brugia malayi outgroup"]
    with open(os.path.join(out, "itol_especie_mitogenomas.txt"), "w", encoding="utf-8") as f:
        f.write("DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\tEspecie\nCOLOR\t#000000\n")
        f.write("LEGEND_TITLE\tEspecie\n")
        f.write("LEGEND_SHAPES\t" + "\t".join("1" for _ in ordem) + "\n")
        f.write("LEGEND_COLORS\t" + "\t".join(CORES[e] for e in ordem) + "\n")
        f.write("LEGEND_LABELS\t" + "\t".join(LEGENDA[e] for e in ordem) + "\n")
        f.write("DATA\n")
        for acc, (esp, _) in TAB.items():
            f.write("%s\t%s\n" % (acc, CORES[esp]))
    with open(os.path.join(out, "itol_localidade_mitogenomas.txt"), "w", encoding="utf-8") as f:
        f.write("DATASET_TEXT\nSEPARATOR TAB\nDATASET_LABEL\tLocalidade\nCOLOR\t#000000\nDATA\n")
        for acc, (_, pais) in TAB.items():
            f.write("%s\t%s\t-1\t#000000\tnormal\t1\t0\n" % (acc, pais))
    print("Gerado para %d accessions:" % len(TAB))
    for acc, (esp, pais) in TAB.items():
        print("  %s\t%s\t%s" % (acc, LEGENDA[esp], pais))
if __name__ == "__main__":
    if len(sys.argv) != 2: sys.exit("uso: gerar_itol_mitogenomas.py pasta_saida")
    main(sys.argv[1])
