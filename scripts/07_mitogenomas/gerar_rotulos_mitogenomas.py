#!/usr/bin/env python3
# Gera arquivo LABELS do iTOL: troca o nome exibido de cada ponta para "Accession Especie"
# (especie em italico Unicode), mesmo padrao usado nas arvores de COI e 12S.
# Uso: gerar_rotulos.py pasta_saida
import sys, os
def italico(s):
    out = []
    for c in s:
        if c == "h": out.append("\u210E")
        elif "A" <= c <= "Z": out.append(chr(0x1D434 + ord(c) - ord("A")))
        elif "a" <= c <= "z": out.append(chr(0x1D44E + ord(c) - ord("a")))
        else: out.append(c)
    return "".join(out)
ESPECIE = {
    "KX822021.1":  "Mansonella ozzardi",
    "OR271611.1":  "Mansonella ozzardi",
    "OR288092.1":  "Mansonella ozzardi",
    "OR296456.1":  "Mansonella ozzardi",
    "OR296457.1":  "Mansonella ozzardi",
    "OR296458.1":  "Mansonella ozzardi",
    "OR296459.1":  "Mansonella ozzardi",
    "MT361687.1":  "Mansonella perstans",
    "OQ633017.1":  "Mansonella perstans",
    "OQ633019.1":  "Mansonella perstans",
    "OQ633020.1":  "Mansonella sp. DEUX",
    "OQ633021.1":  "Mansonella sp. DEUX",
    "NC_004298.1": "Brugia malayi",
}
def rotulo(acc, esp):
    if esp == "Mansonella sp. DEUX":
        nome = italico("Mansonella") + ' sp. "DEUX"'
    else:
        nome = italico(esp)
    return "%s %s" % (acc, nome)
def main(out):
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "itol_rotulos_mitogenomas.txt"), "w", encoding="utf-8") as f:
        f.write("LABELS\nSEPARATOR TAB\nDATA\n")
        for acc, esp in ESPECIE.items():
            f.write("%s\t%s\n" % (acc, rotulo(acc, esp)))
    print("Gerado para %d accessions:" % len(ESPECIE))
    for acc, esp in ESPECIE.items():
        print("  %s -> %s" % (acc, rotulo(acc, esp)))
if __name__ == "__main__":
    if len(sys.argv) != 2: sys.exit("uso: gerar_rotulos.py pasta_saida")
    main(sys.argv[1])
