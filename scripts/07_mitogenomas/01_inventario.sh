#!/usr/bin/env bash
# Inventario: mitogenomas (>5 kb) de Mansonella + outgroups Onchocercidae (RefSeq)
set -uo pipefail
export PATH="$HOME/edirect:$PATH"
OUT="$HOME/Mansonella_projeto/07_MITOGENOMAS/00_inventario"
CAMPOS="-element AccessionVersion TaxId Slen Genome Completeness CreateDate Title SubType SubName"

esearch -db nuccore -query '"Mansonella"[Organism] AND mitochondrion[filter] AND 5000:20000[SLEN]' \
  | esummary | xtract -pattern DocumentSummary $CAMPOS \
  | sort -t$'\t' -k1,1 > "$OUT/mansonella_mitogenomas.tsv"

esearch -db nuccore -query '"Onchocercidae"[Organism] NOT "Mansonella"[Organism] AND mitochondrion[filter] AND complete genome[Title] AND refseq[filter]' \
  | esummary | xtract -pattern DocumentSummary $CAMPOS \
  | sort -t$'\t' -k7,7 > "$OUT/outgroups_onchocercidae_refseq.tsv"

wc -l "$OUT"/*.tsv
echo "=== Mansonella (accession / tamanho / titulo) ==="
cut -f1,3,7 "$OUT/mansonella_mitogenomas.tsv"
echo "=== Outgroups ==="
cut -f1,3,7 "$OUT/outgroups_onchocercidae_refseq.tsv"
