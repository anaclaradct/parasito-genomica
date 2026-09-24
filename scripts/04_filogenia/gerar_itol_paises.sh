#!/bin/bash
# Extrai accessions do treefile, busca pais no GenBank, gera dataset iTOL (colorstrip)

TREEFILE="04_TREES/12S_final/12S_com_outgroup_brugia.treefile"
OUTFILE="04_TREES/12S_final/itol_paises_dataset.txt"

# extrair accessions unicos da arvore (remove sufixo de coordenada tipo _9215-9879)
grep -oE '[A-Z]{2}[0-9]+\.[0-9]+(_[0-9]+-[0-9]+)?' "$TREEFILE" | sed 's/_[0-9]*-[0-9]*$//' | sed 's/\.[0-9]*$//' | sort -u > /tmp/tips_accessions.txt
wc -l /tmp/tips_accessions.txt

efetch -db nuccore -input /tmp/tips_accessions.txt -format gb > /tmp/tips_gb.gb

# extrair pais por accession (pega so a primeira palavra antes de ':' quando houver)
awk '
  /^ACCESSION/{acc=$2}
  /\/geo_loc_name="/{
    gsub(/.*\/geo_loc_name="|"/,"",$0)
    split($0, a, ":")
    print acc"\t"a[1]
  }
' /tmp/tips_gb.gb > /tmp/acc_pais.tsv

wc -l /tmp/acc_pais.tsv
