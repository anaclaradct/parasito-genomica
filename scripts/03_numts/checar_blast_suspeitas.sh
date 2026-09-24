#!/bin/bash
set -euo pipefail

REF_ACC="KX822021.1"
OUT_DIR="numt_check"
mkdir -p "$OUT_DIR"

if ! command -v blastn &> /dev/null; then
    echo "ERRO: blastn nao encontrado. Instale com: mamba install -c bioconda blast"
    exit 1
fi

echo ">> Baixando mitogenoma de referencia ($REF_ACC)..."
efetch -db nuccore -id "$REF_ACC" -format fasta > "$OUT_DIR/ref_mitogenoma_ozzardi.fasta"

echo ">> Extraindo sequencias suspeitas..."
python3 scripts/03_numts/extrair_suspeitas.py

echo ">> Rodando blastn..."
blastn -query "$OUT_DIR/suspeitas_blast.fasta" -subject "$OUT_DIR/ref_mitogenoma_ozzardi.fasta" \
    -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -out "$OUT_DIR/resultado_blast.tsv"

echo ""
cat "$OUT_DIR/resultado_blast.tsv"
