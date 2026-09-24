#!/bin/bash
set -euo pipefail

ALIGNED="03_ALIGNMENTS/COI_aligned.fasta"
OUTGROUP_SRC_ACC="AF015193"
OUTGROUP_CDS_FASTA="03_ALIGNMENTS/Ovolvulus_all_cds.fasta"
OUTGROUP_FASTA="03_ALIGNMENTS/outgroup_Ovolvulus_COI.fasta"
ALIGNED_ROOTED="03_ALIGNMENTS/COI_aligned_rooted_Ovolvulus.fasta"
OUT_PREFIX="04_TREES/COI_rooted_Ovolvulus"

mkdir -p 04_TREES

if [ ! -f "$ALIGNED" ]; then
    echo "ERRO: nao encontrei $ALIGNED. Rode scripts/02_alinhamento/01_alinhar_COI.sh primeiro." >&2
    exit 1
fi

echo ">> Baixando todas as CDS do mitogenoma de referencia de O. volvulus ($OUTGROUP_SRC_ACC)..."
efetch -db nuccore -id "$OUTGROUP_SRC_ACC" -format fasta_cds_na > "$OUTGROUP_CDS_FASTA"

if [ ! -s "$OUTGROUP_CDS_FASTA" ]; then
    echo "ERRO: download do mitogenoma de referencia falhou ou veio vazio." >&2
    exit 1
fi

echo ">> Extraindo so o gene COI/cox1 desse conjunto..."
awk '
    BEGIN { keep = 0 }
    /^>/ {
        keep = (tolower($0) ~ /cox1|subunit 1\]|coi\b/) ? 1 : 0
    }
    keep { print }
' "$OUTGROUP_CDS_FASTA" > "$OUTGROUP_FASTA"

if [ ! -s "$OUTGROUP_FASTA" ]; then
    echo "ERRO: nao achei um gene com nome compativel a COI/cox1 nas CDS baixadas."
    echo "Cabecalhos disponiveis:"
    grep "^>" "$OUTGROUP_CDS_FASTA"
    exit 1
fi
echo "   OK: $(grep '^>' "$OUTGROUP_FASTA")"

echo ">> Adicionando outgroup ao alinhamento existente (mafft --add)..."
mafft --add "$OUTGROUP_FASTA" --keeplength --thread -1 "$ALIGNED" > "$ALIGNED_ROOTED" 2> "$OUT_PREFIX.mafft_add.log"

n_seqs=$(grep -c "^>" "$ALIGNED_ROOTED")
echo "   Alinhamento agora tem $n_seqs sequencias."

OUTGROUP_HEADER=$(grep "^>" "$ALIGNED_ROOTED" | tail -1 | sed 's/^>//' | awk '{print $1}')

IQTREE_BIN=""
for candidate in iqtree3 iqtree2 iqtree; do
    if command -v "$candidate" &> /dev/null; then
        IQTREE_BIN="$candidate"
        break
    fi
done
if [ -z "$IQTREE_BIN" ]; then
    echo "ERRO: nenhum binario iqtree encontrado no PATH." >&2
    exit 1
fi

echo ">> Rodando $IQTREE_BIN com outgroup explicito ($OUTGROUP_HEADER)..."
"$IQTREE_BIN" -s "$ALIGNED_ROOTED" -m MFP -alrt 1000 -bb 1000 -nt AUTO -o "$OUTGROUP_HEADER" -pre "$OUT_PREFIX"

echo ""
echo ">> Pronto. Arvore em: ${OUT_PREFIX}.treefile"
echo "   Compare a posicao de AM749265.1 (M. perforata) com a arvore anterior"
echo "   (04_TREES/COI_rooted.treefile, com outgroup Brugia malayi)."
