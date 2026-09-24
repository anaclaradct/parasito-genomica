#!/bin/bash
set -euo pipefail

ALIGNED="03_ALIGNMENTS/COI_aligned.fasta"
OUTGROUP_ACC="MN598543.1"
OUTGROUP_FASTA="03_ALIGNMENTS/outgroup_Bmalayi_COI.fasta"
ALIGNED_ROOTED="03_ALIGNMENTS/COI_aligned_rooted.fasta"
OUT_PREFIX="04_TREES/COI_rooted"

mkdir -p 04_TREES

if [ ! -f "$ALIGNED" ]; then
    echo "ERRO: nao encontrei $ALIGNED. Rode scripts/02_alinhamento/01_alinhar_COI.sh primeiro." >&2
    exit 1
fi

echo ">> Baixando outgroup $OUTGROUP_ACC (Brugia malayi COI)..."
efetch -db nuccore -id "$OUTGROUP_ACC" -format fasta > "$OUTGROUP_FASTA"

if [ ! -s "$OUTGROUP_FASTA" ]; then
    echo "ERRO: download do outgroup falhou ou veio vazio." >&2
    exit 1
fi
echo "   OK: $(head -1 "$OUTGROUP_FASTA")"

echo ">> Adicionando outgroup ao alinhamento existente (mafft --add)..."
mafft --add "$OUTGROUP_FASTA" --keeplength --thread -1 "$ALIGNED" > "$ALIGNED_ROOTED" 2> "$OUT_PREFIX.mafft_add.log"

n_seqs=$(grep -c "^>" "$ALIGNED_ROOTED")
echo "   Alinhamento agora tem $n_seqs sequencias (68 + 1 outgroup)."

OUTGROUP_HEADER=$(grep "^>$OUTGROUP_ACC" "$ALIGNED_ROOTED" | head -1 | sed 's/^>//' | awk '{print $1}')

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
echo ">> Pronto. Arvore enraizada em: ${OUT_PREFIX}.treefile"
echo "   Compare a topologia com a de 04_TREES/COI.treefile (sem outgroup) --"
echo "   os agrupamentos internos devem ser praticamente os mesmos; o que muda"
echo "   e a raiz/direcao da arvore ficar biologicamente correta agora."
