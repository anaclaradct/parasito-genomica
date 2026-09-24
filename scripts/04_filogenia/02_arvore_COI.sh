#!/bin/bash
set -euo pipefail

ALIGNED="03_ALIGNMENTS/COI_aligned.fasta"
OUT_PREFIX="04_TREES/COI"

mkdir -p 04_TREES

if [ ! -f "$ALIGNED" ]; then
    echo "ERRO: nao encontrei $ALIGNED. Rode scripts/02_alinhamento/01_alinhar_COI.sh primeiro." >&2
    exit 1
fi

IQTREE_BIN=""
for candidate in iqtree3 iqtree2 iqtree; do
    if command -v "$candidate" &> /dev/null; then
        IQTREE_BIN="$candidate"
        break
    fi
done

if [ -z "$IQTREE_BIN" ]; then
    echo "ERRO: nenhum binario iqtree/iqtree2/iqtree3 encontrado no PATH." >&2
    exit 1
fi

echo ">> Usando binario: $IQTREE_BIN ($(command -v $IQTREE_BIN))"
echo ">> Rodando em $ALIGNED..."
"$IQTREE_BIN" -s "$ALIGNED" -m MFP -alrt 1000 -bb 1000 -nt AUTO -pre "$OUT_PREFIX"

echo ""
echo ">> Pronto. Arquivos principais:"
echo "   ${OUT_PREFIX}.treefile   -- arvore final (newick), com suporte SH-aLRT/UFBoot nos nos"
echo "   ${OUT_PREFIX}.iqtree     -- relatorio completo (modelo escolhido, log-likelihood, etc.)"
echo "   ${OUT_PREFIX}.log        -- log da execucao"
echo ""
echo "Para visualizar a arvore: iTOL (itol.embl.de, upload do .treefile) ou FigTree."
