#!/bin/bash
# 01_alinhar_COI.sh
# Roda a partir de ~/Mansonella_projeto
#
# Uso: bash 01_alinhar_COI.sh
#
# O que faz:
#   1. Junta todos os 01_CURATED/<especie>/COI.fasta (conjunto principal,
#      SEM quarentena) num unico arquivo multi-especie.
#   2. Confere se ha accessions duplicados no arquivo final (nao deveria,
#      ja que a curadoria ja deduplicou, mas mais vale garantir).
#   3. Roda o MAFFT (--auto: MAFFT escolhe a estrategia de alinhamento mais
#      adequada pelo numero/tamanho das sequencias).
#   4. Reporta estatisticas basicas do alinhamento (numero de sequencias,
#      comprimento do alinhamento, % de gaps).

set -euo pipefail

RAW_BASE="01_CURATED"
OUT_DIR="03_ALIGNMENTS"
mkdir -p "$OUT_DIR"

COMBINED="$OUT_DIR/COI_all.fasta"
ALIGNED="$OUT_DIR/COI_aligned.fasta"

echo ">> Consolidando COI de todas as especies (exceto quarentena)..."
: > "$COMBINED"
find "$RAW_BASE" -name "COI.fasta" -not -path "*_QUARANTINE*" | sort | while read -r f; do
    especie=$(basename "$(dirname "$f")")
    n=$(grep -c "^>" "$f" || true)
    echo "   $especie: $n sequencias ($f)"
    cat "$f" >> "$COMBINED"
done

n_total=$(grep -c "^>" "$COMBINED" || true)
echo ">> Total consolidado: $n_total sequencias em $COMBINED"

echo ">> Checando duplicatas de accession no arquivo consolidado..."
n_unique=$(grep "^>" "$COMBINED" | awk '{print $1}' | sort -u | wc -l)
if [ "$n_total" != "$n_unique" ]; then
    echo "   ATENCAO: $n_total sequencias mas so $n_unique accessions unicos -- ha duplicata!"
    grep "^>" "$COMBINED" | awk '{print $1}' | sort | uniq -d
    exit 1
else
    echo "   OK: todos os $n_unique accessions sao unicos."
fi

if ! command -v mafft &> /dev/null; then
    echo ""
    echo "ERRO: mafft nao encontrado no PATH."
    echo "Instale com uma das opcoes:"
    echo "  mamba install -c bioconda mafft"
    echo "  sudo apt install mafft"
    exit 1
fi

echo ">> Rodando MAFFT (--auto)..."
mafft --auto --thread -1 "$COMBINED" > "$ALIGNED" 2> "$OUT_DIR/COI_mafft.log"
echo ">> Alinhamento salvo em: $ALIGNED"

echo ""
echo ">> Estatisticas do alinhamento:"
n_aligned=$(grep -c "^>" "$ALIGNED")
aln_len=$(awk '/^>/{if(seq)print length(seq); seq=""} !/^>/{seq=seq$0} END{print length(seq)}' "$ALIGNED" | head -1)
echo "   Sequencias alinhadas: $n_aligned"
echo "   Comprimento do alinhamento: $aln_len pb"
echo ""
echo "Proximo passo: inspecionar visualmente o alinhamento antes do IQ-TREE2."
echo "  - AliView (GUI, se tiver X11/WSLg): aliview $ALIGNED"
echo "  - ou so conferir zonas de gap excessivo com trimal -in $ALIGNED -out ${ALIGNED%.fasta}_summary.html -htmlout"
