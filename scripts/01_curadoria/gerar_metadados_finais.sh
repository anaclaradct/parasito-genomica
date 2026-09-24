#!/bin/bash
# gerar_metadados_finais.sh
# Roda a partir de ~/Mansonella_projeto, DEPOIS de reorganizar_curadoria.py
#
# Uso: bash gerar_metadados_finais.sh > metadados_finais.tsv
#
# Junta: accession + especie/marcador (pela pasta em 01_CURATED, ja curada)
#        + host/country/isolate/collection_date/tamanho (via efetch, do GenBank)

set -euo pipefail

echo -e "accession\tespecie_pasta\tmarcador_pasta\tstatus\torganism_genbank\thost\tcountry\tisolate\tcollection_date\tlength_bp"

process_dir() {
    local base="$1"
    local status="$2"
    find "$base" -name "*.fasta" -not -path "*_QUARANTINE*" | while read -r f; do
        # pasta relativa dentro de 01_CURATED define especie/marcador
        rel="${f#01_CURATED/}"
        if [ "$status" = "curado" ]; then
            especie=$(dirname "$rel")
            marcador=$(basename "$rel" .fasta)
        else
            especie="_QUARANTINE"
            marcador="-"
        fi

        grep "^>" "$f" | sed 's/^>//' | awk '{print $1}' | while read -r acc; do
            GB=$(efetch -db nuccore -id "$acc" -format gb 2>/dev/null || echo "")
            if [ -z "$GB" ]; then
                echo -e "${acc}\t${especie}\t${marcador}\t${status}\tERRO_EFETCH\tNA\tNA\tNA\tNA\tNA"
                continue
            fi
            ORG=$(echo "$GB" | grep -m1 "ORGANISM" | sed 's/.*ORGANISM *//' | xargs || true)
            LEN=$(echo "$GB" | grep -m1 "^LOCUS" | awk '{print $3}' || true)
            HOST=$(echo "$GB" | grep -oP '/host="\K[^"]+' | head -1 || true)
            COUNTRY=$(echo "$GB" | grep -oP '/country="\K[^"]+' | head -1 || true)
            [ -z "$COUNTRY" ] && COUNTRY=$(echo "$GB" | grep -oP '/geo_loc_name="\K[^"]+' | head -1 || true)
            ISOLATE=$(echo "$GB" | grep -oP '/isolate="\K[^"]+' | head -1 || true)
            COLDATE=$(echo "$GB" | grep -oP '/collection_date="\K[^"]+' | head -1 || true)
            echo -e "${acc}\t${especie}\t${marcador}\t${status}\t${ORG:-NA}\t${HOST:-NA}\t${COUNTRY:-NA}\t${ISOLATE:-NA}\t${COLDATE:-NA}\t${LEN:-NA}"
            sleep 0.4
        done
    done
}

process_dir "01_CURATED" "curado"
find 01_CURATED/_QUARANTINE -name "*.fasta" 2>/dev/null | while read -r f; do
    acc=$(grep "^>" "$f" | sed 's/^>//' | awk '{print $1}')
    GB=$(efetch -db nuccore -id "$acc" -format gb 2>/dev/null || echo "")
    ORG=$(echo "$GB" | grep -m1 "ORGANISM" | sed 's/.*ORGANISM *//' | xargs || true)
    LEN=$(echo "$GB" | grep -m1 "^LOCUS" | awk '{print $3}' || true)
    HOST=$(echo "$GB" | grep -oP '/host="\K[^"]+' | head -1 || true)
    COUNTRY=$(echo "$GB" | grep -oP '/country="\K[^"]+' | head -1 || true)
            [ -z "$COUNTRY" ] && COUNTRY=$(echo "$GB" | grep -oP '/geo_loc_name="\K[^"]+' | head -1 || true)
    ISOLATE=$(echo "$GB" | grep -oP '/isolate="\K[^"]+' | head -1 || true)
    COLDATE=$(echo "$GB" | grep -oP '/collection_date="\K[^"]+' | head -1 || true)
    echo -e "${acc}\t_QUARANTINE\t-\tquarentena\t${ORG:-NA}\t${HOST:-NA}\t${COUNTRY:-NA}\t${ISOLATE:-NA}\t${COLDATE:-NA}\t${LEN:-NA}"
    sleep 0.4
done
