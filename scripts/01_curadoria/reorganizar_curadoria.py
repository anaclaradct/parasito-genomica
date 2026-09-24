#!/usr/bin/env python3
"""
Reorganiza e deduplica os FASTAs de Mansonella a partir do CABECALHO (nunca da pasta),
porque a pasta demonstrou estar errada em pelo menos um caso (LT623914.1).

REGRAS APLICADAS:
1. Deduplicação por accession (com versao, ex. JF412322.1). Se o mesmo accession
   aparecer em mais de um arquivo/pasta, mantem 1 copia e loga o resto como duplicata.
2. Classificacao de especie e marcador SEMPRE a partir do texto do cabecalho, nunca
   do nome da pasta/arquivo.
3. Sequencias "Onchocercidae cf. Mansonella ozzardi" (identificacao incerta) vao para
   uma pasta separada 01_CURATED/Mansonella_ozzardi_cf_LAM2011/, isoladas do dataset
   "limpo" de M. ozzardi -- conforme principio de curadoria do projeto (nao assumir
   cf. como confirmado).
4. Pares que precisam SUA decisao manual (nao sao auto-resolvidos):
   a) RefSeq (prefixo NC_) vs accession original com mesmo isolate -> normalmente é
      a MESMA sequencia biologica sob dois accessions; script sugere manter so o RefSeq
      mas NAO decide sozinho.
   b) Mesmo isolate com genoma "complete" E "partial" -> possivel redundancia; listado
      para voce conferir (geralmente mantem so o complete).
   c) Organismo no cabecalho não bate com o esperado pela pasta/isolate (ex.: isolate
      "mdeux*" rotulado como Mansonella perstans) -> listado para checagem via efetch
      (taxid) antes de decidir a especie real.
   d) Mesmo voucher em accessions diferentes do MESMO marcador (nao mesmo caso de
      COI+12S do mesmo especime, que é esperado e correto).

SAIDA:
  01_CURATED/<Especie>/<marcador>.fasta      -> pronto para MAFFT
  01_CURATED/Mansonella_ozzardi_cf_LAM2011/  -> cf., separado
  02_LOGS/duplicatas_removidas.tsv
  02_LOGS/revisao_manual.tsv                 -> DECIDA antes de rodar MAFFT
"""

import os
import re
import sys
from collections import defaultdict

RAW_DIR = sys.argv[1] if len(sys.argv) > 1 else "00_RAW"
OUT_DIR = sys.argv[2] if len(sys.argv) > 2 else "01_CURATED"
LOG_DIR = sys.argv[3] if len(sys.argv) > 3 else "02_LOGS"

# ---------------------------------------------------------------------------
# DECISOES MANUAIS JA TOMADAS (bloco M. perstans / Mansonella sp. "DEUX"):
# Resolvido via TaxId (efetch docsum), nao apenas pelo texto do isolado --
# o nome do isolado ("mdeux2/3") e um codigo de amostra reaproveitado entre
# especies simpatricas no estudo original, NAO indica a especie sozinho.
#
# Confirmado por TaxId + tamanho (Slen):
#   OQ633017.1  TaxId 42231 (M. perstans)      13616 pb  completo
#   OQ633019.1  TaxId 42231 (M. perstans)      13617 pb  completo  <- 2o genoma real, nao duplicata
#   OQ633020.1  TaxId 1719275 (DEUX)           13619 pb  completo
#   OQ633021.1  TaxId 1719275 (DEUX)           13621 pb  completo
#   NC_077638.1 TaxId 42231, identico a OQ633017 (13616 pb)   -> duplicata exata
#   NC_077639.1 TaxId 1719275, identico a OQ633020 (13619 pb) -> duplicata exata
#   OQ633018.1  TaxId 42231 (M. perstans), mas PARCIAL (8978 pb) -> fora do
#               conjunto de genomas completos (procedencia ok, so nao e completo)
#   OQ633022.1  TaxId 2756192 -- DIFERENTE dos dois acima, nao identificado
#               nesta revisao -> excluido por principio (sem procedencia clara
#               confirmada). Checar antes de reincluir:
#               efetch -db taxonomy -id 2756192 -format docsum | \
#                 xtract -pattern DocumentSummary -element TaxId,ScientificName,Rank
#
#   INCLUIR: OQ633017.1, OQ633019.1, OQ633020.1, OQ633021.1
#   EXCLUIR: NC_077638.1, NC_077639.1 (duplicata exata)
#            OQ633018.1 (parcial -- procedencia ok, mas nao e genoma completo)
#            OQ633022.1 (TaxId nao resolvido nesta revisao)
# ---------------------------------------------------------------------------
EXCLUDE_ACCESSIONS = {
    "NC_077638.1", "NC_077639.1",   # duplicata EXATA confirmada por TaxId + tamanho (Slen) -- apagar
}

# Nao sao duplicatas, mas tambem nao entram no conjunto principal "limpo" --
# vao para uma pasta de quarentena em vez de serem apagados.
QUARANTINE_ACCESSIONS = {
    "OQ633018.1": "M. perstans (TaxId 42231) confirmado, mas genoma PARCIAL (8978 pb) -- "
                  "nao e duplicata de nada, so nao e um genoma completo",
    "OQ633022.1": "TaxId 2756192 = 'Mansonella sp.' (rank species), distinto de M. perstans "
                  "(42231) e de DEUX (1719275), mas sem nome formal/Lineage no NCBI -- taxon "
                  "legitimo porem nao identificavel com confianca; isolado para nao contaminar "
                  "nem M. perstans nem DEUX na arvore principal",
}

# Accessions cuja heuristica de "isolate mdeux* mas rotulado M. perstans"
# ja foi investigada e CONFIRMADA correta via TaxId (efetch docsum) --
# suprime o alerta repetido no log de revisao manual.
CONFIRMED_NOT_CONFLICTING = {
    "OQ633019.1",  # TaxId 42231 = M. perstans, confirmado; nome do isolado
                    # e so um codigo de amostra reaproveitado entre especies
                    # simpatricas no estudo original, nao indica organismo.
}

MARKER_RULES = [
    (r"cytochrome c? ?oxidase subunit ?1?I\b|COX1|\bCOI\b|coxI", "COI"),
    (r"small subunit ribosomal RNA|12S ribosomal RNA|containing 12S region|\b12S rRNA\b|partial 12S", "12S"),
    (r"5S ribosomal RNA", "5S"),
    (r"mitochondrion,? (complete|partial) genome", "MITOGENOME"),
]

def classify_marker(header):
    for pattern, marker in MARKER_RULES:
        if re.search(pattern, header, re.IGNORECASE):
            if marker == "MITOGENOME":
                completeness = "complete" if re.search(r"complete genome", header, re.IGNORECASE) else "partial"
                return f"MITOGENOME_{completeness}"
            return marker
    return "UNKNOWN_MARKER"

def classify_species(header):
    m = re.search(r"Onchocercidae cf\. (Mansonella \w+) ([A-Z0-9\-]+)", header)
    if m:
        return f"{m.group(1).replace(' ', '_')}_cf_{m.group(2)}"
    m = re.search(r"(Mansonella (?:sp\.? ?['\"]?[A-Z]+['\"]?|[a-z]+))", header)
    if m:
        sp = m.group(1).strip().replace("'", "").replace('"', "")
        return sp.replace(" ", "_")
    return "UNKNOWN_SPECIES"

def extract_isolate(header):
    m = re.search(r"isolate ([\w\-\./]+)", header)
    if m:
        return m.group(1)
    m = re.search(r"voucher ([\w\-\. ]+?)(?=\s+(?:cytochrome|12S|5S|small|mitochondrial|mitochondrion))", header)
    if m:
        return m.group(1).strip()
    return None

def parse_fasta(path):
    with open(path, "r", errors="replace") as f:
        header, seq = None, []
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header = line[1:]
                seq = []
            else:
                seq.append(line.strip())
        if header is not None:
            yield header, "".join(seq)

def main():
    records = []  # (accession, header, seq, source_file)
    for root, _, files in os.walk(RAW_DIR):
        for fn in files:
            if fn.lower().endswith((".fasta", ".fa", ".fna")):
                path = os.path.join(root, fn)
                for header, seq in parse_fasta(path):
                    acc = header.split()[0]
                    records.append((acc, header, seq, path))

    seen = {}
    dup_log = []
    excluded_log = []
    quarantined = {}
    for acc, header, seq, src in records:
        if acc in EXCLUDE_ACCESSIONS:
            excluded_log.append((acc, src, "exclusao definitiva documentada (ver EXCLUDE_ACCESSIONS no topo do script)"))
            continue
        if acc in seen or acc in quarantined:
            prior_src = seen.get(acc, quarantined.get(acc))[2]
            dup_log.append((acc, src, prior_src, "accession repetido -- mantida 1a ocorrencia"))
            continue
        if acc in QUARANTINE_ACCESSIONS:
            quarantined[acc] = (header, seq, src)
            continue
        seen[acc] = (header, seq, src)

    by_species_marker = defaultdict(list)  # (species, marker) -> [(acc, header, seq)]
    isolate_index = defaultdict(list)      # (species, isolate) -> [(acc, marker, header)]
    review = []

    for acc, (header, seq, src) in seen.items():
        species = classify_species(header)
        marker = classify_marker(header)
        isolate = extract_isolate(header)
        by_species_marker[(species, marker)].append((acc, header, seq))
        if isolate:
            isolate_index[(species, isolate)].append((acc, marker, header))

    # heuristica: acessos RefSeq (NC_) x nao-RefSeq com mesmo (species, isolate)
    for (species, isolate), entries in isolate_index.items():
        accs = [e[0] for e in entries]
        refseqs = [e for e in entries if e[0].startswith("NC_")]
        non_refseqs = [e for e in entries if not e[0].startswith("NC_")]
        markers_seen = defaultdict(list)
        for e in entries:
            markers_seen[e[1]].append(e[0])
        for marker, acclist in markers_seen.items():
            if len(acclist) > 1:
                if refseqs and non_refseqs and marker.startswith("MITOGENOME"):
                    review.append((species, isolate, marker, ";".join(acclist),
                                    "RefSeq (NC_) e accession original com mesmo isolate/marcador -- "
                                    "provavel MESMA sequencia; normalmente manter so o RefSeq"))
                elif marker in ("MITOGENOME_complete", "MITOGENOME_partial") or \
                     ("MITOGENOME_complete" in markers_seen and "MITOGENOME_partial" in markers_seen):
                    review.append((species, isolate, marker, ";".join(acclist),
                                    "mesmo isolate com >1 accession no mesmo marcador -- conferir se e "
                                    "complete vs partial redundante ou espécimes distintos"))
                else:
                    review.append((species, isolate, marker, ";".join(acclist),
                                    "mesmo isolate, mesmo marcador, accessions diferentes -- conferir manualmente"))

    # heuristica: organismo do cabecalho nao bate com isolate padrao de outra especie
    # (ex.: isolate "mdeux*" classificado como Mansonella_perstans)
    for (species, marker), entries in by_species_marker.items():
        for acc, header, seq in entries:
            isolate = extract_isolate(header) or ""
            if species == "Mansonella_perstans" and isolate.lower().startswith("mdeux") \
                    and acc not in CONFIRMED_NOT_CONFLICTING:
                review.append((species, isolate, marker, acc,
                                "isolate 'mdeux*' mas organismo no cabecalho = Mansonella perstans -- "
                                "confirmar taxid real via efetch antes de manter aqui "
                                "(suspeita: deveria ser Mansonella sp. 'DEUX')"))
            if species == "Mansonella_sp" and isolate.lower().startswith("mdeux"):
                review.append((species, isolate, marker, acc,
                                "cabecalho diz 'Mansonella sp.' (sem tag DEUX) mas isolate 'mdeux*' -- "
                                "provavelmente o MESMO isolado que aparece em Mansonella_sp._DEUX como "
                                "genoma completo; checar se esta sequencia parcial e redundante"))

    # heuristica extra: mesmo isolate aparecendo em BUCKETS DE ESPECIE diferentes
    # (cobre o caso Mansonella_sp x Mansonella_sp._DEUX, que o check acima por
    # (species, isolate) sozinho nao pega, pois a chave de especie difere)
    isolate_species_map = defaultdict(set)
    for (species, marker), entries in by_species_marker.items():
        for acc, header, seq in entries:
            isolate = extract_isolate(header)
            if isolate:
                isolate_species_map[isolate].add(species)
    for isolate, species_set in isolate_species_map.items():
        if len(species_set) > 1:
            if isolate.lower() == "mdeux3" and species_set == {"Mansonella_perstans", "Mansonella_sp._DEUX"}:
                continue  # confirmado via TaxId: sao 2 especimes reais, nao conflito de rotulagem
            review.append(("|".join(sorted(species_set)), isolate, "-", "-",
                            "mesmo isolate classificado em mais de uma especie (nomenclatura "
                            "inconsistente entre registros) -- confirmar taxonomia real via efetch"))

    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)

    for (species, marker), entries in by_species_marker.items():
        sp_dir = os.path.join(OUT_DIR, species)
        os.makedirs(sp_dir, exist_ok=True)
        out_path = os.path.join(sp_dir, f"{marker}.fasta")
        with open(out_path, "w") as out:
            for acc, header, seq in sorted(entries):
                out.write(f">{header}\n{seq}\n")

    with open(os.path.join(LOG_DIR, "duplicatas_removidas.tsv"), "w") as f:
        f.write("accession\tarquivo_duplicado\tarquivo_mantido\tmotivo\n")
        for row in dup_log:
            f.write("\t".join(row) + "\n")

    with open(os.path.join(LOG_DIR, "excluidos_manualmente.tsv"), "w") as f:
        f.write("accession\tarquivo\tmotivo\n")
        for row in excluded_log:
            f.write("\t".join(row) + "\n")

    quarantine_dir = os.path.join(OUT_DIR, "_QUARANTINE")
    os.makedirs(quarantine_dir, exist_ok=True)
    with open(os.path.join(LOG_DIR, "quarentena.tsv"), "w") as f:
        f.write("accession\tespecie\tmarcador\tmotivo\n")
        for acc, (header, seq, src) in quarantined.items():
            species = classify_species(header)
            marker = classify_marker(header)
            f.write(f"{acc}\t{species}\t{marker}\t{QUARANTINE_ACCESSIONS[acc]}\n")
            out_path = os.path.join(quarantine_dir, f"{acc}.fasta")
            with open(out_path, "w") as out:
                out.write(f">{header}\n{seq}\n")

    with open(os.path.join(LOG_DIR, "revisao_manual.tsv"), "w") as f:
        f.write("especie\tisolate\tmarcador\taccessions\tmotivo\n")
        for row in review:
            f.write("\t".join(row) + "\n")

    total_in = len(records)
    total_unique = len(seen)
    print(f"Registros lidos:        {total_in}")
    print(f"Sequencias unicas (set principal): {total_unique}")
    print(f"Em quarentena:          {len(quarantined)}  -> {OUT_DIR}/_QUARANTINE/ (log: {LOG_DIR}/quarentena.tsv)")
    print(f"Excluidas definitivamente: {len(excluded_log)}  -> {LOG_DIR}/excluidos_manualmente.tsv")
    print(f"Duplicatas removidas:   {len(dup_log)}  -> {LOG_DIR}/duplicatas_removidas.tsv")
    print(f"Casos p/ revisao:       {len(review)}   -> {LOG_DIR}/revisao_manual.tsv")
    print(f"FASTAs curados em:      {OUT_DIR}/<especie>/<marcador>.fasta")

if __name__ == "__main__":
    main()
