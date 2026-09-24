#!/usr/bin/env python3
"""
buscar_numts.py
Varredura sistematica de possiveis NUMTs/pseudogenes em TODAS as sequencias
de COI curadas, combinando 3 sinais independentes (nenhum sozinho é prova
definitiva -- por isso os tres juntos):

  1. STOP CODON INTERNO: traduz nos 3 frames (codigo genetico mitocondrial
     de invertebrados, tabela 5 do NCBI) e reporta o menor numero de stops
     internos entre os frames. Sequencia funcional real nao devia ter
     nenhum, em nenhum frame de leitura correto.

  2. COMPOSICAO DE BASES (%AT): mtDNA de nematoide e tipicamente rico em
     A+T. Sequencias com %AT muito abaixo da media do proprio dataset
     podem ter vindo (ou estar mais proximas) de DNA nuclear.

  3. DISTANCIA ANORMAL DENTRO DA PROPRIA ESPECIE: calcula a distancia p
     media de cada sequencia para as outras da MESMA especie (deleção
     par-a-par). Sequencias com distancia media muito maior que as demais
     da mesma especie (outlier) sao suspeitas -- e exatamente o padrao que
     apareceu na rede de M. ozzardi.

Uso:
    python3 buscar_numts.py <alinhamento.fasta> <metadados.tsv>

Saida: buscar_numts_resultado.tsv com uma linha por sequencia e uma coluna
de flags (pode ter 0, 1, 2 ou 3 flags acesas). Sequencias com 2+ flags
merecem atencao prioritaria.

IMPORTANTE: isso e uma triagem, nao um veredito automatico. Sequencia
flagada deve ser investigada manualmente (ex: comparar com o registro
GenBank completo, checar se ha nota do depositante, ou fazer BLAST contra
um mitogenoma de referencia) antes de decidir excluir.
"""
import sys
import csv
import statistics as st
from itertools import combinations
from collections import defaultdict

CODON_TABLE_5 = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'M','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'W','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'S','AGG':'S','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

def parse_fasta(path):
    records = []
    header, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq).upper()))
                header = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if header is not None:
            records.append((header, "".join(seq).upper()))
    return records

def min_internal_stops(seq_nogap):
    best = None
    for frame in (0, 1, 2):
        seq = seq_nogap[frame:]
        n_codons = len(seq) // 3
        stops = 0
        for i in range(n_codons):
            codon = seq[i*3:i*3+3]
            if 'N' in codon or len(codon) < 3:
                continue
            aa = CODON_TABLE_5.get(codon, 'X')
            if aa == '*' and i < n_codons - 2:  # ignora stop nos 2 ultimos codons (pode ser fim real)
                stops += 1
        if best is None or stops < best:
            best = stops
    return best

def at_content(seq_nogap):
    if not seq_nogap:
        return 0
    at = sum(1 for c in seq_nogap if c in "AT")
    return 100 * at / len(seq_nogap)

def pairwise_pdist(a, b):
    sites = [(x, y) for x, y in zip(a, b) if x != '-' and y != '-']
    if not sites:
        return None
    diffs = sum(1 for x, y in sites if x != y)
    return diffs / len(sites)

def main():
    if len(sys.argv) < 3:
        print("Uso: python3 buscar_numts.py <alinhamento.fasta> <metadados.tsv>")
        sys.exit(1)

    fasta_path, meta_path = sys.argv[1], sys.argv[2]
    records = parse_fasta(fasta_path)
    meta = {}
    with open(meta_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            meta[row["accession"]] = row

    # organiza por especie
    by_species = defaultdict(list)  # especie -> [(acc, seq)]
    for acc, seq in records:
        row = meta.get(acc)
        if row is None:
            continue
        sp = row.get("especie_pasta", "NA")
        by_species[sp].append((acc, seq))

    results = []

    # --- sinal 1 e 2: por sequencia individual ---
    at_values_all = []
    per_seq = {}
    for acc, seq in records:
        seq_nogap = seq.replace("-", "")
        stops = min_internal_stops(seq_nogap)
        at = at_content(seq_nogap)
        per_seq[acc] = {"stops": stops, "at": at}
        at_values_all.append(at)

    at_mean_all = st.mean(at_values_all)
    at_sd_all = st.stdev(at_values_all) if len(at_values_all) > 1 else 0

    # --- sinal 3: distancia media dentro da propria especie ---
    dist_by_acc = {}
    for sp, seqs in by_species.items():
        if len(seqs) < 3:
            continue  # precisa de pelo menos 3 pra "media do grupo" fazer sentido
        accs = [a for a, _ in seqs]
        seqmap = dict(seqs)
        avg_dist = {a: [] for a in accs}
        for a1, a2 in combinations(accs, 2):
            d = pairwise_pdist(seqmap[a1], seqmap[a2])
            if d is not None:
                avg_dist[a1].append(d)
                avg_dist[a2].append(d)
        means = {a: (st.mean(v) if v else None) for a, v in avg_dist.items()}
        valid_means = [m for m in means.values() if m is not None]
        if len(valid_means) < 3:
            continue
        group_mean = st.mean(valid_means)
        group_sd = st.stdev(valid_means)
        for a, m in means.items():
            if m is None:
                continue
            dist_by_acc[a] = {"mean_dist": m, "group_mean": group_mean,
                               "group_sd": group_sd, "species": sp}

    # --- monta resultado final ---
    for acc, seq in records:
        row = meta.get(acc, {})
        sp = row.get("especie_pasta", "NA")
        info = per_seq[acc]
        flags = []

        if info["stops"] and info["stops"] > 0:
            flags.append(f"STOP_INTERNO({info['stops']})")

        z_at = (info["at"] - at_mean_all) / at_sd_all if at_sd_all else 0
        if z_at < -2:
            flags.append(f"AT_BAIXO(z={z_at:.1f})")

        d = dist_by_acc.get(acc)
        if d and d["group_sd"] > 0:
            z_dist = (d["mean_dist"] - d["group_mean"]) / d["group_sd"]
            if z_dist > 2:
                flags.append(f"DISTANCIA_OUTLIER(z={z_dist:.1f})")

        results.append({
            "accession": acc, "especie": sp,
            "pct_AT": f"{info['at']:.1f}",
            "stops_internos_melhor_frame": info["stops"],
            "dist_media_propria_especie": f"{d['mean_dist']:.4f}" if d else "NA",
            "n_flags": len(flags),
            "flags": ";".join(flags) if flags else "",
        })

    results.sort(key=lambda r: -r["n_flags"])

    out_path = "buscar_numts_resultado.tsv"
    with open(out_path, "w") as f:
        cols = ["accession", "especie", "pct_AT", "stops_internos_melhor_frame",
                "dist_media_propria_especie", "n_flags", "flags"]
        f.write("\t".join(cols) + "\n")
        for r in results:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    n_flagged = sum(1 for r in results if r["n_flags"] > 0)
    n_priority = sum(1 for r in results if r["n_flags"] >= 2)
    print(f"Total de sequencias analisadas: {len(results)}")
    print(f"Com pelo menos 1 sinal de alerta: {n_flagged}")
    print(f"Com 2+ sinais (prioridade de revisao): {n_priority}")
    print(f"\nResultado completo em: {out_path}")
    if n_priority:
        print("\n--- Sequencias com 2+ sinais ---")
        for r in results:
            if r["n_flags"] >= 2:
                print(f"  {r['accession']} ({r['especie']}): {r['flags']}")

if __name__ == "__main__":
    main()
