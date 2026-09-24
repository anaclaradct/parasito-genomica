#!/usr/bin/env python3
"""
diversidade_coi.py
Calcula estatisticas de diversidade a partir de um alinhamento FASTA --
so Python padrao, sem dependencias externas.

Uso:
    python3 diversidade_coi.py 03_ALIGNMENTS/COI_aligned.fasta [metadados.tsv]

Se um TSV de metadados for passado (mesmo formato do metadados_finais.tsv,
com colunas "accession" e "especie_pasta"), o script tambem quebra a
analise por especie, alem do resultado global.

METODOLOGIA (para reportar nos metodos):
  - Sitios conservados/variaveis/parsimony-informative: contados coluna a
    coluna, ignorando gaps na contagem de estados de cada coluna (um gap
    numa sequencia nao invalida a coluna, so nao entra na contagem de
    alelos daquela sequencia).
  - Haplotipos, diversidade haplotipica (Hd) e diversidade nucleotidica
    (pi): calculados sobre a regiao "core" -- colunas SEM NENHUM gap em
    NENHUMA sequencia do grupo analisado (seja o conjunto global, seja um
    subgrupo de especie). Isso e o equivalente ao "complete deletion" do
    DnaSP/MEGA. O tamanho dessa regiao core e reportado, porque datasets
    com sequencias de comprimentos bem diferentes (seu caso) reduzem essa
    regiao -- e isso deve ser explicito no texto de metodos.
  - Distancia p (p-distance) e Kimura 2-parametros (K2P) sao calculadas
    par a par, usando SOMENTE POSICOES SEM GAP NESSE PAR ESPECIFICO
    (seleção por pares/"pairwise deletion") -- diferente do calculo de
    haplotipos, que usa deleção completa. Isso preserva mais informação
    para a distância genética geral, ao custo de pares diferentes
    comparando números de sítios ligeiramente diferentes (padrão em
    MEGA/DnaSP).
"""
import sys
import csv
import math
from itertools import combinations
from collections import Counter

PURINES = set("AG")
PYRIMIDINES = set("CT")


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


def basic_site_stats(seqs):
    """seqs: lista de strings (mesmo comprimento, alinhadas)."""
    if not seqs:
        return 0, 0, 0, 0
    length = len(seqs[0])
    conserved = 0
    variable = 0
    parsimony_informative = 0
    for i in range(length):
        col = [s[i] for s in seqs if s[i] != "-"]
        if not col:
            continue
        counts = Counter(col)
        n_states = len(counts)
        if n_states <= 1:
            conserved += 1
        else:
            variable += 1
            n_allele_ge2 = sum(1 for c in counts.values() if c >= 2)
            if n_allele_ge2 >= 2:
                parsimony_informative += 1
    return length, conserved, variable, parsimony_informative


def complete_deletion_core(seqs):
    """Retorna as sequencias restritas as colunas sem NENHUM gap em NENHUMA sequencia,
    e TAMBEM os indices (posicoes no alinhamento original) mantidos -- para
    diagnostico de contiguidade e para permitir reuso do mesmo conjunto de
    colunas em outro grupo (core "comum")."""
    if not seqs:
        return [], 0, []
    length = len(seqs[0])
    keep_cols = [i for i in range(length) if all(s[i] != "-" for s in seqs)]
    core_seqs = ["".join(s[i] for i in keep_cols) for s in seqs]
    return core_seqs, len(keep_cols), keep_cols


def describe_contiguity(keep_cols):
    """Resume os blocos contiguos de posicoes mantidas -- para diagnostico."""
    if not keep_cols:
        return "nenhuma posicao"
    blocks = []
    start = prev = keep_cols[0]
    for c in keep_cols[1:]:
        if c == prev + 1:
            prev = c
            continue
        blocks.append((start, prev))
        start = prev = c
    blocks.append((start, prev))
    if len(blocks) == 1:
        return f"1 bloco continuo: posicoes {blocks[0][0]+1}-{blocks[0][1]+1} (1-indexed)"
    total_gap_internal = sum(blocks[i+1][0] - blocks[i][1] - 1 for i in range(len(blocks)-1))
    return (f"{len(blocks)} blocos (NAO continuo) -- ex: {blocks[0][0]+1}-{blocks[0][1]+1}, "
            f"{blocks[1][0]+1}-{blocks[1][1]+1}{'...' if len(blocks) > 2 else ''} "
            f"(total {total_gap_internal} pb de 'buracos' entre blocos)")


def core_seqs_from_cols(seqs, cols):
    """Extrai as mesmas colunas (indices) de um outro grupo, para comparacao
    direta de pi sobre a MESMA janela do gene."""
    return ["".join(s[i] for i in cols if i < len(s)) for s in seqs]


def haplotypes_and_diversity(core_seqs):
    n = len(core_seqs)
    if n < 2:
        return 0, 0.0, 0.0
    counts = Counter(core_seqs)
    n_hap = len(counts)
    freqs = [c / n for c in counts.values()]
    hd = (n / (n - 1)) * (1 - sum(f * f for f in freqs)) if n > 1 else 0.0

    total_pdist = 0.0
    n_pairs = 0
    L = len(core_seqs[0]) if core_seqs else 0
    for a, b in combinations(core_seqs, 2):
        diffs = sum(1 for x, y in zip(a, b) if x != y)
        total_pdist += diffs / L if L else 0
        n_pairs += 1
    pi = total_pdist / n_pairs if n_pairs else 0.0
    return n_hap, hd, pi


def pairwise_p_and_k2p(seqs):
    """p-distance e K2P medios, com delecao par-a-par (so posicoes sem gap NESSE par)."""
    p_total, k2p_total, n_pairs = 0.0, 0.0, 0
    for a, b in combinations(seqs, 2):
        sites = [(x, y) for x, y in zip(a, b) if x != "-" and y != "-"]
        L = len(sites)
        if L == 0:
            continue
        diffs = sum(1 for x, y in sites if x != y)
        p = diffs / L
        p_total += p

        transitions = sum(1 for x, y in sites if x != y and
                           ((x in PURINES and y in PURINES) or (x in PYRIMIDINES and y in PYRIMIDINES)))
        transversions = sum(1 for x, y in sites if x != y and
                             not ((x in PURINES and y in PURINES) or (x in PYRIMIDINES and y in PYRIMIDINES)))
        P = transitions / L
        Q = transversions / L
        try:
            k2p = -0.5 * math.log((1 - 2 * P - Q) * math.sqrt(1 - 2 * Q))
        except ValueError:
            k2p = float("nan")  # saturacao -- distancia grande demais pra formula K2P
        k2p_total += k2p
        n_pairs += 1
    if n_pairs == 0:
        return 0.0, 0.0, 0
    return p_total / n_pairs, k2p_total / n_pairs, n_pairs


def report_group(name, seqs, global_core_cols=None):
    print(f"\n=== {name} (n={len(seqs)}) ===")
    if len(seqs) < 2:
        print("  Menos de 2 sequencias -- pulando (nao da pra calcular diversidade).")
        return None
    length, conserved, variable, pinf = basic_site_stats(seqs)
    print(f"  Comprimento do alinhamento: {length} pb")
    print(f"  Sitios conservados:         {conserved} ({100*conserved/length:.1f}%)")
    print(f"  Sitios variaveis:           {variable} ({100*variable/length:.1f}%)")
    print(f"  Sitios parsimony-informative:{pinf} ({100*pinf/length:.1f}%)")

    core_seqs, core_len, keep_cols = complete_deletion_core(seqs)
    n_hap, hd, pi_core = haplotypes_and_diversity(core_seqs)
    print(f"  Regiao core p/ haplotipos (deletion completa, PROPRIA do grupo): {core_len} pb")
    print(f"  Contiguidade do core proprio: {describe_contiguity(keep_cols)}")
    print(f"  Numero de haplotipos (core proprio): {n_hap}")
    print(f"  Diversidade haplotipica (Hd, core proprio): {hd:.4f}")
    print(f"  Diversidade nucleotidica (pi, core PROPRIO -- NAO comparar direto entre grupos): {pi_core:.5f}")

    p_mean, k2p_mean, n_pairs = pairwise_p_and_k2p(seqs)
    print(f"  Distancia p media (delecao par-a-par, {n_pairs} pares): {p_mean:.5f}")
    print(f"  Distancia K2P media:        {k2p_mean:.5f}")

    pi_common = None
    if global_core_cols is not None:
        common_seqs = core_seqs_from_cols(seqs, global_core_cols)
        # se alguma sequencia tiver gap dentro das colunas do core global
        # (pode acontecer se o core global veio de OUTRO subconjunto), avisa.
        if any("-" in s for s in common_seqs):
            print(f"  [core COMUM ({len(global_core_cols)} pb)]: ATENCAO -- ainda ha gap dentro dessas "
                  f"colunas para este grupo; pi sobre core comum nao calculado (nao e deletion completa valida aqui).")
        else:
            _, _, pi_common = haplotypes_and_diversity(common_seqs)
            print(f"  Diversidade nucleotidica (pi, core COMUM de {len(global_core_cols)} pb -- "
                  f"COMPARAVEL entre grupos): {pi_common:.5f}")

    return {"core_len": core_len, "keep_cols": keep_cols, "pi_core": pi_core, "pi_common": pi_common,
            "n_hap": n_hap, "hd": hd}


def main():
    if len(sys.argv) < 2:
        print("Uso: python3 diversidade_coi.py <alinhamento.fasta> [metadados.tsv]")
        sys.exit(1)

    fasta_path = sys.argv[1]
    meta_path = sys.argv[2] if len(sys.argv) > 2 else None

    records = parse_fasta(fasta_path)
    accs = [h for h, s in records]
    seqs = [s for h, s in records]

    global_result = report_group("GLOBAL (todas as sequencias)", seqs)
    global_core_cols = global_result["keep_cols"] if global_result else None

    if meta_path:
        acc2species = {}
        with open(meta_path, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                acc2species[row["accession"]] = row["especie_pasta"]

        by_species = {}
        for acc, seq in zip(accs, seqs):
            sp = acc2species.get(acc, "DESCONHECIDA")
            by_species.setdefault(sp, []).append(seq)

        for sp, sp_seqs in sorted(by_species.items(), key=lambda x: -len(x[1])):
            report_group(sp, sp_seqs, global_core_cols=global_core_cols)


if __name__ == "__main__":
    main()
