#!/usr/bin/env python3
"""
analises_complementares_coi.py

Quatro analises complementares sobre o alinhamento de COI ja curado:
  1. Barcoding gap (distancia max intraespecifica vs min interespecifica)
  2. Testes de neutralidade: Tajima's D e Fu's Fs + distribuicao de
     diferencas pareadas (mismatch distribution)
  3. Grafico de saturacao (transicoes/transversoes vs distancia)
  4. AMOVA / Phi_ST entre grupos geograficos (com teste de permutacao)

Uso:
    python3 analises_complementares_coi.py <alinhamento.fasta> <metadados.tsv>

Todas as funcoes usam SO Python padrao (sem numpy/scipy), pensado pra
rodar em qualquer maquina sem instalar nada.
"""
import sys
import os
import csv
import math
import random
from itertools import combinations
from collections import defaultdict, Counter

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from fasta_utils import parse_fasta

PURINES = set("AG")
PYRIMIDINES = set("CT")


def pairwise_diff_count(a, b):
    """Numero de diferencas (contagem bruta, NAO normalizada) usando
    delecao par-a-par (so posicoes sem gap em nenhum dos dois)."""
    sites = [(x, y) for x, y in zip(a, b) if x != '-' and y != '-']
    if not sites:
        return None, 0
    diffs = sum(1 for x, y in sites if x != y)
    return diffs, len(sites)


def pairwise_pdist(a, b):
    diffs, L = pairwise_diff_count(a, b)
    if diffs is None or L == 0:
        return None
    return diffs / L


def transitions_transversions(a, b):
    """Retorna (transicoes, transversoes, total_sitios_comparados) -- o total
    inclui posicoes IGUAIS tambem (nao so as diferentes), necessario para
    calcular corretamente p-distancia = (ts+tv)/total."""
    compared = [(x, y) for x, y in zip(a, b) if x != '-' and y != '-']
    diffs = [(x, y) for x, y in compared if x != y]
    ts = sum(1 for x, y in diffs if (x in PURINES and y in PURINES) or (x in PYRIMIDINES and y in PYRIMIDINES))
    tv = len(diffs) - ts
    return ts, tv, len(compared)


def load_by_group(fasta_path, meta_path, group_col, species_filter=None):
    records = parse_fasta(fasta_path)
    meta = {}
    with open(meta_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            meta[row["accession"]] = row

    by_group = defaultdict(list)
    for acc, seq in records:
        row = meta.get(acc)
        if row is None:
            continue
        if species_filter and row.get("especie_pasta") != species_filter:
            continue
        g_raw = row.get(group_col, "NA")
        g = g_raw.split(":")[0].strip() if g_raw else "NA"
        if not g or g == "NA":
            g = "Sem_dado"
        by_group[g].append((acc, seq))
    return by_group, records, meta


# ---------------------------------------------------------------------------
# 1. BARCODING GAP
# ---------------------------------------------------------------------------

def barcoding_gap(records, meta):
    by_species = defaultdict(list)
    for acc, seq in records:
        row = meta.get(acc)
        if row is None:
            continue
        by_species[row.get("especie_pasta", "NA")].append((acc, seq))

    species_list = [s for s in by_species if len(by_species[s]) >= 2]
    results = []
    for sp in species_list:
        seqs_sp = by_species[sp]
        # maior distancia intraespecifica
        max_intra = 0.0
        for (a1, s1), (a2, s2) in combinations(seqs_sp, 2):
            d = pairwise_pdist(s1, s2)
            if d is not None and d > max_intra:
                max_intra = d

        # menor distancia interespecifica (contra qualquer outra especie)
        min_inter = None
        min_inter_sp = None
        others = [(acc, seq, osp) for osp, seqs in by_species.items() if osp != sp for acc, seq in seqs]
        for a1, s1 in seqs_sp:
            for a2, s2, osp in others:
                d = pairwise_pdist(s1, s2)
                if d is not None and (min_inter is None or d < min_inter):
                    min_inter = d
                    min_inter_sp = osp

        gap = (min_inter - max_intra) if min_inter is not None else None
        results.append({
            "especie": sp, "n": len(seqs_sp),
            "max_dist_intra": max_intra,
            "min_dist_inter": min_inter,
            "especie_mais_proxima": min_inter_sp,
            "gap": gap,
            "ha_gap": (gap is not None and gap > 0),
        })
    return results


# ---------------------------------------------------------------------------
# 2a. TAJIMA'S D
# ---------------------------------------------------------------------------

def count_segregating_sites(seqs):
    """S = numero de posicoes variaveis, considerando so colunas SEM gap em
    nenhuma sequencia do grupo (deleção completa) -- necessario pro Tajima's D
    classico, que assume dados completos."""
    if not seqs:
        return 0, 0
    length = len(seqs[0])
    keep_cols = [i for i in range(length) if all(s[i] != '-' for s in seqs)]
    S = 0
    for i in keep_cols:
        states = set(s[i] for s in seqs)
        if len(states) > 1:
            S += 1
    return S, len(keep_cols)


def mean_pairwise_diff_raw(seqs, core_len):
    """k = media do NUMERO BRUTO de diferencas par a par (nao normalizado
    por comprimento) -- e o 'k' classico usado no Tajima's D, calculado
    sobre a mesma regiao core usada para S."""
    length = len(seqs[0])
    keep_cols = [i for i in range(length) if all(s[i] != '-' for s in seqs)]
    core_seqs = ["".join(s[i] for i in keep_cols) for s in seqs]
    total = 0
    n_pairs = 0
    for a, b in combinations(core_seqs, 2):
        diffs = sum(1 for x, y in zip(a, b) if x != y)
        total += diffs
        n_pairs += 1
    return (total / n_pairs) if n_pairs else 0.0


def tajimas_d(seqs):
    """Formula classica de Tajima (1989). Retorna (D, S, k, n) ou
    (None, S, k, n) se D nao puder ser calculado (S=0 ou variancia<=0)."""
    n = len(seqs)
    if n < 4:
        return None, None, None, n  # Tajima's D exige n razoavel

    S, core_len = count_segregating_sites(seqs)
    if S == 0:
        return None, S, 0.0, n

    k = mean_pairwise_diff_raw(seqs, core_len)

    a1 = sum(1 / i for i in range(1, n))
    a2 = sum(1 / (i * i) for i in range(1, n))
    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n * n + n + 3) / (9 * n * (n - 1))
    c1 = b1 - 1 / a1
    c2 = b2 - (n + 2) / (a1 * n) + a2 / (a1 * a1)
    e1 = c1 / a1
    e2 = c2 / (a1 * a1 + a2)

    var_d = e1 * S + e2 * S * (S - 1)
    if var_d <= 0:
        return None, S, k, n

    D = (k - S / a1) / math.sqrt(var_d)
    return D, S, k, n


# ---------------------------------------------------------------------------
# 2b. FU'S FS (via formula de amostragem de Ewens)
# ---------------------------------------------------------------------------

def unsigned_stirling_first_kind_row(n):
    """Retorna lista [s(n,0), s(n,1), ..., s(n,n)] dos numeros de Stirling
    de primeira especie NAO-SINALIZADOS, via recorrencia:
    s(n,k) = s(n-1,k-1) + (n-1)*s(n-1,k)
    Usa inteiros grandes (Python trata automaticamente)."""
    row = [1]  # s(0,0) = 1
    for m in range(1, n + 1):
        new_row = [0] * (m + 1)
        for k in range(0, m + 1):
            term1 = row[k - 1] if k - 1 >= 0 and k - 1 < len(row) else 0
            term2 = (m - 1) * row[k] if k < len(row) else 0
            new_row[k] = term1 + term2
        row = new_row
    return row


def prob_k_alleles_ewens(n, theta, stirling_row):
    """P(K=j | n, theta) para j=1..n, via formula de amostragem de Ewens:
    P(K=j) = |s(n,j)| * theta^j / theta_(n)
    onde theta_(n) = theta*(theta+1)*...*(theta+n-1) (fatorial ascendente).
    Trabalha em log-space para evitar overflow com n e theta maiores."""
    log_theta = math.log(theta)
    log_rising_fact = sum(math.log(theta + i) for i in range(n))  # log(theta_(n))
    probs = []
    for j in range(1, n + 1):
        s_nj = stirling_row[j]
        if s_nj == 0:
            probs.append(0.0)
            continue
        log_p = math.log(s_nj) + j * log_theta - log_rising_fact
        probs.append(math.exp(log_p))
    return probs  # indice 0 -> j=1, indice n-1 -> j=n


def fus_fs(seqs):
    """Fu's Fs (Fu 1997). Retorna (Fs, k_haplotipos_observados, theta_pi, n)."""
    n = len(seqs)
    if n < 3:
        return None, None, None, n

    S, core_len = count_segregating_sites(seqs)
    length = len(seqs[0])
    keep_cols = [i for i in range(length) if all(s[i] != '-' for s in seqs)]
    core_seqs = ["".join(s[i] for i in keep_cols) for s in seqs]

    n_haplotypes = len(set(core_seqs))
    k = mean_pairwise_diff_raw(seqs, core_len)  # theta_pi = k (numero medio de diferencas)
    theta_pi = k
    if theta_pi <= 0:
        return None, n_haplotypes, theta_pi, n

    stirling_row = unsigned_stirling_first_kind_row(n)
    probs = prob_k_alleles_ewens(n, theta_pi, stirling_row)  # P(K=1..n)

    # S' = P(K >= k_obs) = soma das probabilidades de j=k_obs ate n
    s_prime = sum(probs[j - 1] for j in range(n_haplotypes, n + 1))
    s_prime = min(max(s_prime, 1e-300), 1 - 1e-15)  # evita log(0) ou log(negativo)

    Fs = math.log(s_prime / (1 - s_prime))
    return Fs, n_haplotypes, theta_pi, n


# ---------------------------------------------------------------------------
# 2c. MISMATCH DISTRIBUTION
# ---------------------------------------------------------------------------

def mismatch_distribution(seqs):
    diffs_list = []
    for a, b in combinations(seqs, 2):
        d, L = pairwise_diff_count(a, b)
        if d is not None:
            diffs_list.append(d)
    if not diffs_list:
        return {}, 0.0
    hist = Counter(diffs_list)
    mean_diff = sum(diffs_list) / len(diffs_list)
    return dict(sorted(hist.items())), mean_diff


# ---------------------------------------------------------------------------
# 3. SATURACAO (transicoes/transversoes vs distancia)
# ---------------------------------------------------------------------------

def saturation_data(records_subset):
    """records_subset: lista de (acc, seq). Retorna lista de tuplas
    (p_distancia, k2p_distancia, transicoes, transversoes, n_sitios)
    para TODOS os pares -- para plotar Ts e Tv contra a distancia K2P."""
    rows = []
    for (a1, s1), (a2, s2) in combinations(records_subset, 2):
        ts, tv, L = transitions_transversions(s1, s2)
        if L == 0:
            continue
        p = (ts + tv) / L
        P = ts / L
        Q = tv / L
        try:
            k2p = -0.5 * math.log((1 - 2 * P - Q) * math.sqrt(1 - 2 * Q))
        except ValueError:
            k2p = float("nan")
        rows.append((p, k2p, ts / L, tv / L))
    return rows


# ---------------------------------------------------------------------------
# 4. AMOVA / Phi_ST (com teste de permutacao)
# ---------------------------------------------------------------------------

def amova_phi_st(by_group, n_permutations=1000, seed=42):
    """AMOVA simplificado (Excoffier et al. 1992) sobre distancias par a par
    (numero bruto de diferencas, delecao par-a-par). Retorna
    (phi_st, p_valor, detalhes)."""
    all_seqs = []  # (seq, grupo)
    for g, seqs in by_group.items():
        for acc, seq in seqs:
            all_seqs.append((seq, g))

    N = len(all_seqs)
    groups = list(by_group.keys())
    k = len(groups)
    if N < 3 or k < 2:
        return None, None, {"erro": "dados insuficientes (precisa >=2 grupos, >=3 sequencias)"}

    # SSD total independe da atribuicao de grupo (so depende do conjunto
    # fixo de sequencias), entao e calculado uma unica vez aqui fora,
    # em vez de dentro de compute_phi() a cada uma das n_permutations chamadas.
    ssd_total = 0.0
    for (s1, _), (s2, _) in combinations(all_seqs, 2):
        d, L = pairwise_diff_count(s1, s2)
        if d is not None:
            ssd_total += d * d
    ssd_total /= N

    def compute_phi(assignment):
        # assignment: lista paralela a all_seqs, com o grupo de cada um

        # SSD dentro dos grupos
        by_g = defaultdict(list)
        for (seq, _), g in zip(all_seqs, assignment):
            by_g[g].append(seq)

        ssd_within = 0.0
        group_sizes = {}
        for g, seqs in by_g.items():
            group_sizes[g] = len(seqs)
            if len(seqs) < 2:
                continue
            s = 0.0
            for s1, s2 in combinations(seqs, 2):
                d, L = pairwise_diff_count(s1, s2)
                if d is not None:
                    s += d * d
            ssd_within += s / len(seqs)

        ssd_among = ssd_total - ssd_within

        df_among = k - 1
        df_within = N - k
        if df_among <= 0 or df_within <= 0:
            return None

        msd_among = ssd_among / df_among
        msd_within = ssd_within / df_within

        sum_ng2_over_N = sum((ng * ng) / N for ng in group_sizes.values())
        n0 = (N - sum_ng2_over_N) / (k - 1)
        if n0 <= 0:
            return None

        var_within = msd_within
        var_among = (msd_among - msd_within) / n0

        denom = var_among + var_within
        if denom <= 0:
            return 0.0
        phi = var_among / denom
        return phi

    original_assignment = [g for _, g in all_seqs]
    phi_obs = compute_phi(original_assignment)
    if phi_obs is None:
        return None, None, {"erro": "nao foi possivel calcular Phi_ST (verifique tamanho dos grupos)"}

    random.seed(seed)
    count_ge = 0
    for _ in range(n_permutations):
        perm = original_assignment[:]
        random.shuffle(perm)
        phi_perm = compute_phi(perm)
        if phi_perm is not None and phi_perm >= phi_obs:
            count_ge += 1
    p_value = count_ge / n_permutations

    details = {"N": N, "k_grupos": k, "grupos": {g: len(s) for g, s in by_group.items()}}
    return phi_obs, p_value, details


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) < 3:
        print("Uso: python3 analises_complementares_coi.py <alinhamento.fasta> <metadados.tsv>")
        sys.exit(1)

    fasta_path, meta_path = sys.argv[1], sys.argv[2]
    records = parse_fasta(fasta_path)
    meta = {}
    with open(meta_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            meta[row["accession"]] = row

    by_species = defaultdict(list)
    for acc, seq in records:
        row = meta.get(acc)
        if row is None:
            continue
        by_species[row.get("especie_pasta", "NA")].append((acc, seq))

    # ---------------- 1. BARCODING GAP ----------------
    print("=" * 70)
    print("1. BARCODING GAP")
    print("=" * 70)
    gap_results = barcoding_gap(records, meta)
    for r in sorted(gap_results, key=lambda x: -x["n"]):
        gap_str = f"{r['gap']:.4f}" if r["gap"] is not None else "NA"
        status = "GAP OK" if r["ha_gap"] else "SEM GAP (sobreposicao!)"
        print(f"  {r['especie']:35s} n={r['n']:3d}  max_intra={r['max_dist_intra']:.4f}  "
              f"min_inter={r['min_dist_inter']:.4f} (vs {r['especie_mais_proxima']})  "
              f"gap={gap_str}  [{status}]")

    # ---------------- 2. NEUTRALIDADE ----------------
    print()
    print("=" * 70)
    print("2. TESTES DE NEUTRALIDADE E MISMATCH DISTRIBUTION (por especie, n>=4)")
    print("=" * 70)
    for sp, seqs_sp in sorted(by_species.items(), key=lambda x: -len(x[1])):
        seqs = [s for _, s in seqs_sp]
        n = len(seqs)
        if n < 4:
            print(f"\n  {sp} (n={n}): pulado (Tajima's D precisa n>=4)")
            continue
        D, S, k, _ = tajimas_d(seqs)
        Fs, n_hap, theta_pi, _ = fus_fs(seqs)
        hist, mean_diff = mismatch_distribution(seqs)

        print(f"\n  --- {sp} (n={n}) ---")
        print(f"  Sitios segregantes (S, regiao core): {S}")
        print(f"  Tajima's D: {D:.4f}" if D is not None else "  Tajima's D: nao calculavel")
        print(f"  Fu's Fs: {Fs:.4f} (haplotipos={n_hap}, theta_pi={theta_pi:.3f})" if Fs is not None else "  Fu's Fs: nao calculavel")
        print(f"  Mismatch distribution (nº diferencas : frequencia): {hist}")
        print(f"  Media de diferencas pareadas: {mean_diff:.2f}")
        print(f"  Interpretacao rapida: D e Fs muito negativos + distribuicao unimodal")
        print(f"    = assinatura de expansao populacional recente e/ou selecao purificadora.")

    # ---------------- 3. SATURACAO ----------------
    print()
    print("=" * 70)
    print("3. DADOS DE SATURACAO (todas as 68 sequencias, entre especies)")
    print("=" * 70)
    sat_rows = saturation_data(records)
    out_sat = "saturacao_coi.tsv"
    with open(out_sat, "w") as f:
        f.write("p_distancia\tk2p_distancia\tprop_transicoes\tprop_transversoes\n")
        for p_d, k2p, prop_ts, prop_tv in sat_rows:
            f.write(f"{p_d:.5f}\t{k2p:.5f}\t{prop_ts:.5f}\t{prop_tv:.5f}\n")
    print(f"  {len(sat_rows)} pares calculados. Dados salvos em: {out_sat}")
    print(f"  Para visualizar: plote 'prop_transicoes' e 'prop_transversoes' (eixo Y)")
    print(f"  contra 'k2p_distancia' (eixo X). Se a curva de transicoes ACHATAR/CAIR")
    print(f"  nas distancias maiores (em vez de continuar subindo), e sinal de saturacao.")

    # ---------------- 4. AMOVA ----------------
    print()
    print("=" * 70)
    print("4. AMOVA / Phi_ST por pais (M. perstans e M. ozzardi)")
    print("=" * 70)
    for sp in ["Mansonella_perstans", "Mansonella_ozzardi"]:
        if sp not in by_species:
            continue
        by_country = defaultdict(list)
        for acc, seq in by_species[sp]:
            row = meta.get(acc, {})
            c_raw = row.get("country", "NA")
            c = c_raw.split(":")[0].strip() if c_raw and c_raw != "NA" else "Sem_dado"
            by_country[c].append((acc, seq))
        # remove grupos com so 1 sequencia isolada demais? AMOVA tolera, mas avisa
        print(f"\n  --- {sp}: grupos por pais ---")
        phi, pval, details = amova_phi_st(by_country, n_permutations=1000)
        if phi is None:
            print(f"  Nao calculavel: {details.get('erro')}")
            continue
        print(f"  N={details['N']}, grupos={details['k_grupos']} {details['grupos']}")
        print(f"  Phi_ST = {phi:.4f}  (p={pval:.4f}, {1000} permutacoes)")
        if pval < 0.05:
            print(f"  -> Estatisticamente significativo: HA estruturacao genetica por pais.")
        else:
            print(f"  -> NAO significativo: sem evidencia de estruturacao por pais neste teste.")

    print()
    print("Concluido.")


if __name__ == "__main__":
    main()
