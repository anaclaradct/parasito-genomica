#!/usr/bin/env python3
"""
nucleo_comum_alinhamento.py - recorta um alinhamento FASTA para o "nucleo comum" por
criterio OBJETIVO de cobertura (sem olhar resultados de saturacao/filogenia).

Criterio (2 passos):
  1) mantem as colunas cobertas (base, nao gap/?/N) por >= --min-col das sequencias;
  2) remove as sequencias que cobrem < --min-seq das colunas mantidas;
  depois apaga colunas que ficaram so com gaps.
  --maior-bloco: fica so com a maior corrida CONTIGUA de colunas do passo 1
  (util se quiser posicoes de codon / regiao continua).

Boas praticas contra vies:
  * rode primeiro com --varredura (nao escreve nada) e escolha o corte olhando a
    tabela de colunas/sequencias que sobram, NAO o resultado do DAMBE;
  * registre o corte escolhido e reporte tambem o alinhamento completo;
  * confira no relatorio (--relatorio) quais sequencias sairam (nenhuma especie deve
    desaparecer sem voce saber);
  * lembre: menos sitios = menos poder para detectar saturacao.

Exemplos:
  python3 nucleo_comum_alinhamento.py --fasta COI_aligned_sem_NUMT.fasta --varredura
  python3 nucleo_comum_alinhamento.py --fasta COI_aligned_sem_NUMT.fasta \
      --min-col 0.9 --min-seq 0.9 --saida COI_nucleo.fasta --relatorio COI_nucleo_seqs.tsv
"""
import argparse
import sys
from filo_utils import ler_fasta, escrever_fasta

GAPS = set("-?N")
EPS = 1e-9


def cobertura_colunas(seqs):
    N = len(seqs)
    L = len(seqs[0])
    cov = []
    for j in range(L):
        n = sum(1 for s in seqs if s[j] not in GAPS)
        cov.append(n / N)
    return cov


def maior_corrida(indices):
    melhor, atual = [], []
    for j in indices:
        if atual and j == atual[-1] + 1:
            atual.append(j)
        else:
            atual = [j]
        if len(atual) > len(melhor):
            melhor = list(atual)
    return melhor


def recortar(regs, min_col, min_seq, maior_bloco=False):
    seqs = [s for _, s in regs]
    cov = cobertura_colunas(seqs)
    cols = [j for j, c in enumerate(cov) if c >= min_col - EPS]
    if maior_bloco and cols:
        cols = maior_corrida(cols)
    if not cols:
        return [], [], []
    cov_seq = [sum(1 for j in cols if s[j] not in GAPS) / len(cols) for s in seqs]
    manter = [i for i, c in enumerate(cov_seq) if c >= min_seq - EPS]
    # apaga colunas que ficaram so com gaps entre as sequencias mantidas
    cols = [j for j in cols if any(seqs[i][j] not in GAPS for i in manter)]
    return cols, manter, cov_seq


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--min-col", type=float, default=0.9)
    ap.add_argument("--min-seq", type=float, default=0.9)
    ap.add_argument("--maior-bloco", action="store_true")
    ap.add_argument("--varredura", action="store_true",
                    help="so imprime a tabela de sequencias/colunas que sobram em varios cortes")
    ap.add_argument("--saida")
    ap.add_argument("--relatorio")
    a = ap.parse_args()

    regs = ler_fasta(a.fasta)
    if not regs:
        sys.exit("FASTA vazio")
    L = {len(s) for _, s in regs}
    if len(L) != 1:
        sys.exit("sequencias com comprimentos diferentes - nao e um alinhamento")
    print("alinhamento: %d sequencias x %d colunas" % (len(regs), L.pop()))

    if a.varredura:
        cols_t = [0.7, 0.8, 0.9, 0.95, 1.0]
        seq_t = [0.8, 0.9, 0.95, 1.0]
        print("\nCada celula: colunas / sequencias que SOBRAM  (linhas = --min-col, colunas = --min-seq)")
        print("%-10s" % "min-col" + "".join("%16s" % ("seq>=%.2f" % t) for t in seq_t))
        for mc in cols_t:
            linha = "%-10.2f" % mc
            for ms in seq_t:
                cols, manter, _ = recortar(regs, mc, ms, a.maior_bloco)
                linha += "%16s" % ("%d / %d" % (len(cols), len(manter)))
            print(linha)
        print("\nEscolha o corte olhando esta tabela (nao o resultado da saturacao) e "
              "registre a escolha.")
        return

    if not a.saida:
        sys.exit("informe --saida (ou use --varredura)")
    cols, manter, cov_seq = recortar(regs, a.min_col, a.min_seq, a.maior_bloco)
    if not cols or not manter:
        sys.exit("nenhuma coluna/sequencia atende ao corte - afrouxe --min-col/--min-seq")
    seqs = [s for _, s in regs]
    novo = [(regs[i][0], "".join(seqs[i][j] for j in cols)) for i in manter]
    escrever_fasta(novo, a.saida)
    runs = 1 + sum(1 for x, y in zip(cols, cols[1:]) if y != x + 1)
    print("nucleo: %d colunas (em %d bloco(s) contiguo(s)), %d de %d sequencias -> %s" %
          (len(cols), runs, len(manter), len(regs), a.saida))
    print("intervalo de colunas (1-based): %d a %d" % (cols[0] + 1, cols[-1] + 1))
    fora = [regs[i][0].split()[0] for i in range(len(regs)) if i not in set(manter)]
    if fora:
        print("REMOVIDAS (%d): %s" % (len(fora), ", ".join(fora)))
    if a.relatorio:
        ms = set(manter)
        with open(a.relatorio, "w") as fh:
            fh.write("sequencia\tcobertura_do_nucleo\tstatus\tcabecalho\n")
            for i, (cab, _) in enumerate(regs):
                fh.write("%s\t%.3f\t%s\t%s\n" % (cab.split()[0], cov_seq[i],
                                                  "mantida" if i in ms else "removida", cab))
        print("relatorio por sequencia -> " + a.relatorio)
    print("Lembrete: nucleo menor = menos poder do teste; reporte tambem o alinhamento completo.")


if __name__ == "__main__":
    main()
