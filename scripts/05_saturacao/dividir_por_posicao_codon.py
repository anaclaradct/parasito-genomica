#!/usr/bin/env python3
"""
dividir_por_posicao_codon.py

Divide o alinhamento de COI em dois novos FASTAs -- 1a+2a posicao do
codon (juntas) e 3a posicao (sozinha) -- e calcula a proporcao de sitios
invariaveis de cada um, pronta pra colar no campo "Proportion of invariant
sites" do DAMBE (evita usar 0,676 errado, que e o valor do alinhamento
INTEIRO, nao de cada subconjunto por posicao).

O quadro de leitura (frame 0, 1 ou 2) e escolhido automaticamente pelo
mesmo criterio ja validado no checar_pseudogene_coi.py: o frame com MENOS
stop codons internos, testado numa sequencia de referencia do alinhamento.

Uso:
    python3 dividir_por_posicao_codon.py <alinhamento.fasta>

Gera:
    <nome>_pos12.fasta  (1a+2a posicao)
    <nome>_pos3.fasta   (3a posicao)
"""
import sys
import os

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


def count_internal_stops(seq_nogap, frame):
    seq = seq_nogap[frame:]
    n_codons = len(seq) // 3
    stops = 0
    for i in range(n_codons):
        codon = seq[i*3:i*3+3]
        if 'N' in codon or '-' in codon or len(codon) < 3:
            continue
        aa = CODON_TABLE_5.get(codon, 'X')
        if aa == '*' and i < n_codons - 2:
            stops += 1
    return stops


def detect_best_frame(records):
    """Usa a primeira sequencia sem gap como referencia, testa os 3 frames,
    escolhe o de menos stop codons internos (mesmo criterio ja validado)."""
    ref_seq = None
    for _, seq in records:
        if '-' not in seq:
            ref_seq = seq
            break
    if ref_seq is None:
        ref_seq = records[0][1].replace('-', '')

    best_frame, best_stops = 0, None
    for frame in (0, 1, 2):
        stops = count_internal_stops(ref_seq.replace('-', ''), frame)
        if best_stops is None or stops < best_stops:
            best_frame, best_stops = frame, stops
    return best_frame, best_stops


def split_positions(records, frame):
    """Retorna (records_pos12, records_pos3), preservando gaps (mantem o
    alinhamento -- posicoes de codon contadas sobre a sequencia COM gap,
    a partir do frame detectado)."""
    pos12, pos3 = [], []
    for header, seq in records:
        seq_frame = seq[frame:]
        p12_chars, p3_chars = [], []
        for i, ch in enumerate(seq_frame):
            codon_pos = i % 3
            if codon_pos in (0, 1):
                p12_chars.append(ch)
            else:
                p3_chars.append(ch)
        pos12.append((header, "".join(p12_chars)))
        pos3.append((header, "".join(p3_chars)))
    return pos12, pos3


def pct_invariant(records):
    if not records:
        return 0.0
    length = len(records[0][1])
    n_invariant = 0
    n_valid_cols = 0
    for i in range(length):
        col = [s[i] for _, s in records if i < len(s) and s[i] != '-']
        if not col:
            continue
        n_valid_cols += 1
        if len(set(col)) == 1:
            n_invariant += 1
    return (100 * n_invariant / n_valid_cols) if n_valid_cols else 0.0


def write_fasta(records, path):
    with open(path, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n{seq}\n")


def main():
    if len(sys.argv) < 2:
        print("Uso: python3 dividir_por_posicao_codon.py <alinhamento.fasta>")
        sys.exit(1)

    fasta_path = sys.argv[1]
    records = parse_fasta(fasta_path)

    frame, stops = detect_best_frame(records)
    print(f"Quadro de leitura detectado: frame {frame} (0-indexed) -- {stops} stop codon(s) interno(s) na sequencia de referencia")
    if stops > 0:
        print("  AVISO: mesmo o melhor frame teve stop codon interno na sequencia de referencia.")
        print("  Confira manualmente antes de confiar no resultado -- pode ser sequencia parcial")
        print("  comecando fora do inicio do codon (comum em fragmentos de PCR).")

    pos12, pos3 = split_positions(records, frame)

    base = os.path.splitext(fasta_path)[0]
    out12 = f"{base}_pos12.fasta"
    out3 = f"{base}_pos3.fasta"
    write_fasta(pos12, out12)
    write_fasta(pos3, out3)

    pct12 = pct_invariant(pos12)
    pct3 = pct_invariant(pos3)

    print(f"\nGerado: {out12}  ({len(pos12[0][1])} pb por sequencia)")
    print(f"  Proporcao de sitios invariaveis (1a+2a posicao): {pct12/100:.4f}  <- cole isso no DAMBE")
    print(f"\nGerado: {out3}  ({len(pos3[0][1])} pb por sequencia)")
    print(f"  Proporcao de sitios invariaveis (3a posicao):    {pct3/100:.4f}  <- cole isso no DAMBE")
    print(f"\nEsperado biologicamente: 3a posicao deve ter proporcao de invariaveis")
    print(f"MENOR que 1a+2a (ela muda mais livremente, sofre menos selecao purificadora).")
    if pct3 >= pct12:
        print("AVISO: o resultado saiu ao contrario do esperado -- vale conferir o frame detectado.")


if __name__ == "__main__":
    main()
