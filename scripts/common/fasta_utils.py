"""Utilitarios de FASTA compartilhados entre os scripts do pipeline.

parse_fasta e CODON_TABLE_5 estavam duplicados (copia-e-cola identica) em
varios scripts de 03_numts/, 05_saturacao/ e 06_diversidade_redes/; ficam
centralizados aqui para que uma correcao futura (ex: lidar com CRLF ou
headers duplicados) se propague para todos os usos de uma vez.
"""


def parse_fasta(path):
    """Le um FASTA e retorna lista de (header_ate_primeiro_espaco, sequencia_maiuscula)."""
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


CODON_TABLE_5 = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'M', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': 'W', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'S', 'AGG': 'S',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}
