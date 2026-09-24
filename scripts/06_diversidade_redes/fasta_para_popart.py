#!/usr/bin/env python3
"""
fasta_para_popart.py
Converte um alinhamento FASTA + metadados em um arquivo NEXUS pronto pra
abrir direto no PopART (File > Open), ja com o bloco TRAITS para colorir
a rede por grupo (especie, pais, etc).

Uso:
    python3 fasta_para_popart.py <alinhamento.fasta> <metadados.tsv> <coluna_grupo> <saida.nex> [--filtro coluna=valor] [--core]

Exemplos:
    # Rede com TODAS as sequencias, coloridas por especie
    python3 fasta_para_popart.py 03_ALIGNMENTS/COI_aligned.fasta metadados_finais.tsv especie_pasta rede_todas_especies.nex

    # Rede SO de M. perstans, coloridas por pais, usando a regiao core
    # definida por delecao completa. Mantem todas as sequencias originais;
    # sequencias identicas na regiao core sao representadas pelo mesmo
    # haplotipo no PopART.
    python3 fasta_para_popart.py 03_ALIGNMENTS/COI_aligned.fasta metadados_finais.tsv country rede_perstans_por_pais_core.nex --filtro especie_pasta=Mansonella_perstans --core

NOTA: sem --core, o script manda os 689 caracteres completos (com gap) pro
PopART -- util como exploracao inicial ("rede exploratoria"), mas os
gaps das sequencias mais curtas podem gerar nos artificiais e o numero de
"haplotipos" da rede pode nao bater com o numero reportado no diversidade_coi.py
(que usa deleção completa). Use --core para a versao comparavel/definitiva.
"""
import sys
import csv

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

def complete_deletion_core(seqs):
    """Mesma logica do diversidade_coi.py: mantem so colunas sem NENHUM gap
    em NENHUMA sequencia do conjunto atual (apos filtro de especie)."""
    if not seqs:
        return seqs, 0, []
    length = len(seqs[0])
    keep_cols = [i for i in range(length) if all(s[i] != "-" for s in seqs)]
    core_seqs = ["".join(s[i] for i in keep_cols) for s in seqs]
    return core_seqs, len(keep_cols), keep_cols

def main():
    if len(sys.argv) < 5:
        print("Uso: python3 fasta_para_popart.py <alinhamento.fasta> <metadados.tsv> <coluna_grupo> <saida.nex> [--filtro coluna=valor] [--core]")
        sys.exit(1)

    fasta_path, meta_path, group_col, out_path = sys.argv[1:5]
    resto = sys.argv[5:]
    filtro = None
    use_core = "--core" in resto
    resto = [a for a in resto if a != "--core"]
    if len(resto) >= 2 and resto[0] == "--filtro":
        fcol, fval = resto[1].split("=", 1)
        filtro = (fcol, fval)

    records = parse_fasta(fasta_path)
    meta = {}
    with open(meta_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            meta[row["accession"]] = row

    # aplica filtro (ex: so uma especie) e pega o grupo de cada sequencia mantida
    kept = []
    for acc, seq in records:
        row = meta.get(acc)
        if row is None:
            continue
        if filtro and row.get(filtro[0]) != filtro[1]:
            continue
        group_raw = row.get(group_col, "NA")
        # simplifica pais tipo "Peru: Loreto" -> "Peru"
        group = group_raw.split(":")[0].strip() if group_raw else "NA"
        if not group or group == "NA":
            group = "Sem_dado"
        group = group.replace(" ", "_")
        kept.append((acc, seq, group))

    if not kept:
        print("ERRO: nenhuma sequencia sobrou apos aplicar metadados/filtro. Confira os nomes de coluna/valores.")
        sys.exit(1)

    lengths = set(len(s) for _, s, _ in kept)
    if len(lengths) > 1:
        print(f"AVISO: sequencias com comprimentos diferentes {sorted(lengths)} -- "
              f"confirme que o arquivo de entrada e mesmo um alinhamento (mafft), nao fasta cru.")

    if use_core:
        seqs_only = [s for _, s, _ in kept]
        core_seqs, core_len, keep_cols = complete_deletion_core(seqs_only)
        kept = [(acc, core_seqs[i], g) for i, (acc, _, g) in enumerate(kept)]
        n_unique_core = len(set(core_seqs))
        print(f">> Modo --core ativado: alinhamento recortado para {core_len} pb "
              f"(deleção completa, sem gap em nenhuma sequência deste subconjunto)")
        print(f"   Sequências idênticas nessa região formam {n_unique_core} grupos distintos "
              f"(equivalente ao número de haplótipos que o PopART deve mostrar; "
              f"confira se bate com o 'core próprio' do diversidade_coi.py para este grupo)")

    nchar = max(len(s) for _, s, _ in kept)

    groups = sorted(set(g for _, _, g in kept))
    n_taxa = len(kept)

    with open(out_path, "w") as f:
        f.write("#NEXUS\n\n")
        f.write("BEGIN TAXA;\n")
        f.write(f"DIMENSIONS NTAX={n_taxa};\n")
        f.write("TAXLABELS\n")
        for acc, _, _ in kept:
            f.write(f"{acc}\n")
        f.write(";\nEND;\n\n")

        f.write("BEGIN CHARACTERS;\n")
        f.write(f"DIMENSIONS NCHAR={nchar};\n")
        f.write("FORMAT DATATYPE=DNA MISSING=? GAP=- ;\n")
        f.write("MATRIX\n")
        for acc, seq, _ in kept:
            f.write(f"{acc}  {seq}\n")
        f.write(";\nEND;\n\n")

        f.write("BEGIN TRAITS;\n")
        f.write(f"Dimensions NTRAITS={len(groups)};\n")
        f.write("Format labels=yes missing=? separator=Comma;\n")
        f.write("TraitLabels " + " ".join(groups) + ";\n")
        f.write("Matrix\n")
        for acc, _, g in kept:
            row = ",".join("1" if g == gg else "0" for gg in groups)
            f.write(f"{acc} {row}\n")
        f.write(";\nEND;\n")

    print(f"Gerado: {out_path}")
    print(f"  Sequencias: {n_taxa}")
    print(f"  Grupos ({group_col}): {groups}")

if __name__ == "__main__":
    main()
