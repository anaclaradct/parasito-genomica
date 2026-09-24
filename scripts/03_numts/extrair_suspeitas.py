targets = {"LT623910.1", "KP760195.1"}
records = {}
header, seq = None, []

with open("03_ALIGNMENTS/COI_aligned.fasta") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if header is not None:
                records[header] = "".join(seq)
            header = line[1:].split()[0]
            seq = []
        else:
            seq.append(line.strip())
    if header is not None:
        records[header] = "".join(seq)

with open("numt_check/suspeitas_blast.fasta", "w") as out:
    for acc in targets:
        if acc in records:
            seq_nogap = records[acc].replace("-", "")
            out.write(f">{acc}\n{seq_nogap}\n")
            print(f"   {acc}: {len(seq_nogap)} pb extraidos")
        else:
            print(f"   AVISO: {acc} nao encontrado no alinhamento!")
