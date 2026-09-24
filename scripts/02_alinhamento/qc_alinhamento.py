#!/usr/bin/env python3
"""
qc_alinhamento.py
Gera um relatorio de qualidade de um alinhamento FASTA -- sem depender de
trimal, AliView ou qualquer coisa alem do Python padrao (ja instalado).

Uso:
    python3 qc_alinhamento.py 03_ALIGNMENTS/COI_aligned.fasta

Gera, ao lado do arquivo de entrada:
    <nome>_qc.tsv    -- uma linha por sequencia, com % de gap, gaps nas
                        pontas (terminal), etc.
    <nome>_qc.html   -- relatorio visual (abre no navegador), com barras
                        coloridas por sequencia e um mapa de gaps por coluna

O que observar no relatorio:
    - Sequencias com % de gap muito acima da media -> provavelmente
      cobrem uma regiao menor do gene, ou tem baixa qualidade.
    - Gaps nas PONTAS (terminal_gap_5/3) grandes -> normal quando a
      sequencia so cobre parte do marcador (comum com sequencias parciais
      do GenBank de diferentes primers); gaps no MEIO da sequencia sao
      mais preocupantes (podem indicar erro de alinhamento ou indel real).
    - Colunas do alinhamento com >50% de gap no mapa de calor -> regioes
      pouco informativas, candidatas a trimming antes do IQ-TREE2.
"""
import sys
import os

def parse_fasta(path):
    records = []
    header, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line.strip())
        if header is not None:
            records.append((header, "".join(seq)))
    return records

def gap_stats(seq):
    length = len(seq)
    n_gap = seq.count("-")
    pct_gap = 100 * n_gap / length if length else 0
    lead = len(seq) - len(seq.lstrip("-"))
    trail = len(seq) - len(seq.rstrip("-"))
    return n_gap, pct_gap, lead, trail

def main():
    if len(sys.argv) < 2:
        print("Uso: python3 qc_alinhamento.py <alinhamento.fasta>")
        sys.exit(1)

    path = sys.argv[1]
    records = parse_fasta(path)
    if not records:
        print("Nenhuma sequencia encontrada no arquivo.")
        sys.exit(1)

    lengths = set(len(seq) for _, seq in records)
    if len(lengths) > 1:
        print(f"AVISO: sequencias com comprimentos diferentes ({sorted(lengths)}) "
              f"-- isso NAO deveria acontecer num alinhamento. Confira se o "
              f"arquivo e realmente o output do MAFFT.")
    aln_len = max(lengths)

    base = os.path.splitext(path)[0]
    tsv_path = base + "_qc.tsv"
    html_path = base + "_qc.html"

    rows = []
    for header, seq in records:
        acc = header.split()[0]
        n_gap, pct_gap, lead, trail = gap_stats(seq)
        internal_gap_pct = 100 * (n_gap - lead - trail) / aln_len if aln_len else 0
        rows.append({
            "accession": acc, "header": header, "seq": seq,
            "length_no_gap": len(seq) - n_gap, "pct_gap": pct_gap,
            "terminal_gap_5": lead, "terminal_gap_3": trail,
            "internal_gap_pct": internal_gap_pct,
        })

    mean_gap = sum(r["pct_gap"] for r in rows) / len(rows)
    outlier_thresh = mean_gap + 20  # 20 pontos percentuais acima da media

    with open(tsv_path, "w") as f:
        f.write("accession\tcomprimento_sem_gap\tpct_gap_total\tgap_5prime\tgap_3prime\tpct_gap_interno\tflag\n")
        for r in sorted(rows, key=lambda x: -x["pct_gap"]):
            flag = ""
            if r["pct_gap"] > outlier_thresh:
                flag = "OUTLIER_gap_total"
            elif r["internal_gap_pct"] > 5:
                flag = "gap_interno_relevante"
            f.write(f"{r['accession']}\t{r['length_no_gap']}\t{r['pct_gap']:.1f}\t"
                    f"{r['terminal_gap_5']}\t{r['terminal_gap_3']}\t{r['internal_gap_pct']:.1f}\t{flag}\n")

    # coluna a coluna: % de gap por posicao (para o mapa de calor)
    col_gap_pct = []
    for i in range(aln_len):
        n = sum(1 for r in rows if i < len(r["seq"]) and r["seq"][i] == "-")
        col_gap_pct.append(100 * n / len(rows))

    n_outliers = sum(1 for r in rows if r["pct_gap"] > outlier_thresh)
    n_internal = sum(1 for r in rows if r["internal_gap_pct"] > 5)

    # HTML
    html = []
    html.append("<html><head><meta charset='utf-8'><title>QC Alinhamento</title>")
    html.append("""<style>
    body{font-family:Arial,sans-serif;margin:24px;background:#fafafa}
    h1{font-size:20px} h2{font-size:16px;margin-top:28px}
    table{border-collapse:collapse;width:100%;font-size:12px}
    th,td{border:1px solid #ddd;padding:4px 8px;text-align:left}
    th{background:#1F4E78;color:white;position:sticky;top:0}
    tr:nth-child(even){background:#f2f2f2}
    .bar-container{background:#eee;width:200px;height:12px;display:inline-block;vertical-align:middle}
    .bar{background:#c0392b;height:12px}
    .flag{color:#c0392b;font-weight:bold}
    .heatmap{display:flex;height:24px;border:1px solid #ccc;margin-top:8px}
    .col{flex:1;min-width:1px}
    .summary{background:white;border:1px solid #ddd;padding:12px;border-radius:6px;margin-bottom:16px}
    </style></head><body>""")
    html.append(f"<h1>QC do alinhamento: {os.path.basename(path)}</h1>")
    html.append(f"<div class='summary'>")
    html.append(f"<b>Sequencias:</b> {len(rows)} &nbsp; | &nbsp; <b>Comprimento do alinhamento:</b> {aln_len} pb<br>")
    html.append(f"<b>% de gap medio por sequencia:</b> {mean_gap:.1f}% &nbsp; | &nbsp; ")
    html.append(f"<b>Outliers de gap total (&gt;{outlier_thresh:.0f}%):</b> {n_outliers} &nbsp; | &nbsp; ")
    html.append(f"<b>Com gap interno relevante (&gt;5%):</b> {n_internal}")
    html.append("</div>")

    html.append("<h2>Mapa de gap por coluna do alinhamento (vermelho = mais gap nessa posicao)</h2>")
    html.append("<div class='heatmap'>")
    for i, pct in enumerate(col_gap_pct):
        g = int(255 * (1 - pct / 100))
        html.append(f"<div class='col' style='background:rgb(255,{g},{g})' title='pos {i+1}: {pct:.0f}% gap'></div>")
    html.append("</div>")

    html.append("<h2>Por sequencia (ordenado por % de gap total, maior primeiro)</h2>")
    html.append("<table><tr><th>Accession</th><th>Header</th><th>Comprimento sem gap</th>"
                 "<th>% gap total</th><th>Gap 5'</th><th>Gap 3'</th><th>% gap interno</th><th>Flag</th></tr>")
    for r in sorted(rows, key=lambda x: -x["pct_gap"]):
        flag = ""
        if r["pct_gap"] > outlier_thresh:
            flag = "OUTLIER (gap total)"
        elif r["internal_gap_pct"] > 5:
            flag = "gap interno"
        bar_width = min(200, int(2 * r["pct_gap"]))
        html.append(f"<tr><td>{r['accession']}</td><td>{r['header'][:70]}</td>"
                     f"<td>{r['length_no_gap']}</td>"
                     f"<td><div class='bar-container'><div class='bar' style='width:{bar_width}px'></div></div> {r['pct_gap']:.1f}%</td>"
                     f"<td>{r['terminal_gap_5']}</td><td>{r['terminal_gap_3']}</td>"
                     f"<td>{r['internal_gap_pct']:.1f}%</td>"
                     f"<td class='flag'>{flag}</td></tr>")
    html.append("</table></body></html>")

    with open(html_path, "w") as f:
        f.write("\n".join(html))

    print(f"Sequencias: {len(rows)} | Comprimento do alinhamento: {aln_len} pb")
    print(f"% de gap medio: {mean_gap:.1f}%")
    print(f"Outliers de gap total (>{outlier_thresh:.0f}%): {n_outliers}")
    print(f"Com gap interno relevante (>5%): {n_internal}")
    print(f"Relatorio TSV:  {tsv_path}")
    print(f"Relatorio HTML: {html_path}")

if __name__ == "__main__":
    main()
