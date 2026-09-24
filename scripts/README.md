# Scripts do pipeline — Mansonella (COI / 12S)

Todos os scripts usam caminhos relativos à raiz do projeto: rode sempre a
partir de `~/Mansonella_projeto`, ex.:

    bash scripts/02_alinhamento/01_alinhar_COI.sh
    python3 scripts/02_alinhamento/qc_alinhamento.py 03_ALIGNMENTS/COI_aligned.fasta

| Etapa | Scripts |
|---|---|
| `01_curadoria` | `reorganizar_curadoria.py` → `gerar_metadados_finais.sh` |
| `02_alinhamento` | `01_alinhar_COI.sh`, `qc_alinhamento.py`, `remover_outliers_gap.py`, `gerar_nucleo_comum_12S.py` |
| `03_numts` | `buscar_numts.py`, `checar_blast_suspeitas.sh` (chama `extrair_suspeitas.py`) |
| `04_filogenia` | `02_arvore_COI.sh`, `03_adicionar_outgroup_COI.sh`, `04_testar_outgroup2_COI.sh`, `gerar_itol_paises.sh` |
| `05_saturacao` | `dividir_por_posicao_codon.py`, `preparar_dambe_12S.py` |
| `06_diversidade_redes` | `diversidade_coi.py`, `analises_complementares_coi.py`, `fasta_para_popart.py`, `Contar_diferencas.py` |
