# Scripts do pipeline — Mansonella

Os scripts das pastas `01` a `09` usam caminhos relativos à raiz do projeto
(a pasta que contém 01_CURATED, 03_ALIGNMENTS etc.). Rode sempre a partir dela, ex.:

    bash scripts/02_alinhamento/01_alinhar_COI.sh
    python3 scripts/02_alinhamento/qc_alinhamento.py 03_ALIGNMENTS/COI_aligned.fasta

## Marcadores mitocondriais e mitogenomas

| Etapa | Scripts |
|---|---|
| `01_curadoria` | `reorganizar_curadoria.py` → `gerar_metadados_finais.sh`, `padronizar_paises_12S.py` |
| `02_alinhamento` | `01_alinhar_COI.sh`, `qc_alinhamento.py`, `remover_outliers_gap.py`, `gerar_nucleo_comum_12S.py` |
| `03_numts` | `buscar_numts.py`, `checar_blast_suspeitas.sh` (chama `extrair_suspeitas.py`) |
| `04_filogenia` | `02_arvore_COI.sh`, `03_adicionar_outgroup_COI.sh`, `04_testar_outgroup2_COI.sh`, `gerar_itol_paises.sh`, `gerar_itol_paises_12S_en.py` |
| `05_saturacao` | `dividir_por_posicao_codon.py`, `preparar_dambe_12S.py` |
| `06_diversidade_redes` | `diversidade_coi.py`, `analises_complementares_coi.py`, `fasta_para_popart.py`, `Contar_diferencas.py`, `pares_proximos.py` |
| `07_mitogenomas` | `01_inventario.sh`, `extrair_genes_mito.py`, `checar_frameshift_cox1.py`, `traduzir_3_quadros.py`, `traduzir_para_alinhar.py`, `retrotraduzir.py`, `concatenar_genes.py`, `gerar_rotulos_mitogenomas.py`, `gerar_itol_mitogenomas.py` |
| `08_comparacao_topologias` | `rf_distancia.py` |
| `09_complementos` | `filo_utils.py` (funções compartilhadas), `dividir_por_particoes.py`, `metricas_por_gene.py`, `resolucao_arvores.py`, `nucleo_comum_alinhamento.py`, `podar_para_dambe.py` |
| `common` | `fasta_utils.py` (leitura de FASTA e tabela genética 5) |

## Assemblies nucleares

| Etapa | Scripts |
|---|---|
| `10_eixo2_assemblies` | `n50.py` (contiguidade), `analise_busco.py` (completude), `classifica.py` (classificação dos genes candidatos a partir do miniprot), `testa62.py` (busca por HMMER dos genes BUSCO ausentes) |
