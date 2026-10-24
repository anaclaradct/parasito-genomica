# parasito-genomica

Scripts de curadoria e análise de sequências de filárias do gênero *Mansonella*.

## O que este repositório contém

**Marcadores mitocondriais e mitogenomas** (pastas `01` a `09`): curadoria de
sequências de COI e 12S, alinhamento, controle de qualidade, busca de
NUMTs/pseudogenes, filogenia, saturação, diversidade genética, análise de
mitogenomas completos e comparação de topologias.

**Assemblies nucleares** (pasta `10`): contiguidade (N50), completude por genes
BUSCO e busca de genes candidatos por alinhamento proteína-genoma.

Todas as análises aqui partem de dados públicos (GenBank/NCBI).

## Organização

    scripts/
      01_curadoria/          02_alinhamento/            03_numts/
      04_filogenia/          05_saturacao/              06_diversidade_redes/
      07_mitogenomas/        08_comparacao_topologias/  09_complementos/
      10_eixo2_assemblies/   common/  (funções compartilhadas)

A descrição de cada script está em [`scripts/README.md`](scripts/README.md).

Os arquivos de dados (FASTA, alinhamentos, árvores, genomas) não são
versionados; ver `.gitignore`.

## Programas utilizados

Python 3, MAFFT (v7.525), IQ-TREE (v3.1.3), BLAST+, Entrez Direct, compleasm,
miniprot e HMMER.

## Autora

Ana Clara da Costa Teixeira

## Licença

Ver o arquivo [`LICENSE`](LICENSE).
