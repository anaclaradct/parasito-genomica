# parasito-genomica

Scripts de curadoria e análise de sequências de filárias do gênero *Mansonella*.

## O que este repositório contém

Etapas de curadoria e análise de marcadores mitocondriais (COI e 12S) de
*Mansonella* obtidos de bancos públicos (GenBank/NCBI): curadoria, alinhamento,
controle de qualidade, busca de NUMTs/pseudogenes, filogenia, saturação e
diversidade genética.

Todas as análises aqui partem de sequências públicas.

## Organização

    scripts/
      01_curadoria/          02_alinhamento/       03_numts/
      04_filogenia/          05_saturacao/         06_diversidade_redes/
      common/                (funções compartilhadas)

A descrição de cada script e a ordem de execução estão em
[`scripts/README.md`](scripts/README.md).

Os arquivos de dados (FASTA, alinhamentos, árvores) não são versionados;
ver `.gitignore`.

## Programas utilizados

Python 3, MAFFT (v7.525), IQ-TREE (v3.1.3), BLAST+ e Entrez Direct.

## Autora

Ana Clara da Costa Teixeira

## Licença

Ver o arquivo [`LICENSE`](LICENSE).
