#!/usr/bin/env python3
"""Compara cada sequencia de COI com o cox1 de mitogenomas completos (referencias).
Uso: python3 checar_frameshift_cox1.py alinhamento.fasta REF1,REF2,... > saida.tsv
Para cada sequencia nao-referencia, informa:
  - referencia mais proxima e distancia p (so posicoes sem gap nas duas)
  - indels internos em relacao a referencia mais proxima cujo tamanho NAO e multiplo de 3
    (um indel assim num gene codificador quebra o quadro de leitura = sinal de pseudogene/NUMT
     ou de erro de sequenciamento; o cox1 do mitogenoma e codificador validado)
Sem dependencias externas.
"""
import sys

def ler_fasta(caminho):
    seqs, nome = {}, None
    for linha in open(caminho):
        linha = linha.strip()
        if linha.startswith(">"):
            nome = linha[1:].split()[0]
            seqs[nome] = []
        elif nome:
            seqs[nome].append(linha.upper())
    return {k: "".join(v) for k, v in seqs.items()}

def comparar(a, b):
    """a = consulta, b = referencia (ja alinhadas). Retorna p, n_sitios, indels_fora_de_quadro"""
    cob = [i for i in range(len(a)) if a[i] != "-" and b[i] != "-"]
    if not cob:
        return None, 0, []
    ini, fim = cob[0], cob[-1]
    dif = sum(1 for i in cob if a[i] != b[i])
    p = dif / len(cob)
    # indels internos: blocos de coluna onde exatamente uma das duas tem gap
    indels, i = [], ini
    while i <= fim:
        if (a[i] == "-") != (b[i] == "-"):
            quem = "del_na_consulta" if a[i] == "-" else "ins_na_consulta"
            j = i
            while j <= fim and (a[j] == "-") != (b[j] == "-") and \
                    ((a[j] == "-") == (quem == "del_na_consulta")):
                j += 1
            tam = j - i
            if tam % 3 != 0:
                indels.append("%s:%dpb@col%d" % (quem, tam, i + 1))
            i = j
        else:
            i += 1
    return p, len(cob), indels

def main(aln, refs):
    seqs = ler_fasta(aln)
    refs = [r for r in refs.split(",") if r]
    faltam = [r for r in refs if r not in seqs]
    if faltam:
        sys.exit("referencias nao encontradas no alinhamento: " + ", ".join(faltam))
    print("sequencia\tref_mais_proxima\tp_dist\tsitios_comparados\tn_indels_fora_de_quadro\tdetalhe")
    for nome, s in seqs.items():
        if nome in refs:
            continue
        melhor = None
        for r in refs:
            p, n, ind = comparar(s, seqs[r])
            if p is not None and n >= 100 and (melhor is None or p < melhor[1]):
                melhor = (r, p, n, ind)
        if melhor is None:
            print("%s\tNA\tNA\t0\tNA\tsobreposicao<100pb" % nome)
            continue
        r, p, n, ind = melhor
        print("%s\t%s\t%.4f\t%d\t%d\t%s" % (nome, r, p, n, len(ind), ";".join(ind) or "-"))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
