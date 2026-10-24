#!/usr/bin/env python3
"""Calcula a distancia de Robinson-Foulds entre duas arvores Newick, restrita a um
conjunto comum de taxons (com possibilidade de renomear rotulos antes de comparar).
Uso: rf_distancia.py arvoreA.nwk arvoreB.nwk mapaA.tsv mapaB.tsv
Cada mapaX.tsv tem duas colunas (rotulo_original<TAB>rotulo_comum). Só as folhas
listadas no mapa sao mantidas; as demais sao podadas. Os rotulos comuns devem
bater entre os dois mapas para as folhas que devem ser comparadas.
Sem dependencias externas (sem ete3/dendropy/biopython).
"""
import sys

def parse_newick(s):
    s = s.strip()
    if s.endswith(";"):
        s = s[:-1]
    pos = [0]
    def parse_clade():
        node = {"children": [], "name": None}
        if s[pos[0]] == "(":
            pos[0] += 1
            while True:
                node["children"].append(parse_clade())
                if s[pos[0]] == ",":
                    pos[0] += 1
                    continue
                elif s[pos[0]] == ")":
                    pos[0] += 1
                    break
            # nome interno opcional apos ')'
            start = pos[0]
            while pos[0] < len(s) and s[pos[0]] not in ",():;":
                pos[0] += 1
            node["name"] = s[start:pos[0]] or None
        else:
            start = pos[0]
            while pos[0] < len(s) and s[pos[0]] not in ",():;":
                pos[0] += 1
            node["name"] = s[start:pos[0]]
        # remove e consome o comprimento de ramo (":0.123"), se presente
        if node["name"] and ":" in node["name"]:
            node["name"] = node["name"].split(":")[0]
        if pos[0] < len(s) and s[pos[0]] == ":":
            pos[0] += 1
            while pos[0] < len(s) and s[pos[0]] not in ",();":
                pos[0] += 1
        return node
    return parse_clade()

def get_leaves(node):
    if not node["children"]:
        return [node["name"]] if node["name"] else []
    out = []
    for c in node["children"]:
        out.extend(get_leaves(c))
    return out

def relabel_and_prune(node, mapa):
    """Retorna None se a subarvore nao tem nenhuma folha do mapa, ou a subarvore podada."""
    if not node["children"]:
        novo = mapa.get(node["name"])
        if novo is None:
            return None
        return {"children": [], "name": novo}
    filhos = [relabel_and_prune(c, mapa) for c in node["children"]]
    filhos = [f for f in filhos if f is not None]
    if not filhos:
        return None
    if len(filhos) == 1:
        return filhos[0]  # colapsa no unario
    return {"children": filhos, "name": None}

def bipartitions(node, all_leaves, ref):
    """Retorna o conjunto de biparticoes nao triviais, cada uma normalizada como
    frozenset do lado que NAO contem a folha de referencia 'ref'."""
    biparts = set()
    def rec(n):
        if not n["children"]:
            return {n["name"]}
        leafset = set()
        for c in n["children"]:
            leafset |= rec(c)
        return leafset
    def walk(n, is_root):
        if not n["children"]:
            return
        leafset = rec(n)
        if not is_root and 2 <= len(leafset) <= len(all_leaves) - 2:
            lado = leafset if ref not in leafset else (all_leaves - leafset)
            biparts.add(frozenset(lado))
        for c in n["children"]:
            walk(c, False)
    walk(node, True)
    return biparts

def ler_mapa(caminho):
    mapa = {}
    for linha in open(caminho):
        linha = linha.rstrip("\n")
        if not linha or linha.startswith("#"):
            continue
        orig, novo = linha.split("\t")
        mapa[orig] = novo
    return mapa

def rf(arqA, arqB, mapaA_p, mapaB_p):
    mapaA = ler_mapa(mapaA_p)
    mapaB = ler_mapa(mapaB_p)
    comuns = sorted(set(mapaA.values()) & set(mapaB.values()))
    if len(comuns) < 4:
        sys.exit("ERRO: menos de 4 taxons comuns entre os dois mapas (%d) — RF nao definido" % len(comuns))
    mapaA = {k: v for k, v in mapaA.items() if v in comuns}
    mapaB = {k: v for k, v in mapaB.items() if v in comuns}
    treeA = relabel_and_prune(parse_newick(open(arqA).read()), mapaA)
    treeB = relabel_and_prune(parse_newick(open(arqB).read()), mapaB)
    leavesA = set(get_leaves(treeA))
    leavesB = set(get_leaves(treeB))
    if leavesA != set(comuns) or leavesB != set(comuns):
        sys.exit("ERRO: folhas apos poda nao batem com os taxons comuns esperados.\nA: %s\nB: %s\nesperado: %s" % (sorted(leavesA), sorted(leavesB), comuns))
    ref = comuns[0]
    bipsA = bipartitions(treeA, set(comuns), ref)
    bipsB = bipartitions(treeB, set(comuns), ref)
    rf_bruto = len(bipsA ^ bipsB)
    maximo = 2 * (len(comuns) - 3)  # duas arvores binarias totalmente resolvidas, sem taxons em comum nenhuma biparticao
    print("taxons comparados: %d" % len(comuns))
    print("biparticoes internas em A: %d" % len(bipsA))
    print("biparticoes internas em B: %d" % len(bipsB))
    print("RF (bruto, |diferenca simetrica|): %d" % rf_bruto)
    print("RF maximo possivel (arvores totalmente resolvidas): %d" % maximo)
    if maximo > 0:
        print("RF normalizado: %.3f" % (rf_bruto / maximo))
    so_A = bipsA - bipsB
    so_B = bipsB - bipsA
    if so_A:
        print("biparticoes so em A:")
        for b in sorted(so_A, key=lambda x: (len(x), sorted(x))):
            print("  ", sorted(b))
    if so_B:
        print("biparticoes so em B:")
        for b in sorted(so_B, key=lambda x: (len(x), sorted(x))):
            print("  ", sorted(b))

if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    rf(*sys.argv[1:5])
