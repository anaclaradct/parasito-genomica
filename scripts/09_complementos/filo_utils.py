#!/usr/bin/env python3
"""
filo_utils.py - funcoes comuns aos scripts de complemento do Eixo 1
(FASTA, Newick, normalizacao de IDs, splits e distancia de Robinson-Foulds).
Sem dependencias externas (so biblioteca padrao do Python 3).
"""
import re

# accession sem versao e sem sufixo de coordenadas: KX822021.1:9215-9879 -> KX822021
ACC_RE = r"[A-Z]{1,2}_?\d{5,8}"


# ----------------------------------------------------------------- FASTA
def ler_fasta(caminho):
    """Retorna lista de (cabecalho_sem_>, sequencia_maiuscula)."""
    regs, cab, seq = [], None, []
    with open(caminho) as fh:
        for linha in fh:
            linha = linha.rstrip("\n\r")
            if linha.startswith(">"):
                if cab is not None:
                    regs.append((cab, "".join(seq).upper()))
                cab, seq = linha[1:].strip(), []
            elif linha.strip():
                seq.append(linha.strip())
    if cab is not None:
        regs.append((cab, "".join(seq).upper()))
    return regs


def escrever_fasta(regs, caminho, largura=80):
    with open(caminho, "w") as fh:
        for cab, seq in regs:
            fh.write(">" + cab + "\n")
            for i in range(0, len(seq), largura):
                fh.write(seq[i:i + largura] + "\n")


def acc_key(nome, regex=ACC_RE):
    """Extrai o accession (sem versao/coordenadas). Se nao casar, devolve o nome."""
    m = re.search(regex, nome)
    return m.group(0) if m else nome


# ---------------------------------------------------------------- Newick
class No:
    __slots__ = ("name", "length", "label", "children", "parent")

    def __init__(self, parent=None):
        self.name = None
        self.length = None
        self.label = None
        self.children = []
        self.parent = parent


def ler_newick(texto):
    s = re.sub(r"\[[^\]]*\]", "", texto).strip()
    n = len(s)

    def ler_rotulo(i, paradas):
        if i < n and s[i] == "'":
            j = s.index("'", i + 1)
            return s[i + 1:j], j + 1
        j = i
        while j < n and s[j] not in paradas:
            j += 1
        return s[i:j].strip(), j

    raiz, cur, ultimo, i = None, None, None, 0
    while i < n:
        c = s[i]
        if c == "(":
            if cur is None:             # primeiro parentese = raiz
                cur = No(None)
                raiz = cur
            slot = No(cur)              # o no atual vira interno; abre o slot do 1o filho
            cur.children.append(slot)
            cur, ultimo, i = slot, None, i + 1
        elif c == ",":
            pai = cur.parent
            slot = No(pai)
            pai.children.append(slot)
            cur, ultimo, i = slot, None, i + 1
        elif c == ")":
            cur = cur.parent            # volta ao no interno que acabou de fechar
            i += 1
            rot, i = ler_rotulo(i, ":,();")
            cur.label = rot or None
            ultimo = cur
        elif c == ":":
            j = i + 1
            while j < n and s[j] not in ",);":
                j += 1
            alvo = ultimo if ultimo is not None else cur
            try:
                alvo.length = float(s[i + 1:j])
            except ValueError:
                alvo.length = None
            i = j
        elif c == ";":
            break
        elif c.isspace():
            i += 1
        else:
            nome, i = ler_rotulo(i, ":,();")
            cur.name = nome
            ultimo = cur
    if raiz is None:
        raise ValueError("Newick sem parenteses")
    raiz.parent = None
    return raiz


def ler_arvore_arquivo(caminho, linha=0):
    with open(caminho) as fh:
        linhas = [l for l in fh.read().split("\n") if l.strip()]
    return ler_newick(linhas[linha])


def pos_ordem(raiz):
    saida, pilha = [], [(raiz, False)]
    while pilha:
        no, feito = pilha.pop()
        if feito:
            saida.append(no)
        else:
            pilha.append((no, True))
            for c in no.children:
                pilha.append((c, False))
    return saida


def folhas(raiz):
    return [n.name for n in pos_ordem(raiz) if not n.children]


def renomear_folhas(raiz, fn):
    vistos = set()
    for n in pos_ordem(raiz):
        if not n.children:
            n.name = fn(n.name)
            if n.name in vistos:
                raise ValueError("nome duplicado apos normalizacao: " + n.name)
            vistos.add(n.name)


def _fmt_nome(nome):
    if re.search(r"[\s(),:;\[\]']", nome):
        return "'" + nome + "'"
    return nome


def para_newick(raiz):
    def rec(x):
        if x.children:
            s = "(" + ",".join(rec(c) for c in x.children) + ")"
            if x.label:
                s += x.label
        else:
            s = _fmt_nome(x.name)
        if x.length is not None:
            s += ":%.10g" % x.length
        return s
    return rec(raiz) + ";"


def podar_arvore(raiz, remover):
    """Remove folhas cujo nome satisfaz remover(nome); colapsa nos unarios somando comprimentos."""
    def rec(x):
        if not x.children:
            return None if remover(x.name) else x
        filhos = [k for k in (rec(c) for c in x.children) if k is not None]
        if not filhos:
            return None
        x.children = filhos
        for k in filhos:
            k.parent = x
        if len(filhos) == 1:
            f = filhos[0]
            if f.length is not None or x.length is not None:
                f.length = (f.length or 0.0) + (x.length or 0.0)
            return f
        return x
    nova = rec(raiz)
    if nova is not None:
        nova.parent = None
    return nova


# --------------------------------------------------- suporte e splits
def valores_suporte(rotulo):
    """'84.5/97' -> [84.5, 97.0]; sem rotulo -> None."""
    if not rotulo:
        return None
    try:
        return [float(p) for p in re.split(r"[/,]", rotulo)]
    except ValueError:
        return None


def passa(vals, sh_min=None, ufb_min=None):
    """True/False se o suporte atinge os limiares; None se nao ha informacao."""
    if vals is None:
        return None
    if len(vals) >= 2:
        if sh_min is not None and vals[0] < sh_min:
            return False
        if ufb_min is not None and vals[1] < ufb_min:
            return False
        return True
    lim = ufb_min if ufb_min is not None else sh_min
    return True if lim is None else vals[0] >= lim


def _fundir(a, b):
    if a is None:
        return b
    if b is None:
        return a
    if len(a) == len(b):
        return [min(x, y) for x, y in zip(a, b)]  # conservador
    return a


def obter_splits(raiz, manter=None, min_ramo=0.0):
    """
    Bipartições nao triviais (arvore tratada como nao enraizada).
    Retorna (dict {frozenset(lado sem o menor taxon): suporte|None}, frozenset(taxons)).
    manter: conjunto de taxons a considerar (restringe a arvore SEM recalcular suporte;
    se duas arestas colapsam no mesmo split fica o MENOR suporte).
    min_ramo: ramos internos com comprimento < min_ramo sao tratados como colapsados
    (politomia). Necessario porque o IQ-TREE resolve arbitrariamente ramos de
    comprimento ~0 (minimo ~1e-6); use 0 para desligar.
    """
    ordem = pos_ordem(raiz)
    conj = {}
    for n in ordem:
        if not n.children:
            ok = manter is None or n.name in manter
            conj[n] = frozenset([n.name]) if ok else frozenset()
        else:
            s = set()
            for c in n.children:
                s |= conj[c]
            conj[n] = frozenset(s)
    todos = conj[raiz]
    N = len(todos)
    if N == 0:
        return {}, todos
    ref = min(todos)
    saida = {}
    # raiz bifurcante com dois filhos internos: as duas arestas formam UMA so (a raiz e artificial)
    raiz_bif = len(raiz.children) == 2 and all(c.children for c in raiz.children)
    for n in ordem:
        if n is raiz or not n.children:
            continue
        if raiz_bif and n.parent is raiz:
            continue
        if n.length is not None and n.length < min_ramo:
            continue
        s = conj[n]
        if len(s) < 2 or len(s) > N - 2:
            continue
        if ref in s:
            s = todos - s
        v = valores_suporte(n.label)
        saida[s] = _fundir(saida[s], v) if s in saida else v
    if raiz_bif:
        c1, c2 = raiz.children
        s = conj[c1]
        comp = None if (c1.length is None and c2.length is None) else (c1.length or 0.0) + (c2.length or 0.0)
        if comp is not None and comp < min_ramo:
            s = frozenset()
        if 2 <= len(s) <= N - 2:
            if ref in s:
                s = todos - s
            v1, v2 = valores_suporte(c1.label), valores_suporte(c2.label)
            if v1 is None or v2 is None:
                v = v1 if v2 is None else v2
            elif len(v1) == len(v2):
                v = [max(x, y) for x, y in zip(v1, v2)]
            else:
                v = v1
            saida[s] = _fundir(saida[s], v) if s in saida else v
    return saida, todos


def filtrar_suportados(splits, sh_min=None, ufb_min=None):
    return {s: v for s, v in splits.items() if passa(v, sh_min, ufb_min) is True}


def rf(a, b, N):
    """RF nao enraizado. Retorna dict com bruto, maximo 2(N-3), normalizado, compartilhados."""
    A, B = set(a), set(b)
    d = len(A - B) + len(B - A)
    mx = 2 * (N - 3) if N >= 4 else 0
    return {"rf": d, "rf_max": mx, "rf_norm": (d / mx if mx else float("nan")),
            "compartilhados": len(A & B), "n_a": len(A), "n_b": len(B)}


def restringir_split(lado, todos_ref, T):
    """Restringe um split (lado sem o menor taxon de todos_ref) ao conjunto T.
    Devolve o split canonico em T ou None se ficar trivial."""
    s = frozenset(lado) & T
    if len(s) < 2 or len(T - s) < 2:
        return None
    if min(T) in s:
        s = frozenset(T - s)
    return s


def carregar_mapa(caminho):
    """TSV nome_original<TAB>id_comum"""
    m = {}
    with open(caminho) as fh:
        for l in fh:
            p = l.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0] and not p[0].startswith("#"):
                m[p[0]] = p[1]
    return m


def fabrica_normalizador(regex=ACC_RE, mapa=None):
    def fn(nome):
        if mapa and nome in mapa:
            return mapa[nome]
        return acc_key(nome, regex)
    return fn
