"""
froger_viz — logique pure (sans Streamlit)
==========================================
Construction d'un stemma à partir d'une sélection de groupes + dessins
matplotlib. Importé par app.py, mais testable seul.
"""
from __future__ import annotations
from pyfroger.group_method import UnorientedStemma
from pyfroger.list_apparatus import ListApparatus
import networkx as nx
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch
import matplotlib.pyplot as plt

from collections import defaultdict

import matplotlib
matplotlib.use("Agg")


PALETTE = {
    "paper":   "#faf7f2",
    "ink":     "#16181d",
    "indigo":  "#3b4ba8",
    "indigo2": "#9aa3d6",
    "rubric":  "#b5341f",
    "laminar": "#3f7d5c",
    "muted":   "#8a8378",
    "rule":    "#e3dccf",
}

HYP_CHARS = set("αβγδεζηθικλμνξοπρστυφχψωאבגדהוזחטיכלמנסעפצקרשת")


def sigla(group: frozenset[str]) -> str:
    return ", ".join(sorted(group))


def build_stemma(app: ListApparatus, selected_keys, full_witnesses, ref: str,
                 root: str | None = None):
    """Construit un UnorientedStemma à partir des seuls groupes sélectionnés.

    L'univers complet des témoins est réinjecté pour garder un groupe-racine
    correct, même si un témoin n'apparaît que dans des groupes décochés.

    Args:
        root: témoin placé à la racine de l'arbre. Par défaut, la référence de
            l'apparat (dans la méthode par les groupes, la racine *est* le
            témoin de référence contre lequel les fautes sont relevées).
    """
    sel_groups = [g for g in app.groups
                  if g.is_error and g.witnesses in selected_keys]
    sub = ListApparatus(groups=sel_groups, ref_witness=ref)
    sub.witnesses = set(full_witnesses)
    stemma = UnorientedStemma(sub)
    stemma.compute_edges()
    stemma.build_tree(root_witness=root)
    return stemma


def build_stemma_from_faults(app: ListApparatus, faulty_sides, full_witnesses,
                             ref: str, root: str | None = None):
    """Construit un UnorientedStemma à partir des seuls *côtés fautifs* choisis.

    Chaque élément de ``faulty_sides`` est un ensemble de témoins représentant
    le côté innovant d'une variation ; il peut contenir la référence (côté
    complémentaire). Le treillis de couverture est reconstruit à partir de ces
    côtés, si bien que la racine peut être un témoin autre que la référence.
    """
    sub = ListApparatus(groups=[g for g in app.groups if g.is_error],
                        ref_witness=ref)
    sub.witnesses = set(full_witnesses)
    stemma = UnorientedStemma(sub)
    stemma.members = set(full_witnesses)
    stemma.build_from_faults(faulty_sides, root_witness=root)
    return stemma


def _new_axes(figsize):
    fig, ax = plt.subplots(figsize=figsize)
    fig.patch.set_facecolor(PALETTE["paper"])
    ax.set_facecolor(PALETTE["paper"])
    ax.axis("off")
    return fig, ax


def draw_stemma(tree: dict) -> plt.Figure:
    edges, children = [], defaultdict(list)
    leaf_x, depth_of, label_of, nodes = {}, {}, {}, []
    counter, nid = [0], [0]

    def rec(label, kids, depth, parent):
        nid[0] += 1
        me = nid[0]
        label_of[me] = label
        depth_of[me] = depth
        if parent is not None:
            children[parent].append(me)
            edges.append((parent, me))
        nodes.append(me)
        if not kids:
            leaf_x[me] = counter[0]
            counter[0] += 1
        else:
            for c, gc in kids.items():
                rec(c, gc, depth + 1, me)
        return me

    root_label = next(iter(tree))
    rec(root_label, tree[root_label], 0, None)

    def x_of(me):
        if me in leaf_x:
            return leaf_x[me]
        return sum(x_of(c) for c in children[me]) / len(children[me])

    pos = {me: (x_of(me), -depth_of[me]) for me in nodes}
    max_depth = max(depth_of.values()) if depth_of else 0
    width = max(counter[0], 1)

    fig, ax = _new_axes((max(5.0, width * 1.35),
                         max(3.2, (max_depth + 1) * 1.25)))

    for p, c in edges:
        ax.plot([pos[p][0], pos[c][0]], [pos[p][1], pos[c][1]],
                color=PALETTE["indigo"], lw=1.4, zorder=1)

    for me in nodes:
        x, y = pos[me]
        lab = label_of[me]
        if lab in HYP_CHARS:
            ax.add_patch(Circle((x, y), 0.18, facecolor=PALETTE["paper"],
                                edgecolor=PALETTE["muted"], lw=1.4, zorder=2))
            ax.text(x, y, lab, ha="center", va="center", zorder=3,
                    fontsize=12, color=PALETTE["muted"], style="italic")
        else:
            w = 0.16 + 0.10 * len(lab)
            ax.add_patch(FancyBboxPatch((x - w / 2, y - 0.16), w, 0.32,
                         boxstyle="round,pad=0.02,rounding_size=0.08",
                         facecolor="white", edgecolor=PALETTE["indigo"],
                         lw=1.4, zorder=2))
            ax.text(x, y, lab, ha="center", va="center", zorder=3,
                    fontsize=12, color=PALETTE["ink"],
                    fontweight="500", family="monospace")

    rx, ry = pos[1]
    ax.text(rx, ry + 0.42, "racine / référence", ha="center", va="bottom",
            fontsize=8.5, color=PALETTE["muted"], style="italic")

    ax.set_xlim(-1, width)
    ax.set_ylim(-max_depth - 0.8, 0.9)
    fig.tight_layout()
    return fig


def draw_lattice(stemma: UnorientedStemma) -> plt.Figure:
    children_of = stemma.children_of or {}
    overlaps = stemma.contamination_overlap or []

    G = nx.DiGraph()
    for parent, kids in children_of.items():
        G.add_node(parent)
        for k in kids:
            G.add_node(k)
            G.add_edge(parent, k)
    for g1, g2 in overlaps:
        G.add_node(g1)
        G.add_node(g2)

    if not G.nodes:
        fig, ax = _new_axes((5, 2))
        ax.text(0.5, 0.5, "Aucun groupe sélectionné", ha="center",
                va="center", color=PALETTE["muted"])
        return fig

    universe_size = max(len(n) for n in G.nodes)
    layers = defaultdict(list)
    for n in G.nodes:
        layers[len(n)].append(n)
    pos = {}
    for size, ns in layers.items():
        ns = sorted(ns, key=lambda s: sorted(s))
        for i, n in enumerate(ns):
            x = i - (len(ns) - 1) / 2
            pos[n] = (x * 1.6, size)

    conflict_nodes = {n for pair in overlaps for n in pair}
    span = max((p[0] for p in pos.values()), default=1)
    fig, ax = _new_axes((max(5.5, span * 1.4 + 3),
                         max(3.2, universe_size * 0.95)))

    for p, c in G.edges:
        ax.annotate("", xy=pos[c], xytext=pos[p],
                    arrowprops=dict(arrowstyle="-|>", color=PALETTE["indigo"],
                                    lw=1.3, shrinkA=14, shrinkB=14, alpha=0.85),
                    zorder=1)

    for g1, g2 in overlaps:
        shared = sigla(g1 & g2)
        arc = FancyArrowPatch(pos[g1], pos[g2], connectionstyle="arc3,rad=0.28",
                              arrowstyle="-", linestyle=(0, (4, 3)),
                              color=PALETTE["rubric"], lw=1.8, zorder=2)
        ax.add_patch(arc)
        mx = (pos[g1][0] + pos[g2][0]) / 2
        my = (pos[g1][1] + pos[g2][1]) / 2 + 0.18
        ax.text(mx, my, shared, ha="center", va="center", fontsize=8.5,
                color=PALETTE["rubric"], style="italic", zorder=4,
                bbox=dict(boxstyle="round,pad=0.15", fc=PALETTE["paper"],
                          ec="none"))

    for n in G.nodes:
        x, y = pos[n]
        is_universe = len(n) == universe_size
        in_conflict = n in conflict_nodes
        label = "Ω" if is_universe else sigla(n)
        edge = PALETTE["rubric"] if in_conflict else (
            PALETTE["muted"] if is_universe else PALETTE["indigo"])
        face = "#fbeeea" if in_conflict else "white"
        w = 0.22 + 0.16 * len(label)
        ax.add_patch(FancyBboxPatch((x - w / 2, y - 0.2), w, 0.4,
                     boxstyle="round,pad=0.02,rounding_size=0.1",
                     facecolor=face, edgecolor=edge,
                     lw=1.8 if in_conflict else 1.4, zorder=3))
        ax.text(x, y, label, ha="center", va="center", fontsize=11,
                color=PALETTE["ink"], family="monospace", zorder=4)

    ax.text(-span - 1.2, universe_size, "plus inclusif", ha="left",
            va="center", fontsize=8, color=PALETTE["muted"], rotation=90)
    ax.set_xlim(-span - 1.6, span + 1.6)
    ax.set_ylim(0.3, universe_size + 0.8)
    fig.tight_layout()
    return fig
