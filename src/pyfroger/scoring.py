"""
stemma_compare.py — score a reconstructed stemma against the true one.

Inner-node handling
-------------------
Two kinds of internal labels, treated differently:

1. Hypothetical ancestors with no cross-tree identity:
     - your Greek letters: "α"/"alpha", "β"/"beta", ...
     - Bio.Phylo's "Inner1", "Inner2", ...
   "alpha" in the truth has nothing to do with "Inner3" in the reconstruction.
   These are anonymised (label -> None) so they never count as taxa.

2. Extant witnesses sitting at an internal position (your p10 = parent of p5,p11):
   real identity, position matters. We re-express them so the metric can see
   them: node X labelled p10 with children {p5,p11} becomes
       X( pendant=p10 , g(p5,p11) )
   i.e. p10 is marked at that level while {p5,p11} stays a clade. This is what
   lets the score separate "p10 is the ancestor of p5,p11" (true) from
   "p10 is a leaf beside p5,p11" (a wrong reconstruction).

Discriminator: membership in the witness universe (real sigla). Inferred by
default as every label that is not a placeholder; override with `witnesses=`.
"""

import re
import dendropy
import networkx as nx
from dendropy.calculate import treecompare
import itertools
import math


_GREEK_CHARS = set("αβγδεζηθικλμνξοπρστυφχψωΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩ")
_GREEK_NAMES = ("alpha beta gamma delta epsilon zeta eta theta iota kappa "
                "lambda mu nu xi omicron pi rho sigma tau upsilon phi chi "
                "psi omega").split()
_PLACEHOLDER_RE = re.compile(
    r"^(inner|node|clade|htu|" + "|".join(_GREEK_NAMES) + r"|root|archetype)\d*$",
    re.IGNORECASE,
)


def _is_placeholder(label):
    if label is None:
        return True
    s = str(label).strip()
    if s == "":
        return True
    if _PLACEHOLDER_RE.match(s):
        return True
    if all(ch in _GREEK_CHARS for ch in s):      # pure Greek unicode
        return True
    return False


def _all_labels(*trees):
    out = set()
    for t in trees:
        for nd in t:
            if nd.taxon is not None:
                out.add(nd.taxon.label)
            if nd.label is not None:
                out.add(nd.label)
    return out


def infer_witnesses(*trees):
    return {lab for lab in _all_labels(*trees) if not _is_placeholder(lab)}


def _normalise(tree, witnesses):
    """Anonymise placeholders; re-express internal witnesses as marked leaves.
    Uses the tree's OWN taxon namespace (no cross-tree mixing here)."""
    tns = tree.taxon_namespace
    for nd in tree:
        lab = nd.taxon.label if nd.taxon is not None else nd.label
        nd._norm = lab if (lab is not None and lab in witnesses) else None

    for nd in list(tree.leaf_node_iter()):
        nd.taxon = tns.require_taxon(label=nd._norm) if nd._norm else None
        nd.label = None

    for nd in list(tree.preorder_internal_node_iter()):
        norm = getattr(nd, "_norm", None)
        nd.label = None
        nd.taxon = None
        if norm:                                  # extant witness as ancestor
            g = nd.new_child()                    # grouping node holds old subtree
            for ch in list(nd.child_nodes()):
                if ch is g:
                    continue
                nd.remove_child(ch)
                g.add_child(ch)
            pendant = nd.new_child(edge_length=0.0)
            pendant.taxon = tns.require_taxon(label=norm)
            nd._norm = None
    tree.purge_taxon_namespace()
    return tree

def _grafted_witness_distances(newick, witnesses):
    """Hop-distances between witnesses AFTER grafting every witness-ancestor to a
    pendant leaf. Grafting is uniform, so an internal witness and a leaf witness
    look identical to the metric: only the *grouping* is scored, never whether a
    witness was realised as an ancestor. Ancestor-agnostic by construction."""
    t = dendropy.Tree.get(data=newick, schema="newick",
                          suppress_internal_node_taxa=False)
    _normalise(t, witnesses)                       # your existing function
    g = nx.Graph(); ids = {}; wit = {}
    for k, nd in enumerate(t.preorder_node_iter()):
        ids[nd] = k; g.add_node(k)
        lab = nd.taxon.label if nd.taxon is not None else nd.label
        if lab in witnesses:
            wit[lab] = k
    for nd in t.preorder_node_iter():
        if nd.parent_node is not None:
            g.add_edge(ids[nd.parent_node], ids[nd])
    d = dict(nx.all_pairs_shortest_path_length(g))
    common = sorted(wit)
    return common, {(a, b): d[wit[a]].get(wit[b]) for a in common for b in common}


def _quartet_topology(D, a, b, c, d):
    s = {("ab", "cd"): D[(a, b)] + D[(c, d)],
         ("ac", "bd"): D[(a, c)] + D[(b, d)],
         ("ad", "bc"): D[(a, d)] + D[(b, c)]}
    m = min(s.values())
    win = [k for k, v in s.items() if v == m]
    return frozenset(win[0]) if len(win) == 1 else None   # None = unresolved/star


def quartet_agreement(true_newick, recon_newick, witnesses):
    """Estabrook-style quartet congruence on the grafted witness graph.
    Unrooted, orientation-free, ancestor-blind."""
    import itertools
    cw, Dt = _grafted_witness_distances(true_newick, witnesses)
    cr, Dr = _grafted_witness_distances(recon_newick, witnesses)
    common = sorted(set(cw) & set(cr))
    agree = res_agree = resolved_true = total = 0
    for q in itertools.combinations(common, 4):
        qt, qr = _quartet_topology(Dt, *q), _quartet_topology(Dr, *q)
        total += 1
        if qt == qr: agree += 1
        if qt is not None:
            resolved_true += 1
            if qt == qr: res_agree += 1
    return {
        "quartet_agreement": round(agree / total, 4) if total else 1.0,
        "quartet_distance": round(1 - agree / total, 4) if total else 0.0,
        "resolved_quartet_recall": round(res_agree / resolved_true, 4) if resolved_true else None,
        "n_quartets": total,
    }


def path_distance_agreement(true_newick, recon_newick, witnesses):
    """Spearman correlation + normalised RMSE of pairwise witness hop-distances.
    Continuous companion to the quartet score; same ancestor-blind grafting."""
    cw, Dt = _grafted_witness_distances(true_newick, witnesses)
    cr, Dr = _grafted_witness_distances(recon_newick, witnesses)
    common = sorted(set(cw) & set(cr))
    pairs = list(itertools.combinations(common, 2))
    xs, ys = [Dt[p] for p in pairs], [Dr[p] for p in pairs]
    n = len(pairs)

    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0] * len(v); i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]: j += 1
            for k in range(i, j + 1): r[order[k]] = (i + j) / 2 + 1
            i = j + 1
        return r

    rx, ry = rank(xs), rank(ys)
    mx, my = sum(rx) / n, sum(ry) / n
    cov = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    sx = math.sqrt(sum((rx[i] - mx) ** 2 for i in range(n)))
    sy = math.sqrt(sum((ry[i] - my) ** 2 for i in range(n)))
    rmse = math.sqrt(sum((xs[i] - ys[i]) ** 2 for i in range(n)) / n)
    denom = math.sqrt(sum(x * x for x in xs) / n) or 1.0
    return {"path_spearman": round(cov / (sx * sy), 4) if sx and sy else 0.0,
            "path_nrmse": round(rmse / denom, 4), "n_pairs": n}

def compare_stemmata(true_newick, recon_newick, *,
                     witnesses=None, rooted=False, restrict_to_common=False):
    t_true = dendropy.Tree.get(data=true_newick, schema="newick",
                               suppress_internal_node_taxa=False)
    t_rec = dendropy.Tree.get(data=recon_newick, schema="newick",
                              suppress_internal_node_taxa=False)

    if witnesses is None:
        witnesses = infer_witnesses(t_true, t_rec)

    _normalise(t_true, witnesses)
    _normalise(t_rec, witnesses)

    leaves_true = {l.taxon.label for l in t_true.leaf_node_iter() if l.taxon}
    leaves_rec = {l.taxon.label for l in t_rec.leaf_node_iter() if l.taxon}
    missing_in_rec = leaves_true - leaves_rec
    extra_in_rec = leaves_rec - leaves_true

    if restrict_to_common and (missing_in_rec or extra_in_rec):
        common = leaves_true & leaves_rec
        t_true.retain_taxa_with_labels(common)
        t_rec.retain_taxa_with_labels(common)

    # one shared namespace for the actual comparison
    shared = dendropy.TaxonNamespace()
    t_true.migrate_taxon_namespace(shared)
    t_rec.migrate_taxon_namespace(shared)
    t_true.is_rooted = t_rec.is_rooted = rooted
    t_true.encode_bipartitions()
    t_rec.encode_bipartitions()

    fp, fn = treecompare.false_positives_and_negatives(t_true, t_rec)
    R = sum(1 for b in t_true.bipartition_encoding if not b.is_trivial())  # ref clades
    rf = fp + fn
    tp = R - fn
    precision = tp / (tp + fp) if (tp + fp) else 1.0
    recall = tp / (tp + fn) if (tp + fn) else 1.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    rf_max = 2 * tp + fp + fn
    norm_rf = rf / rf_max if rf_max else 0.0

    return {
        "rf_distance": rf,
        "normalised_rf": round(norm_rf, 4),
        "similarity": round(1 - norm_rf, 4),
        "clade_precision": round(precision, 4),
        "clade_recall": round(recall, 4),
        "clade_f1": round(f1, 4),
        "true_clades_recovered": f"{tp}/{R}",
        "spurious_clades": fp,
        "n_taxa_compared": len(leaves_true & leaves_rec),
        "witnesses_missing_from_recon": sorted(missing_in_rec),
        "witnesses_extra_in_recon": sorted(extra_in_rec),
        "rooted": rooted,
    }


# --------------------------------------------------------------------------
# Roos & Heikkila (2009) average sign distance  (DSH 24(4):417-433)
# Faithful port of the reference scorer tripled.py, adapted to Newick input.
#   skeleton    : ordered witness triplets (i,j,k); compare
#                 sign(dist(i,j) vs dist(i,k)) in true vs reconstructed graph.
#   orientation : ordered witness pairs (i,j); compare
#                 sign(depth_to_root(i) vs depth_to_root(j)).
# Award 1 / 0.5 / 0 for full / partial / no agreement; report as a percentage.
# Distances and depths traverse hidden nodes (Greek / Inner placeholders);
# the loops range only over real witnesses (the original isalpha() filter).
# Calibration from the paper: ~50% = chance, ~60% = average, ~75% = best in
# practice, 100% = perfect.
# --------------------------------------------------------------------------

def _roos_sign(x, y):                 # hop counts; -1 = unreachable
    if x == -1 and y == -1: return 0
    if x == -1: return -1
    if y == -1: return 1
    if x < y: return -1
    if x > y: return 1
    return 0


def _build_graph(newick, witnesses):
    """Undirected graph over all nodes, {siglum: node_id}, and the root id."""
    t = dendropy.Tree.get(data=newick, schema="newick",
                          suppress_internal_node_taxa=False)
    g = nx.Graph()
    wit, ids = {}, {}
    for k, nd in enumerate(t.preorder_node_iter()):
        ids[nd] = k
        g.add_node(k)
        lab = nd.taxon.label if nd.taxon is not None else nd.label
        if lab is not None and lab in witnesses:
            wit[lab] = k
    for nd in t.preorder_node_iter():
        if nd.parent_node is not None:
            g.add_edge(ids[nd.parent_node], ids[nd])
    return g, wit, ids[t.seed_node]


def average_sign_distance(true_newick, recon_newick, witnesses,
                          do_orientation=True):
    gt, wt, rt = _build_graph(true_newick, witnesses)
    gp, wp, rp = _build_graph(recon_newick, witnesses)
    common = sorted(set(wt) & set(wp))
    dt = dict(nx.all_pairs_shortest_path_length(gt))
    dp = dict(nx.all_pairs_shortest_path_length(gp))

    def d(dd, w, a, b):
        return dd[w[a]].get(w[b], -1)

    s = tot = 0.0
    for i in common:
        for j in common:
            for k in common:
                if i != j and i != k and j != k:
                    st = _roos_sign(d(dt, wt, i, j), d(dt, wt, i, k))
                    sp = _roos_sign(d(dp, wp, i, j), d(dp, wp, i, k))
                    s += 1 - abs(st - sp) / 2.0
                    tot += 1
    out = {"skeleton_sign_distance_pct": round(100.0 * s / tot, 3) if tot else None,
           "n_witnesses_scored": len(common)}

    if do_orientation:
        deptht = nx.single_source_shortest_path_length(gt, rt)
        depthp = nx.single_source_shortest_path_length(gp, rp)
        s2 = tot2 = 0.0
        for i in common:
            for j in common:
                if i != j:
                    st = _roos_sign(deptht.get(wt[i], -1), deptht.get(wt[j], -1))
                    sp = _roos_sign(depthp.get(wp[i], -1), depthp.get(wp[j], -1))
                    s2 += 1 - abs(st - sp) / 2.0
                    tot2 += 1
        out["orientation_sign_distance_pct"] = (
            round(100.0 * s2 / tot2, 3) if tot2 else None)
    return out


def internal_ancestor_identification(true_newick, recon_newick, witnesses):
    """Did the reconstruction get the *internal-ancestor* witnesses right?

    A witness is "internal" if its node has degree >= 2 in the undirected graph,
    i.e. it sits between other witnesses / has descendants (e.g. A in A->B, or
    p10 = parent of p5,p11). A witness is "terminal" if degree 1 (a leaf).
    This is rooting-free, so it works on UnorientedStemma output too, and it
    isolates the eliminatio / false-common-ancestor phenomenon:
      * miss  (FN): internal in truth but terminal in recon  -> A->B turned into
                    alpha->(A,B); a superfluous hyparchetype was introduced.
      * spurious (FP): terminal in truth but internal in recon -> a witness was
                    wrongly promoted to ancestor.
    Restricted to witnesses present in BOTH trees (missing/extra witnesses are
    reported by the Robinson-Foulds block).
    """
    gt, wt, _ = _build_graph(true_newick, witnesses)
    gp, wp, _ = _build_graph(recon_newick, witnesses)
    common = set(wt) & set(wp)
    true_internal = {s for s in common if gt.degree(wt[s]) >= 2}
    recon_internal = {s for s in common if gp.degree(wp[s]) >= 2}
    tp = true_internal & recon_internal
    fn = true_internal - recon_internal     # missed ancestors (false common anc.)
    fp = recon_internal - true_internal     # spuriously promoted
    precision = len(tp) / (len(tp) + len(fp)) if (tp or fp) else 1.0
    recall = len(tp) / (len(tp) + len(fn)) if (tp or fn) else 1.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    return {
        "internal_ancestor_recall": round(recall, 4),
        "internal_ancestor_precision": round(precision, 4),
        "internal_ancestor_f1": round(f1, 4),
        "correctly_identified": f"{len(tp)}/{len(true_internal)}",
        "missed_ancestors": sorted(fn),
        "spurious_ancestors": sorted(fp),
    }


def evaluate(true_newick, recon_newick, *, witnesses=None,
             rooted=False, orientation=False):
    """One call -> both metric families.

    rooted      : rooted (clade) RF if True, else unrooted (topology only).
    orientation : also compute Roos orientation score. Only meaningful when the
                  reconstruction is genuinely rooted/oriented (NOT for raw
                  UnorientedStemma output); leave False otherwise.
    """
    t_true = dendropy.Tree.get(data=true_newick, schema="newick",
                               suppress_internal_node_taxa=False)
    t_rec = dendropy.Tree.get(data=recon_newick, schema="newick",
                              suppress_internal_node_taxa=False)
    if witnesses is None:
        witnesses = infer_witnesses(t_true, t_rec)

    rf_block = compare_stemmata(true_newick, recon_newick,
                                witnesses=witnesses, rooted=rooted)
    roos_block = average_sign_distance(true_newick, recon_newick, witnesses,
                                       do_orientation=orientation)
    anc_block = internal_ancestor_identification(true_newick, recon_newick,
                                                 witnesses)
    anc_blind = {**quartet_agreement(true_newick, recon_newick, witnesses),
                 **path_distance_agreement(true_newick, recon_newick, witnesses)}
    
    return {"robinson_foulds": rf_block,
            "roos_sign_distance": roos_block,
            "internal_ancestors": anc_block,
            "ancestor_blind": anc_blind}


def _load(path_or_str):
    import os
    if os.path.exists(path_or_str):
        with open(path_or_str, encoding="utf-8") as fh:
            return fh.read()
    return path_or_str


# --------------------------------------------------------------------------
# Diagnostic tool: localise and label every discrepancy between true & recon
# --------------------------------------------------------------------------

def _exemplar_map(newick, witnesses):
    """{witness: nearest witness strictly above it, or None for the archetype}.
    Requires a rooted Newick (uses parent pointers from the seed node)."""
    t = dendropy.Tree.get(data=newick, schema="newick",
                          suppress_internal_node_taxa=False)
    nearest, exo = {}, {}
    for nd in t.preorder_node_iter():
        lab = nd.taxon.label if nd.taxon is not None else nd.label
        above = None if nd.parent_node is None else nearest[nd.parent_node]
        if lab in witnesses:
            exo[lab] = above
            nearest[nd] = lab
        else:
            nearest[nd] = above
    return exo


def _clade_sets(newick, witnesses, rooted):
    """Non-trivial clades as frozensets of witness labels (witnesses-as-ancestors
    grafted in, exactly as the RF score does)."""
    t = dendropy.Tree.get(data=newick, schema="newick",
                          suppress_internal_node_taxa=False)
    _normalise(t, witnesses)
    t.is_rooted = rooted
    t.encode_bipartitions()
    tns = t.taxon_namespace
    out = set()
    for b in t.bipartition_encoding:
        if b.is_trivial():
            continue
        leaves = frozenset(x.label for x in b.leafset_taxa(tns))
        if rooted:
            out.add(leaves)
        else:                                  # canonical unordered split
            comp = frozenset(x.label for x in tns) - leaves
            out.add(frozenset([leaves, comp]))
    return out


def diagnose(true_newick, recon_newick, *, witnesses=None, rooted=True):
    """Structured, per-error comparison of true vs reconstructed stemma."""
    tt = dendropy.Tree.get(data=true_newick, schema="newick",
                           suppress_internal_node_taxa=False)
    tr = dendropy.Tree.get(data=recon_newick, schema="newick",
                           suppress_internal_node_taxa=False)
    if witnesses is None:
        witnesses = infer_witnesses(tt, tr)

    exo_t = _exemplar_map(true_newick, witnesses)
    exo_r = _exemplar_map(recon_newick, witnesses)
    common = set(exo_t) & set(exo_r)
    ARCH = "(archétype)"

    def name(x):
        return ARCH if x is None else x

    # reattachments, categorised
    reattach = []
    for w in sorted(common):
        if exo_t[w] != exo_r[w]:
            t_anc, r_anc = exo_t[w], exo_r[w]
            if t_anc is not None and r_anc is None:
                cat = "ancêtre-témoin perdu (faux ancêtre commun introduit)"
            elif t_anc is None and r_anc is not None:
                cat = "exemplaire inventé (rattaché à un témoin à tort)"
            else:
                cat = "réattaché à un autre témoin"
            reattach.append({"témoin": w, "exemplaire_vrai": name(t_anc),
                             "exemplaire_reconstruit": name(r_anc),
                             "catégorie": cat})

    anc = internal_ancestor_identification(true_newick, recon_newick, witnesses)

    ct, cr = _clade_sets(true_newick, witnesses, rooted), \
             _clade_sets(recon_newick, witnesses, rooted)

    def fmt_clade(c):
        if rooted:
            return "{" + ", ".join(sorted(c)) + "}"
        a, b = list(c)
        return "{" + ", ".join(sorted(a)) + "} | {" + ", ".join(sorted(b)) + "}"

    missing_clades = sorted((fmt_clade(c) for c in (ct - cr)))
    spurious_clades = sorted((fmt_clade(c) for c in (cr - ct)))

    # depth shifts
    gt, wt, rt = _build_graph(true_newick, witnesses)
    gp, wp, rp = _build_graph(recon_newick, witnesses)
    dt = nx.single_source_shortest_path_length(gt, rt)
    dp = nx.single_source_shortest_path_length(gp, rp)
    depth_shifts = []
    for w in sorted(common):
        a, b = dt.get(wt[w]), dp.get(wp[w])
        if a is not None and b is not None and a != b:
            depth_shifts.append({"témoin": w, "profondeur_vraie": a,
                                 "profondeur_reconstruite": b})

    return {
        "witnesses_missing_from_recon": sorted(set(exo_t) - set(exo_r)),
        "witnesses_extra_in_recon": sorted(set(exo_r) - set(exo_t)),
        "reattachments": reattach,
        "demoted_ancestors": anc["missed_ancestors"],
        "spuriously_promoted": anc["spurious_ancestors"],
        "missing_true_clades": missing_clades,
        "spurious_clades": spurious_clades,
        "depth_shifts": depth_shifts,
        "ascii_true": tt.as_ascii_plot(plot_metric="level",
                                       show_internal_node_labels=True),
        "ascii_recon": tr.as_ascii_plot(plot_metric="level",
                                        show_internal_node_labels=True),
    }


def format_diagnosis(d):
    """Human-readable French report from diagnose()."""
    L = []
    P = L.append
    P("=" * 64)
    P("DIAGNOSTIC  vrai  vs  reconstruit")
    P("=" * 64)

    if d["witnesses_missing_from_recon"] or d["witnesses_extra_in_recon"]:
        P("\n## Témoins absents / surnuméraires")
        if d["witnesses_missing_from_recon"]:
            P("  absents du reconstruit : " + ", ".join(d["witnesses_missing_from_recon"]))
        if d["witnesses_extra_in_recon"]:
            P("  surnuméraires : " + ", ".join(d["witnesses_extra_in_recon"]))

    P("\n## Réattachements (exemplaire le plus proche)")
    if not d["reattachments"]:
        P("  aucun — toutes les filiations directes sont correctes.")
    for r in d["reattachments"]:
        P("  %-6s : %s  ->  %s   [%s]" % (
            r["témoin"], r["exemplaire_vrai"], r["exemplaire_reconstruit"],
            r["catégorie"]))

    P("\n## Ancêtres internes")
    P("  rétrogradés en feuille : " +
      (", ".join(d["demoted_ancestors"]) if d["demoted_ancestors"] else "aucun"))
    P("  promus à tort : " +
      (", ".join(d["spuriously_promoted"]) if d["spuriously_promoted"] else "aucun"))

    P("\n## Regroupements (clades)")
    P("  clades vrais manqués : " +
      (" ; ".join(d["missing_true_clades"]) if d["missing_true_clades"] else "aucun"))
    P("  clades inventés : " +
      (" ; ".join(d["spurious_clades"]) if d["spurious_clades"] else "aucun"))

    if d["depth_shifts"]:
        P("\n## Décalages de profondeur (distance à la racine)")
        for s in d["depth_shifts"]:
            P("  %-6s : %d -> %d" % (s["témoin"], s["profondeur_vraie"],
                                     s["profondeur_reconstruite"]))

    P("\n## Arbres (les nœuds internes nommés sont des témoins-ancêtres ;")
    P("   les hyparchétypes apparaissent sous leur étiquette cachée)")
    P("\n-- VRAI --\n" + d["ascii_true"])
    P("-- RECONSTRUIT --\n" + d["ascii_recon"])
    return "\n".join(L)



if __name__ == "__main__":
    import argparse, json
    ap = argparse.ArgumentParser(
        description="Score a reconstructed stemma against the true one "
                    "(Robinson-Foulds / clade P-R-F1 + Roos sign distance).")
    ap.add_argument("true", help="Newick file or string for the TRUE tree")
    ap.add_argument("recon", help="Newick file or string for the RECONSTRUCTED tree")
    ap.add_argument("--rooted", action="store_true", help="rooted RF (default unrooted)")
    ap.add_argument("--orientation", action="store_true",
                    help="also compute Roos orientation score (needs a rooted recon)")
    ap.add_argument("--witnesses", default=None, help="comma-separated sigla override")
    ap.add_argument("--diagnose", action="store_true",
                    help="print a per-error diagnostic report instead of the scores")
    a = ap.parse_args()
    w = set(a.witnesses.split(",")) if a.witnesses else None
    if a.diagnose:
        from_diag = diagnose(_load(a.true), _load(a.recon), witnesses=w,
                             rooted=a.rooted or True)
        print(format_diagnosis(from_diag))
    else:
        print(json.dumps(evaluate(_load(a.true), _load(a.recon), witnesses=w,
                                  rooted=a.rooted, orientation=a.orientation),
                         ensure_ascii=False, indent=2))
