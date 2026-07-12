"""
Compute different scoring metrics for comparing two different stemma trees.
"""

import re
from typing import Optional, Union
import dendropy
import networkx as nx
from dendropy.calculate import treecompare


# LIST OF POSSIBLE PLACE HOLDERS FOR LOST HYPERARCHETYPES
_GREEK_CHARS = set("αβγδεζηθικλμνξοπρστυφχψωΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩ")
_GREEK_NAMES = ("alpha beta gamma delta epsilon zeta eta theta iota kappa "
                "lambda mu nu xi omicron pi rho sigma tau upsilon phi chi "
                "psi omega").split()
_PLACEHOLDER_RE = re.compile(
    r"^(inner|node|clade|htu|" + "|".join(_GREEK_NAMES) + r"|root|archetype)\d*$",
    re.IGNORECASE,
)


def _is_placeholder(label):
    """Find out if a node has a placeholder label
    or not, i.e. can be considered as a lost hyperarchetype.

    Args:
        label (str): The label given to the node.

    Returns:
        bool: Whether or not the node can be considered
            as an hyperarchetype.
    """
    if label is None:
        return True
    s = str(label).strip()
    if s == "":
        return True
    if _PLACEHOLDER_RE.match(s):
        return True
    if all(ch in _GREEK_CHARS for ch in s):
        return True
    return False


def _all_labels(*trees):
    """Obtain all the node labels in a list
    of trees.

    Args:
        trees (dendropy.Tree): A list of trees.

    Output:
        set[str]: The list of available labels.
    """
    out = set()
    for t in trees:
        for nd in t:
            if nd.taxon is not None:
                out.add(nd.taxon.label)
            if nd.label is not None:
                out.add(nd.label)
    return out


def infer_witnesses(*trees):
    """Get all the witnesses present in a given tree,
    i.e. the difference between the found labels and the
    labels detected as indicating hyperarchetypes.

    Args:
        trees (dendropy.Tree): A list of trees.

    Returns:
        set[str]: A set of labels.
    """
    return {lab for lab in _all_labels(*trees) if not _is_placeholder(lab)}


def _normalize(tree, witnesses):
    """Normalize a tree for proper computation of metrics
    that measure clade quality (i.e. Robinsons-Foulds).
    Hyperarchetypes are transformed to be considered as leaves, so that
    the bipartitions can be fairly computed between the group method
    and the phylogenetics methods.
    """
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
        if norm:
            g = nd.new_child()
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


# --------------------------------------------------------------------------
# Clade based distances
# --------------------------------------------------------------------------


def clade_based_metrics(true_newick,
                        recon_newick,
                        witnesses=None,
                        rooted=False,
                        restrict_to_common=False):
    """Compute the metrics able to capture if the trees have similar
    clades:
        - Robinson's Foulds
        - Robinson's Foulds normalized by ?
        - 1 - Normalized RF (to obtain a similarity)
        - Precision in detecting clades
        - Recall in detecting clades
        - F1 in detecting clades
        - Number of true clades recovered
        - Wrongly detected caldes
        - Number of witnesses compared
        - Witnesses missing from recon
        - Witnesses extra in recon
        - If the trees should were considered as rooted
            when computing the metrics.

    Args:
        true_newick (str): The ground truth tree.
        recon_newick (str): The reconstructed tree to compare.
        witnesses (set[str]): The list of witnesses.
            If set to None, witnesses are inferred automatically.
        rooted (bool): If the tree should be considered as rooted.
            Defaults to False.
        restrict_to_common (bool): Whether or not to restrict the
            computation of the clades only to common leaves.
    """
    # Load the two trees
    t_true = dendropy.Tree.get(data=true_newick, schema="newick",
                               suppress_internal_node_taxa=False)
    t_rec = dendropy.Tree.get(data=recon_newick, schema="newick",
                              suppress_internal_node_taxa=False)

    if witnesses is None:
        witnesses = infer_witnesses(t_true, t_rec)

    # Normalize the trees to compute metric using witnesses
        # positioned as archetypes.
    _normalize(t_true, witnesses)
    _normalize(t_rec, witnesses)

    # Witnesses in the true tree
    leaves_true = {l.taxon.label for l in t_true.leaf_node_iter() if l.taxon}
    # Witnesses in the reconstructed tree
    leaves_rec = {l.taxon.label for l in t_rec.leaf_node_iter() if l.taxon}

    # Differences in available witnesses between the two trees
    missing_in_rec = leaves_true - leaves_rec
    extra_in_rec = leaves_rec - leaves_true

    # Filter out witnesses if only common witnesses should be retained
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

    # Use tree compare to detect clade detection fp and fn
    fp, fn = treecompare.false_positives_and_negatives(t_true, t_rec)
    # Reference bipartition to compute tp
    R = sum(1 for b in t_true.bipartition_encoding if not b.is_trivial())
    tp = R - fn
    # Compute Robinson fould
    rf = fp + fn
    # Compute precision, recall, n1
    precision = tp / (tp + fp) if (tp + fp) else 1.0
    recall = tp / (tp + fn) if (tp + fn) else 1.0
    f1 = 2 * precision * recall / \
        (precision + recall) if (precision + recall) else 0.0
    rf_max = 2 * tp + fp + fn
    # Normalize RF by the maximum
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
#   skeleton (no orientation)  : for every witness triplets (i,j,k), compare
#                 sign(dist(i,j) vs dist(i,k)) in true vs reconstructed graph.
#   orientation (takes into account orientation) : ordered witness pairs (i,j),
#                   compare sign(depth_to_root(i) vs depth_to_root(j)).
# Award 1 / 0.5 / 0 for full / partial / no agreement;
# --------------------------------------------------------------------------

def _roos_sign(x: int, y: int) -> int:
    """Measure signs difference for computing the metric.

    Args:
        x (int): The sign of the comparison of the minimum paths.
        y (int): The sign of the comparison of the minimum paths.

    Returns:
        int: The sign difference.
    """
    if x == -1 and y == -1:
        return 0
    if x == -1:
        return -1
    if y == -1:
        return 1
    if x < y:
        return -1
    if x > y:
        return 1
    return 0


def _build_graph(newick: str,
                 witnesses: set[str]):
    """Transform a newick string into a networkx graph to compute
    the shortest paths between triplet of witnesses.

    Args:
        newick (str): The tree to convert into a graph as a newick string.
        witnesses (set[str]): The set of witnesses to use.

    Returns:
        graph, dict[str, str], int: The newick tree transformed into a graph,
            the mapping between the labels in the networkx graph and the
            original witness label, the node corresponding to the rooting one.
    """
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
    return g, wit


def witness_dist(dd, w, a, b):
    """Returns the distance between two witnesses given
    a dictionary containing the shortest paths.

    Args:
        dd (dict[int, dict[int, int]]): The shortest path lengths
            between nodes in the graph.
        w (dict[str, int]): The mapping between witness labels and
            node ids.
        a (str): The label of the first witness.
        b (str): The label of the second witness.

    Returns:
        int: The distance between the two witnesses.
    """
    return dd[w[a]].get(w[b], -1)


def _graph_from_newick(newick: str, witnesses: set[str]):
    """Transform a newick string into a networkx graph (class Graph).

    Args:
        newick (str): The newick string to convert.
        witnesses (set[str]): The set of witnesses to use.

    Returns:
        tuple: Graph and witnesses set.
    """
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
    return g, wit


def _graph_from_nx(graph: nx.Graph, witnesses: set[str]):
    """Add the proper label (i.e. witness name) given an input graph.

    Args:
        graph (nx.Grap): The input graph to annotate.
        witnesses (set[str]): The set of witnesses to use.

    Returns:
        tuple: Graph and witnesses set.
    """
    g = graph.to_undirected() if graph.is_directed() else graph
    wit = {n: n for n in g.nodes if n in witnesses}
    return g, wit


def _prepare(source: Union[str, nx.Graph], witnesses: set[str]):
    """Prepare an input for average distance computation,
    depending on if they are a newick string or a networkx graph.

    Args:
        source: The input source, either as a newick string or an already
            constructed network graph.
        witnesses (set[str]): The set of witnesses to use.

    Returns:
        tuple: Graph and witnesses set prepared for computation of distances.
    """
    if isinstance(source, str):
        return _graph_from_newick(source, witnesses)
    if isinstance(source, nx.Graph):
        return _graph_from_nx(source, witnesses)
    raise TypeError(f"Unsupported input format {type(source)!r}")


def average_sign_distance(true_source: Union[str, nx.Graph],
                          recon_source: Union[str, nx.Graph],
                          witnesses: set[str]) -> dict[str, Union[int, float]]:
    """Compute the average sign distance between the reference tree
    represented either as a newick string or a networkx Graph 
    and the reconstructed tree as either a newick string or a networkx Graph.

    Args:
        true_source (Union[str, nx.Graph]): The reference tree
        recon_source (Union[str, nx.Graph]): The reconstructed tree
        witnesses (set[str]): The list of witnesses as a set.

    Returns:
        dict[str, Union[int, float]]: The finale score as a percentage and
        the number of analyzed witnesses.
    """
    # Transform each tree in a graph for computing the shortest
    # path
    gt, wt = _prepare(true_source, witnesses)
    gp, wp = _prepare(recon_source, witnesses)
    common = sorted(set(wt) & set(wp))
    dt = dict(nx.all_pairs_shortest_path_length(gt))
    dp = dict(nx.all_pairs_shortest_path_length(gp))

    s = tot = 0.0
    # Iterate over every triplet to compute the sign
    # Then sum the sign over the whole tree
    for i in common:
        for j in common:
            for k in common:
                if i != j and i != k and j != k:
                    st = _roos_sign(witness_dist(dt, wt, i, j),
                                    witness_dist(dt, wt, i, k))
                    sp = _roos_sign(witness_dist(dp, wp, i, j),
                                    witness_dist(dp, wp, i, k))
                    s += 1 - abs(st - sp) / 2.0
                    tot += 1
    out = {"skeleton_sign_distance": round(100.0 * s / tot, 3)
           if tot else None,
           "n_witnesses_scored": len(common)}
    return out

# -------
# METRICS RELATIVE TO INTERNAL ANCESTOR IDENTIFICATION
# ------


def internal_ancestor_identification(true_newick: str,
                                     recon_newick: str,
                                     witnesses: set[str]):
    """This metric measures if the 'internal' witnesses of a tree
    are properly reconstructed.

    A witness is "internal" if its node has degree >= 2 in the
    undirected graph, i.e it has descendants.
    A witness is "terminal" if degree 1 (a leaf).

    We define:
      * FN: internal in truth but terminal in recon
        -> a superfluous hyparchetype was introduced.
      * FP: terminal in truth but internal in recon
        -> a witness was wrongly promoted to ancestor.

    Args:
        true_newick (str): The reference tree in Newick format.
        recon_newick (str): The reconstructed tree in Newick format.
        witnesses (set[str]): The list of witnesses as a set.

    Returns:
        dict[str, Any]: the metric and their value.
    """
    # Build the graphs
    gt, wt = _build_graph(true_newick, witnesses)
    gp, wp = _build_graph(recon_newick, witnesses)
    # Restrict to common witnesses
    common = set(wt) & set(wp)

    # Count the number of degree per witness
    true_internal = {s for s in common if gt.degree(wt[s]) >= 2}
    recon_internal = {s for s in common if gp.degree(wp[s]) >= 2}

    # Compute TP, FN, FP (see definition above)
    tp = true_internal & recon_internal
    fn = true_internal - recon_internal
    fp = recon_internal - true_internal

    precision = len(tp) / (len(tp) + len(fp)) if (tp or fp) else 1.0
    recall = len(tp) / (len(tp) + len(fn)) if (tp or fn) else 1.0
    f1 = 2 * precision * recall / \
        (precision + recall) if (precision + recall) else 0.0

    return {
        "internal_ancestor_recall": round(recall, 4),
        "internal_ancestor_precision": round(precision, 4),
        "internal_ancestor_f1": round(f1, 4),
        "correctly_identified": f"{len(tp)}/{len(true_internal)}",
        "missed_ancestors": sorted(fn),
        "spurious_ancestors": sorted(fp),
    }

# --- EVALUATE USING ALL METRICS


def _graph_to_newick(g: nx.DiGraph) -> str:
    """Graphe arborescent -> Newick, labels internes préservés."""
    roots = [n for n in g.nodes if g.in_degree(n) == 0]
    if len(roots) != 1:
        raise ValueError(f"racine non unique: {roots}")

    def rec(n):
        ch = sorted(g.successors(n))
        return f"({','.join(rec(c) for c in ch)}){n}" if ch else str(n)

    return rec(roots[0]) + ";"


def _as_newick(src: Union[str, nx.Graph]) -> Optional[str]:
    """Transform an input source (either string or nx.Graphs) into a
    newick file.
    If the file is already in newick format, it is returned
    without transformation.
    Note that if the tree is contaminated, the transformation fails
    as newick does not accept duplicated taxa.

    Args:
        src (Union[str, nx.Graph]): The input source to convert.

    Returns:
        Optional[str]: The converted source. If the input graph cannot
            be transformed into a tree, it returns None.
    """
    if isinstance(src, str):
        return src
    if isinstance(src, nx.Graph):
        if src.is_directed() and nx.is_arborescence(src):
            return _graph_to_newick(src)
        return None
    raise TypeError(f"source non supportée: {type(src)!r}")


def _witnesses_of(src: Union[str, nx.Graph]) -> set[str]:
    """Given a tree represented as a newick or a graph,
    obtain the list of witnesses.
    
    Args:
        src (Union[str, nx.Graph]): The tree or the graph to
            return.
    
    Returns:
        set[str]: The list of witnesses.
    """
    if isinstance(src, nx.Graph):
        return set(src.nodes)
    t = dendropy.Tree.get(data=src, schema="newick",
                          suppress_internal_node_taxa=False)
    return {lab for nd in t.preorder_node_iter()
            if (lab := (nd.taxon.label if nd.taxon is not None else nd.label))}


def evaluate(true_source: Union[str, nx.Graph],
             recon_source: Union[str, nx.Graph],
             witnesses: Optional[set[str]] = None,
             rooted: bool = False):
    """Compute all available metrics and returns them as a
    dictionary. Note that if the tree has reticulation (i.e.
    contamination), the RF distance is set to None.

    Args:
        true_source (Union[str, nx.Graph]):
        recon_source (Union[str, nx.Graph]):
        witnesses (set[str]): The list of witnesses as a set.
            If set to None, it is automatically computed from the
            input witnesses.
        rooted (bool): Whether or not the RF distance should be
            computed on a rooted tree..
        
    Returns:
        dict[str, ?]: The computed scores.
    """
    tn, rn = _as_newick(true_source), _as_newick(recon_source)

    if witnesses is None:
        witnesses = _witnesses_of(true_source) & _witnesses_of(recon_source)

    roos_block = average_sign_distance(true_source, recon_source, witnesses)

    if tn is None or rn is None:
        rf_block = anc_block = None
    else:
        rf_block = clade_based_metrics(tn, rn, witnesses=witnesses,
                                       rooted=rooted)
        anc_block = internal_ancestor_identification(tn, rn, witnesses)

    return {"robinson_foulds": rf_block,
            "roos_sign_distance": roos_block,
            "internal_ancestors": anc_block,
            "tree_representable": tn is not None and rn is not None}
