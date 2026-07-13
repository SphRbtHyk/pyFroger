"""Build an unoriented stemma from a list apparatus (Froger's group method).

The *méthode des groupements* (Dom J. Froger, *La critique des textes et son
automatisation*, 1968) groups the witnesses by the errors they share, orders
those groups by inclusion, and reads the genealogy off that order: a group
strictly containing another sits above it.

Two groups that overlap without either containing the other cannot both be
clean branches of one tree. Such an overlap is the sign of a *contamination*
(horizontal transmission): it is either dismissed as polygenesis, or kept, and
it then shows up as a *reticulation*, an edge a Newick tree cannot hold.
Contamination is the philological fact; reticulation is what it does to the
graph.
"""
import random
import re
from collections import Counter, defaultdict

from rich import box
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from pyfroger.list_apparatus import ListApparatus

ALPHABET = "αβγδεζηθικλμνξοπρστυφχψω"
ANOMALY_METHODS = ("raw_cutoff", "two_by_two", "probabilistic_cutoff", "none")

# How many loci must support a group before it is worth trusting. The bands are
# read in order, the first upper bound reached wins; the last one catches the
# rest. Thresholds are relative to the size of the apparatus, so they are meant
# to be recalibrated: pass your own to UnorientedStemma(support_bands=...).
SUPPORT_BANDS = (
    (10, "🔴", "spurious"),
    (20, "🟠", "doubtful"),
    (50, "🟡", "moderate"),
    (float("inf"), "🟢", "strong"),
)


def natural_key(siglum: str) -> list:
    """Sort key reading embedded digits as numbers, so P9 comes before P10."""
    return [int(part) if part.isdigit() else part
            for part in re.split(r"(\d+)", siglum)]


class UnorientedStemma:
    """Reconstruct an unoriented stemma from a variant apparatus.

    Attributes:
        apparatus (ListApparatus): The variant groups and their readings.
        min_frequency (int): Minimum number of loci at which an error group
            must be attested to enter the lattice.
        reference (str): The reference siglum.
        witnesses (frozenset[str]): Every siglum, reference included.
        support (dict[frozenset[str], int]): Locus count of each error group.
        children_of (dict[frozenset[str], list[frozenset[str]]]): Each group
            mapped to the groups it directly covers.
        overlaps (list[tuple[frozenset[str], frozenset[str]]]): Contaminated
            (non-laminar) group pairs surviving anomaly resolution.
        dismissed (set[frozenset[str]]): Groups dropped as polygenesis.
        ties (list[tuple[frozenset[str], frozenset[str]]]): Conflicting pairs
            whose two groups had equal support, so that the errors did not
            decide between them.
        tree (dict[str, dict] | None): Nested-dict stemma once built.
        reticulations (list[tuple[str, str]]): Contamination edges the tree
            could not hold, as ``(parent_label, witness)``.
        node_groups (dict[str, frozenset[str]]): Each node label mapped to the
            error group it stands for.
    """

    def __init__(self,
                 apparatus: ListApparatus,
                 min_frequency: int = 1,
                 verbose: bool = False,
                 seed: int | None = None,
                 support_bands: tuple = SUPPORT_BANDS) -> None:
        """Initialise the builder from an apparatus.

        Args:
            apparatus (ListApparatus): The variant groups from the collation.
            min_frequency (int): Minimum locus count for a group to be kept.
            verbose (bool): Print the contaminated readings during resolution.
            seed (int): Seed of the tie-breaking draw used by ``two_by_two``
                when both groups of a conflict have the same support. Left to
                None, ties are broken deterministically, by siglum order.
            support_bands (tuple): ``(upper_bound, marker, name)`` triples
                grading how well attested a group is. Defaults to
                :data:`SUPPORT_BANDS`; recalibrate them on the size of your
                apparatus.
        """
        self.apparatus = apparatus
        self.min_frequency = min_frequency
        self.verbose = verbose
        self.seed = seed
        self.rng = random.Random(seed)
        self.support_bands = support_bands
        self.console = Console()

        self.reference = apparatus.reference
        self.witnesses = frozenset(apparatus.witnesses | {apparatus.reference})
        self.group_info = apparatus.group_by_frequency(
            min_frequency=min_frequency)
        self.support = {group: info["frequency"]
                        for group, info in self.group_info.items()}

        self.children_of: dict[frozenset[str], list[frozenset[str]]] = {}
        self.overlaps: list[tuple[frozenset[str], frozenset[str]]] = []
        self.dismissed: set[frozenset[str]] = set()
        self.ties: list[tuple[frozenset[str], frozenset[str]]] = []
        self.tree: dict[str, dict] | None = None
        self.reticulations: list[tuple[str, str]] = []
        self.node_groups: dict[str, frozenset[str]] = {}

    # ------------------------------------------------------------------ #
    # Lattice of error groups
    # ------------------------------------------------------------------ #
    def _compute_lattice(self,
                         dismissed: set[frozenset[str]] = frozenset()) -> None:
        """Order the retained groups by inclusion.

        Fills :attr:`children_of` with the covering relations (a group covers
        another when it strictly contains it with nothing in between) and
        :attr:`overlaps` with the non-laminar pairs.

        Args:
            dismissed: Groups to leave out of the lattice.
        """
        groups = {group for group in self.support if group not in dismissed}
        groups.add(self.witnesses)  # the universe closes the lattice on top

        def covers(child, parent):
            return child < parent and not any(
                child < other < parent for other in groups)

        children_of = defaultdict(list)
        overlapping_pairs = set()
        for child in groups:
            for parent in groups:
                if covers(child, parent):
                    children_of[parent].append(child)
                else:
                    shared = child & parent
                    if shared and shared != child and shared != parent:
                        overlapping_pairs.add(frozenset((child, parent)))

        self.children_of = {parent: sorted(children, key=sorted)
                            for parent, children in children_of.items()}
        self.overlaps = [tuple(sorted(pair, key=sorted))
                         for pair in overlapping_pairs]
        self.dismissed = set(dismissed)

    @property
    def edges(self) -> list[tuple[frozenset[str], frozenset[str]]]:
        """The covering relations, as ``(parent_group, child_group)`` pairs."""
        return sorted((parent, child)
                      for parent, children in self.children_of.items()
                      for child in children)

    @property
    def contaminated_witnesses(self) -> frozenset[str]:
        """The witnesses shared by two overlapping groups."""
        return frozenset(witness
                         for first, second in self.overlaps
                         for witness in first & second)

    def contamination_loci(self) -> set[str]:
        """Return the loci attesting the groups involved in a contamination."""
        return {locus
                for pair in self.overlaps
                for group in pair
                for locus in self.group_info[group]["loci"]}

    # ------------------------------------------------------------------ #
    # Anomaly resolution
    # ------------------------------------------------------------------ #
    def _dismiss_in_pair(self,
                         pair: tuple[frozenset[str], frozenset[str]],
                         anomaly_method: str,
                         threshold: float) -> set[frozenset[str]]:
        """Return which side(s) of a contaminated pair to drop as polygenesis.

        Args:
            pair: Two overlapping error groups.
            anomaly_method (str): See :meth:`build_tree`.
            threshold (float): A locus count, or a probability for
                ``probabilistic_cutoff``.

        Returns:
            set[frozenset[str]]: The groups to dismiss — possibly none, and
            possibly both under ``raw_cutoff``.
        """
        support = {group: self.support[group] for group in pair}
        if anomaly_method == "raw_cutoff":
            return {group for group in pair if support[group] < threshold}
        if anomaly_method == "two_by_two":
            weakest = min(support.values())
            candidates = [group for group in pair if support[group] == weakest]
            if len(candidates) == 1:
                return set(candidates)
            # Equal support: the errors do not decide. Draw one side rather
            # than let siglum order pretend to.
            self.ties.append(pair)
            if self.seed is None:
                return {min(candidates, key=sorted)}
            return {self.rng.choice(candidates)}
        if anomaly_method == "probabilistic_cutoff":
            total = sum(support.values())
            return {group for group in pair
                    if support[group] / total < threshold}
        return set()  # "none": every overlap survives as a reticulation

    def _resolve_anomalies(self,
                           anomaly_method: str,
                           threshold: float) -> None:
        """Dismiss the contaminated groups and rebuild the lattice."""
        if self.verbose and self.overlaps:
            self.report_contamination()

        dismissed = set()
        for pair in self.overlaps:
            dismissed |= self._dismiss_in_pair(pair, anomaly_method, threshold)
        if not dismissed:
            return

        if self.verbose:
            for group in sorted(dismissed, key=sorted):
                self.console.print(
                    f"[red]✗ écarté (polygenèse) : {self._sigla(group)} "
                    f"×{self.support[group]}[/red]")
            self._report_ties()
        self._compute_lattice(dismissed)

    def _report_ties(self) -> None:
        """Warn about the conflicts the errors did not decide."""
        for first, second in self.ties:
            drawn = "tirage au sort" if self.seed is not None \
                else "ordre des sigles"
            self.console.print(
                f"[yellow]≈ égalité : {self._sigla(first)} et "
                f"{self._sigla(second)} à ×{self.support[first]} chacun — "
                f"départagés par {drawn}, non par les fautes.[/yellow]")

    # ------------------------------------------------------------------ #
    # Tree building
    # ------------------------------------------------------------------ #
    def build_tree(self,
                   root_witness: str | None = None,
                   anomaly_threshold: float = 2,
                   anomaly_method: str = "raw_cutoff") -> dict[str, dict]:
        """Reconstruct the stemma as a nested dict.

        Args:
            root_witness (str): Siglum placed at the root. Defaults to the
                reference witness.
            anomaly_threshold (float): Support below which a contaminated group
                is treated as polygenesis. A locus count for ``raw_cutoff``
                (default 2: a one-locus agreement loses a conflict; use 1 to
                keep everything), a probability in ]0, 1[ for
                ``probabilistic_cutoff``, ignored otherwise.
            anomaly_method (str): How contaminated groups are resolved:

                - ``raw_cutoff``: drop every conflicting group attested at
                  fewer than ``anomaly_threshold`` loci;
                - ``two_by_two``: in each conflicting pair, drop the group with
                  fewer loci. On equal support the errors decide nothing: the
                  loser is drawn at random when a ``seed`` was given, and taken
                  in siglum order otherwise. See :attr:`ties`;
                - ``probabilistic_cutoff``: drop a group whose share of the
                  pair's loci falls below ``anomaly_threshold`` (the cut-off
                  Froger describes himself);
                - ``none``: keep every overlap, to be read as a reticulation.

        Returns:
            dict[str, dict]: The stemma, ``{root: {child: {...}, ...}}``.

        Raises:
            ValueError: If ``anomaly_method`` is unknown, or if
                ``probabilistic_cutoff`` is given a threshold outside ]0, 1[.
        """
        if anomaly_method not in ANOMALY_METHODS:
            raise ValueError(f"unknown anomaly_method {anomaly_method!r}; "
                             f"expected one of {ANOMALY_METHODS}.")
        if (anomaly_method == "probabilistic_cutoff"
                and not 0 < anomaly_threshold < 1):
            raise ValueError(
                "probabilistic_cutoff compares a support ratio in ]0, 1[; got "
                f"{anomaly_threshold!r}. Pass a probability, not a count.")

        root_witness = root_witness or self.reference
        self.reticulations = []
        self.node_groups = {}
        self.ties = []

        self._compute_lattice()
        self._resolve_anomalies(anomaly_method, anomaly_threshold)

        if not self.children_of:  # no error group beyond the universe
            self.tree = {root_witness: {}}
            self.node_groups = {root_witness: self.witnesses}
            return self.tree

        root_group = max(self.children_of, key=len)
        self.tree = self._grow_root(root_group, root_witness)
        self._cut_reticulations()
        return self.tree

    def _loose_witnesses(self, group: frozenset[str]) -> frozenset[str]:
        """The witnesses of a group that no child group accounts for."""
        children = self.children_of.get(group, [])
        covered = frozenset().union(*children) if children else frozenset()
        return group - covered

    def _grow_root(self,
                   root_group: frozenset[str],
                   root_witness: str) -> dict[str, dict]:
        """Build the tree from the root, labelled with the reference witness.

        Every other loose witness of the root group therefore hangs off it as a
        leaf instead of naming the node.
        """
        # Bookkeeping filled during the recursion, read by
        # _cut_reticulations().
        self._label_of: dict[frozenset[str], str] = {}
        self._children_of_group: dict[frozenset[str], dict[str, dict]] = {}
        self._group_named_by: dict[str, frozenset[str]] = {}
        self._leaf_copies: dict[str, list[frozenset[str]]] = defaultdict(list)
        self._free_labels = self._label_generator()

        children: dict[str, dict] = {}
        for witness in sorted(self._loose_witnesses(root_group)
                              - {root_witness}):
            children[witness] = {}
            self._leaf_copies[witness].append(root_group)
        for child_group in self.children_of.get(root_group, []):
            children.update(self._grow(child_group))

        self._label_of[root_group] = root_witness
        self._children_of_group[root_group] = children
        self._group_named_by[root_witness] = root_group
        return {root_witness: children}

    def _grow(self, group: frozenset[str]) -> dict[str, dict]:
        """Build the subtree of one error group.

        A single loose witness names the node: it is the group's own ancestor.
        With none or several, no witness can play that part, the node takes a
        hyparchetype letter and the loose witnesses hang off it as leaves.
        """
        if len(group) == 1:
            (witness,) = group
            self._label_of[group] = witness
            self._children_of_group[group] = {}
            self._group_named_by[witness] = group
            return {witness: {}}

        loose = self._loose_witnesses(group)
        children: dict[str, dict] = {}
        if len(loose) != 1:
            for witness in sorted(loose):
                children[witness] = {}
                self._leaf_copies[witness].append(group)
        for child_group in self.children_of.get(group, []):
            children.update(self._grow(child_group))

        if len(loose) == 1:
            (label,) = loose
            self._group_named_by[label] = group
        else:
            label = next(self._free_labels)

        self._label_of[group] = label
        self._children_of_group[group] = children
        return {label: children}

    def _label_generator(self):
        """Yield hyparchetype labels (α, β, …, α2, …) that no siglum uses."""
        index = 0
        while True:
            letter = ALPHABET[index % len(ALPHABET)]
            cycle = index // len(ALPHABET)
            label = letter if cycle == 0 else f"{letter}{cycle + 1}"
            index += 1
            if label not in self.witnesses:
                yield label

    def _cut_reticulations(self) -> None:
        """Keep each contaminated witness at one place, record the other edges.

        A witness left in two overlapping groups is emitted once per group, and
        a tree can hold only one of those positions. We keep the best-attested
        one — or the group the witness names, since it carries its subtree —
        and turn the others into reticulation edges, which :meth:`to_dot`,
        :meth:`to_networkx` and :meth:`to_extended_newick` put back.
        """
        reticulations = []
        for witness, groups in self._leaf_copies.items():
            kept = self._group_named_by.get(witness)
            if kept is None:
                if len(groups) == 1:
                    continue
                kept = max(groups,
                           key=lambda group: (self.support.get(group, 0),
                                              len(group), sorted(group)))
            for group in groups:
                if group != kept:
                    self._children_of_group[group].pop(witness, None)
                    reticulations.append((self._label_of[group], witness))

        self.reticulations = sorted(set(reticulations))
        self.node_groups = {label: group
                            for group, label in self._label_of.items()}

    @property
    def is_reticulate(self) -> bool:
        """True when a contamination survived: the stemma is not a tree.

        :meth:`to_newick` then silently drops one branch; use :meth:`to_dot`,
        :meth:`to_networkx` or :meth:`to_extended_newick` instead.
        """
        return bool(self.reticulations)

    def _ensure_tree(self, root_witness: str | None = None) -> None:
        """Build the tree with the default settings if it is not built yet."""
        if self.tree is None:
            self.build_tree(root_witness)

    def _all_node_groups(self) -> dict[str, frozenset[str]]:
        """Node labels mapped to their group, bare leaf witnesses included."""
        groups = dict(self.node_groups)
        for witness in self.witnesses:
            groups.setdefault(witness, frozenset({witness}))
        return groups

    def _tree_edges(self) -> list[tuple[str, str]]:
        """The tree edges, as ``(parent_label, child_label)`` pairs."""
        edges = []

        def walk(label, children):
            for child, grandchildren in children.items():
                edges.append((label, child))
                walk(child, grandchildren)

        (root_label, root_children), = self.tree.items()
        walk(root_label, root_children)
        return edges

    def _sigla(self, group) -> str:
        """Format a group as a sorted, comma-separated list of sigla."""
        return ", ".join(sorted(group, key=natural_key))

    def strength(self, support: int) -> tuple[str, str]:
        """Grade how well attested a group is.

        Args:
            support (int): The number of loci supporting the group.

        Returns:
            tuple[str, str]: The ``(marker, name)`` of the matching band.
        """
        for upper_bound, marker, name in self.support_bands:
            if support < upper_bound:
                return marker, name
        return self.support_bands[-1][1], self.support_bands[-1][2]

    def support_legend(self) -> str:
        """Return the reading key of the strength markers, e.g. ``🔴 <10``."""
        keys, lower_bound = [], 0
        for upper_bound, marker, name in self.support_bands:
            if upper_bound == float("inf"):
                key = f"{marker} ≥{lower_bound} {name}"
            elif lower_bound:
                key = f"{marker} {lower_bound}–{upper_bound - 1} {name}"
            else:
                key = f"{marker} <{upper_bound} {name}"
            keys.append(key)
            lower_bound = (int(upper_bound)
                           if upper_bound != float("inf") else 0)
        return "   ".join(keys)

    # ------------------------------------------------------------------ #
    # Exports
    # ------------------------------------------------------------------ #
    def to_newick(self, root_witness: str | None = None) -> str:
        """Serialise the stemma to Newick, e.g. ``((C,D)β)A;``.

        Reticulations are dropped: see :meth:`to_extended_newick`.

        Args:
            root_witness (str): Siglum placed at the root.
        """
        self._ensure_tree(root_witness)

        def serialise(label, children):
            if not children:
                return label
            inner = ",".join(serialise(child, grandchildren)
                             for child, grandchildren in children.items())
            return f"({inner}){label}"

        (root_label, root_children), = self.tree.items()
        return serialise(root_label, root_children) + ";"

    def to_extended_newick(self, root_witness: str | None = None) -> str:
        """Serialise to extended Newick, tagging contaminated witnesses.

        A contaminated witness keeps its subtree under its genealogical parent
        and reappears as a bare ``#Hk`` leaf under each contaminating parent,
        e.g. ``((A,B#H1)β,(#H1,C)γ)ref;``.

        Args:
            root_witness (str): Siglum placed at the root.
        """
        self._ensure_tree(root_witness)
        tags = {witness: index for index, witness in enumerate(
            sorted({witness for _, witness in self.reticulations}), start=1)}
        contaminating = defaultdict(list)
        for parent, witness in self.reticulations:
            contaminating[parent].append(witness)

        def serialise(label, children):
            parts = [serialise(child, grandchildren)
                     for child, grandchildren in children.items()]
            parts += [f"#H{tags[witness]}"
                      for witness in sorted(contaminating[label])]
            name = f"{label}#H{tags[label]}" if label in tags else label
            return f"({','.join(parts)}){name}" if parts else name

        (root_label, root_children), = self.tree.items()
        return serialise(root_label, root_children) + ";"

    def to_networkx(self, root_witness: str | None = None):
        """Return the stemma as a :class:`networkx.DiGraph`.

        Nodes carry ``kind`` (``witness`` or ``hyparchetype``), ``group``
        (the witnesses they cover), ``support`` (their locus count, 0 for a
        leaf) and ``contaminated``. Edges carry ``kind``: ``genealogical`` for
        the tree edges, ``contamination`` for the reticulations.

        Args:
            root_witness (str): Siglum placed at the root.

        Returns:
            networkx.DiGraph: The stemma, reticulations included.
        """
        import networkx as nx

        self._ensure_tree(root_witness)
        contaminated = {witness for _, witness in self.reticulations}

        graph = nx.DiGraph()
        for label, group in self._all_node_groups().items():
            graph.add_node(
                label,
                kind="witness" if label in self.witnesses else "hyparchetype",
                group=group,
                support=self.support.get(group, 0),
                contaminated=label in contaminated)
        graph.add_edges_from(self._tree_edges(), kind="genealogical")
        graph.add_edges_from(self.reticulations, kind="contamination")
        return graph

    def to_dot(self,
               root_witness: str | None = None,
               show_support: bool = True) -> str:
        """Serialise the stemma to Graphviz DOT, ready for ``dot -Tsvg``.

        Witnesses are boxes, hyparchetypes circles, contamination edges dashed
        and red.

        Args:
            root_witness (str): Siglum placed at the root.
            show_support (bool): Append ``×N`` (locus count) to the nodes.
        """
        self._ensure_tree(root_witness)
        contaminated = {witness for _, witness in self.reticulations}

        lines = ["digraph stemma {",
                 "    rankdir=TB;",
                 '    node [fontname="serif"];']
        for label, group in self._all_node_groups().items():
            support = self.support.get(group) if show_support else None
            text = f"{label} ×{support}" if support else label
            shape = ("shape=box, style=rounded" if label in self.witnesses
                     else "shape=circle, width=0.3")
            attributes = f'label="{text}", {shape}'
            if label in contaminated:
                attributes += ", color=red, penwidth=2"
            lines.append(f'    "{label}" [{attributes}];')

        for parent, child in self._tree_edges():
            lines.append(f'    "{parent}" -> "{child}";')
        for parent, witness in self.reticulations:
            lines.append(f'    "{parent}" -> "{witness}" '
                         "[style=dashed, color=red, constraint=false];")
        lines.append("}")
        return "\n".join(lines)

    # ------------------------------------------------------------------ #
    # Visualisation
    # ------------------------------------------------------------------ #
    def render_tree(self,
                    root_witness: str | None = None,
                    show_groups: bool = False) -> str:
        """Return the stemma drawn with box-drawing characters, root on top.

        Args:
            root_witness (str): Siglum placed at the root.
            show_groups (bool): Annotate each internal node with its error
                group and the locus count supporting it, as
                ``label (w1, w2, …) ×N``. A node whose group is not an error
                group — the root, which holds the reference — shows no count.
                The group comes from :attr:`node_groups`, so a node keeps its
                full group even when a contaminated witness was moved out of
                its children.

        Returns:
            str: The rendered stemma, one node per line.
        """
        self._ensure_tree(root_witness)

        def draw(label, children, prefix="", is_last=True, is_root=True):
            text = label
            if children and show_groups:
                group = self.node_groups.get(label, frozenset())
                support = self.support.get(group)
                count = ""
                if support:
                    marker, _ = self.strength(support)
                    count = f"  ×{support} {marker}"
                text = f"{label}  ({self._sigla(group)}){count}"
            if is_root:
                lines, child_prefix = [text], ""
            else:
                lines = [prefix + ("└── " if is_last else "├── ") + text]
                child_prefix = prefix + ("    " if is_last else "│   ")
            items = list(children.items())
            for position, (child, grandchildren) in enumerate(items):
                lines += draw(child, grandchildren, child_prefix,
                              position == len(items) - 1, False)
            return lines

        (root_label, root_children), = self.tree.items()
        return "\n".join(draw(root_label, root_children))

    def display_tree(self,
                     root_witness: str | None = None,
                     show_groups: bool = True) -> str:
        """Print the stemma in a panel and return the rendered string.

        Args:
            root_witness (str): Siglum placed at the root.
            show_groups (bool): Annotate nodes with group and support.
        """
        drawing = self.render_tree(root_witness, show_groups)
        self.console.print(Panel(
            drawing,
            title=("Stemma codicum — contaminé" if self.is_reticulate
                   else "Stemma codicum"),
            subtitle=self.support_legend() if show_groups else None,
            border_style="yellow" if self.is_reticulate else "green",
            box=box.ROUNDED))
        if self.is_reticulate:
            witnesses = self._sigla({w for _, w in self.reticulations})
            self.console.print(
                f"[yellow]⚯ {witnesses} : deux ancêtres, stemma non "
                "arborescent. to_newick() en perd une branche ; préférer "
                "to_dot(), to_networkx() ou to_extended_newick().[/yellow]")
        return drawing

    def report_contamination(self) -> None:
        """Print each contaminated pair with its readings, for examination."""
        if not self.overlaps:
            self.console.print(
                "[green]✓ Aucun croisement de groupes.[/green]")
            return

        self.console.print(
            "\n[bold yellow]⚠ Contamination — leçons à examiner[/bold yellow]")
        for first, second in self.overlaps:
            self.console.print(
                f"\n[cyan]Groupes croisés[/cyan] "
                f"(témoin(s) partagé(s) : {self._sigla(first & second)})")
            table = Table(box=box.ROUNDED, header_style="bold magenta")
            table.add_column("Groupe", style="yellow")
            table.add_column("Locus", style="dim")
            table.add_column("Leçon du groupe", style="cyan")
            table.add_column(f"Réf. ({self.reference})", style="green")
            for group in (first, second):
                info = self.group_info[group]
                label = f"{self._sigla(group)} ×{info['frequency']}"
                for locus, reading, reference in zip(
                        info["loci"], info["readings"], info["reference"]):
                    table.add_row(label, locus, reading, reference)
                    label = ""  # the siglum is only printed on the first line
            self.console.print(table)

    def display_contamination(self) -> None:
        """Report the surviving overlaps, the dismissals and the edges."""
        if not (self.overlaps or self.dismissed or self.reticulations
                or self.ties):
            self.console.print(
                "[green]✓ Aucune contamination détectée.[/green]")
            return

        self.report_contamination()
        if self.dismissed:
            self.console.print("\n[red]Écartés (polygenèse) :[/red]")
            for group in sorted(self.dismissed, key=sorted):
                self.console.print(
                    f"  ✗ {self._sigla(group)} ×{self.support[group]}")
        if self.ties:
            self.console.print("\n[yellow]Conflits non tranchés par les "
                               "fautes :[/yellow]")
            self._report_ties()
        if self.reticulations:
            self.console.print(
                "\n[yellow]Réticulations (contamination retenue) :[/yellow]")
            for parent, witness in self.reticulations:
                self.console.print(f"  {parent} ⇢ {witness}")

    def display_variant_groups(self, min_group_size: int = 2) -> None:
        """Tabulate the readings of each error group against the reference.

        Args:
            min_group_size (int): Skip the smaller groups. Defaults to 2, which
                hides the singular errors (*lectiones singulares*).
        """
        for group, info in sorted(self.group_info.items(),
                                  key=lambda item: -item[1]["frequency"]):
            if len(group) < min_group_size:
                continue
            status = " [red](écarté)[/red]" if group in self.dismissed else ""
            marker, name = self.strength(info["frequency"])
            self.console.print(
                f"\n[bold cyan]{self._sigla(group)}[/bold cyan] "
                f"[dim]×{info['frequency']} loci — {marker} {name}[/dim]"
                f"{status}")
            table = Table(box=box.ROUNDED, header_style="bold magenta")
            table.add_column("Locus", style="dim")
            table.add_column("Leçon du groupe", style="cyan")
            table.add_column(f"Réf. ({self.reference})", style="green")
            for locus, reading, reference in zip(
                    info["loci"], info["readings"], info["reference"]):
                table.add_row(locus, reading, reference)
            self.console.print(table)

    def display_summary(self) -> None:
        """Print an overview of the reconstruction."""
        table = Table(title="Stemma — récapitulatif", box=box.ROUNDED)
        table.add_column("Propriété", style="cyan")
        table.add_column("Valeur", style="white")
        table.add_row("Témoins", str(len(self.witnesses)))
        table.add_row("Témoin de référence", self.reference)
        table.add_row("Groupes de fautes", str(len(self.support)))
        if self.dismissed:
            table.add_row("Groupes écartés", str(len(self.dismissed)))
            table.add_row("Groupes retenus",
                          str(len(self.support) - len(self.dismissed)))
        if self.children_of:
            table.add_row("Relations de couverture", str(len(self.edges)))
        if self.ties:
            table.add_row(
                "Conflits départagés au hasard" if self.seed is not None
                else "Conflits départagés arbitrairement",
                f"[yellow]{len(self.ties)}"
                + (f" (seed={self.seed})[/yellow]" if self.seed is not None
                   else "[/yellow]"))
        contaminated = self.contaminated_witnesses
        table.add_row(
            "Témoins contaminés",
            f"[red]{self._sigla(contaminated)}[/red]" if contaminated
            else "[green]aucun[/green]")
        self.console.print(table)

        if self.support:
            self.console.print("\n[bold]Groupes par fréquence[/bold]")
            counts = Counter(self.support.values())
            frequency_table = Table(box=box.MINIMAL)
            frequency_table.add_column("Loci", style="cyan", justify="right")
            frequency_table.add_column("Groupes", style="green",
                                       justify="center")
            for frequency in sorted(counts):
                frequency_table.add_row(str(frequency), str(counts[frequency]))
            self.console.print(frequency_table)
