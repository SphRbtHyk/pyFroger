"""Given the list apparatus, build an unoriented stemma using the group method."""
from collections import Counter
import re
import warnings
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from pyfroger.list_apparatus import ListApparatus
from rich import box


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

    # Additional styles
    INFO = '\033[94m'
    SUCCESS = '\033[92m'
    ERROR = '\033[91m'
    HIGHLIGHT = '\033[103m\033[30m'  # Yellow background, black text


class UnorientedStemma:
    """Build an unoriented stemma using the group method
    ('Méthode par les groupes').

    Attributes:
        list_apparatus: The variant groups and their different information.
        variant_frequencies: Variant groups keyed by frequency.
        members: Set of witness sigla available to build the tree.
        edges: Covering-relation edges once computed, or None.
        children_of: Direct parent→children mapping once computed, or None.
        tree: Nested-dict stemma once built, or None.
        contamination_overlap: Pairs of groups that partially overlap
            (horizontal transmission), or None.
        contaminated_witnesses: Frozenset of witnesses involved in
            contamination, or None.
    """

    def __init__(self,
                 list_apparatus: ListApparatus,
                 min_frequency: int = 1,
                 verbose: bool = False) -> None:
        """Initialise the stemma builder from a variant apparatus.

        Args:
            list_apparatus (ListApparatus): The variant groups and their different information.
            min_frequency (int): The minimum frequency of apparition for the group to be
                used for building the stemma.

        """
        self.min_frequency = min_frequency
        self.list_apparatus: ListApparatus = list_apparatus
        self.variant_frequencies =\
            self.list_apparatus.group_by_frequency(min_frequency=min_frequency)
        self.members: set[str] = list_apparatus.witnesses
        self.edges: list[tuple[frozenset[str], frozenset[str]]] | None = None
        self.children_of: dict[frozenset[str],
                               list[frozenset[str]]] | None = None
        self.tree: dict[str, dict] | None = None
        self.reticulations: list[tuple[str, str]] = []
        self.node_groups: dict[str, frozenset[str]] = {}
        self.contamination_overlap: (
            list[tuple[frozenset[str], ...]] | None
        ) = None
        self.contaminated_witnesses: frozenset[str] | None = None
        self.dismissed: set[frozenset[str]] = set()
        self.verbose = verbose
        # Setup Rich console
        self.console = Console()

    @staticmethod
    def _compute_covering_relations(
        groups: frozenset[frozenset[str]],
    ) -> tuple[
        list[tuple[frozenset[str], frozenset[str]]],
        dict[frozenset[str], list[frozenset[str]]],
        list[tuple[frozenset[str], ...]],
    ]:
        """Compute covering relations between variant groups.

        Two groups are related by a *covering relation* when one strictly
        contains the other with no intermediate group between them in the
        inclusion order. Groups that partially overlap signal potential
        horizontal transmission (contamination).

        Args:
            groups: All variant groups expressed as a frozenset of frozensets
                of witness sigla.

        Returns:
            A three-element tuple:

            - **edges**: Sorted list of ``(parent, child)`` covering pairs.
            - **children_of**: Mapping from each parent group to its list of
              direct child groups, sorted lexicographically.
            - **contamination_overlap**: List of ``(g1, g2)`` pairs whose
              intersection is non-empty but neither group contains the other.
        """
        def covers(c: frozenset[str], p: frozenset[str]) -> bool:
            return c < p and not any(c < u < p for u in groups)

        contamination_overlap: set[frozenset[frozenset[str]]] = set()
        edges: list[tuple[frozenset[str], frozenset[str]]] = []
        children_of: dict[frozenset[str], list[frozenset[str]]] = {}

        for c in groups:
            for p in groups:
                if covers(c, p):
                    edges.append((p, c))
                    children_of.setdefault(p, []).append(c)
                else:
                    shared = c & p
                    if shared and shared != c and shared != p:
                        contamination_overlap.add(frozenset((c, p)))

        return (
            sorted(edges),
            {k: sorted(v, key=lambda s: sorted(s))
             for k, v in children_of.items()},
            [
                tuple(sorted(pair, key=lambda s: sorted(s)))
                for pair in contamination_overlap
            ],
        )

    def _lattice_groups(
        self,
        dismissed: frozenset[frozenset[str]] = frozenset(),
    ) -> frozenset[frozenset[str]]:
        """Witness-set groups used to build the lattice.

        Keeps only error groups attested at >= self.min_frequency loci
        (problematic low-attestation groups are dropped up front), minus any
        explicitly dismissed groups, plus the universe (all witnesses + ref).
        """
        universe = frozenset(self.members | {self.list_apparatus.reference})
        kept = frozenset(
            frozenset(g)
            for g in self.list_apparatus.group_by_frequency(
                min_frequency=self.min_frequency)
            if frozenset(g) not in dismissed
        )
        return kept | {universe}

    def compute_edges(self) -> list[tuple[frozenset[str], frozenset[str]]]:
        """Compute covering relations between error groups and store results."""
        groups = self._lattice_groups()
        self.edges, self.children_of, self.contamination_overlap = (
            self._compute_covering_relations(groups)
        )
        self.contaminated_witnesses = frozenset(
            w for g1, g2 in self.contamination_overlap for w in (g1 & g2)
        )
        return self.edges

    def _rebuild_without(self, dismissed: set[frozenset[str]]) -> None:
        """Recompute covering relations after excluding dismissed groups."""
        groups = self._lattice_groups(dismissed=frozenset(dismissed))
        self.edges, self.children_of, self.contamination_overlap = (
            self._compute_covering_relations(groups)
        )
        self.contaminated_witnesses = frozenset(
            w for g1, g2 in self.contamination_overlap for w in (g1 & g2)
        )
        self.dismissed = set(dismissed)
        self.tree = None
        self.reticulations = []
        self.node_groups = {}

    def _build_tree_structure(
        self,
        children_of: dict[frozenset[str], list[frozenset[str]]],
        root_group: frozenset[str],
        root_witness: str,
        reserved: frozenset[str] = frozenset(),
        support: dict[frozenset[str], int] | None = None,
    ) -> dict[str, dict]:
        """Recursively build a nested-dict stemma from covering relations.

        Each key in the returned dict is a node label (witness siglum or
        Greek-letter hypothetical ancestor); its value is a dict of the same
        shape representing its children.

        In the group method each error group corresponds to a set of witnesses
        that share a common ancestor.  Witnesses in a group that are *not*
        subsumed by any child group sit directly under the node for that group:

        - If exactly one such witness exists it becomes the node label (the
          identifiable ancestor of the group at that level).
        - If zero or more than one exist, no single witness can play that role
          and a Greek-letter hyperarchetype (α, β, γ, …) is introduced instead.

        Args:
            children_of: Mapping from each parent group to its list of direct
                child groups.
            root_group: The largest group, used as the starting point for
                recursion.
            root_witness: Siglum placed at the root of the returned tree
                (in our case the selected reference witness).

        Returns:
            A nested dict ``{root_witness: {child_label: {…}, …}}``.
        """
        support = support or {}
        alphabet = "αβγδεζηθικλμνξοπρστυφχψω"
        idx = 0

        # Bookkeeping needed to turn the tree into a graph: which label stands
        # for which group, which dict holds a node's children, where each
        # witness was emitted as a leaf, and which group a witness labels.
        label_of: dict[frozenset[str], str] = {}
        node_children: dict[frozenset[str], dict[str, dict]] = {}
        emitted: dict[str, list[frozenset[str]]] = {}
        owner: dict[str, frozenset[str]] = {}

        def next_label() -> str:
            """Return a fresh hyperarchetype label that collides neither with
            a witness siglum nor with a previously issued label."""
            nonlocal idx
            while True:
                letter = alphabet[idx % len(alphabet)]
                cycle = idx // len(alphabet)
                label = letter if cycle == 0 else f"{letter}{cycle + 1}"
                idx += 1
                if label not in reserved:
                    return label

        def build(group: frozenset[str]) -> dict[str, dict]:
            if len(group) == 1:
                witness = next(iter(group))
                label_of[group] = witness
                node_children[group] = {}
                owner[witness] = group
                return {witness: {}}

            subs = children_of.get(group, [])
            covered = frozenset().union(*subs) if subs else frozenset()
            complement = group - covered
            merged: dict[str, dict] = {}

            if len(complement) != 1:
                for w in sorted(complement):
                    merged[w] = {}
                    emitted.setdefault(w, []).append(group)

            for sg in subs:
                merged.update(build(sg))

            if len(complement) == 1:
                label = next(iter(complement))
                owner[label] = group
            else:
                label = next_label()

            label_of[group] = label
            node_children[group] = merged
            return {label: merged}

        subs = children_of.get(root_group, [])
        covered = frozenset().union(*subs) if subs else frozenset()
        root_children: dict[str, dict] = {}
        for w in sorted((root_group - covered) - {root_witness}):
            root_children[w] = {}
            emitted.setdefault(w, []).append(root_group)

        for sg in subs:
            root_children.update(build(sg))

        label_of[root_group] = root_witness
        node_children[root_group] = root_children
        owner[root_witness] = root_group

        # A witness sitting in two overlapping (non-laminar) groups is emitted
        # once per group. A tree can hold only one of those positions: keep the
        # best-attested one (or the group the witness labels, which carries its
        # subtree) and record the others as reticulation edges, so that the
        # graph exports can put them back.
        reticulations: list[tuple[str, str]] = []
        for witness, groups in emitted.items():
            primary = owner.get(witness)
            if primary is None:
                if len(groups) == 1:
                    continue
                primary = max(
                    groups,
                    key=lambda g: (support.get(g, 0), len(g), sorted(g)),
                )
            for group in groups:
                if group == primary:
                    continue
                node_children[group].pop(witness, None)
                reticulations.append((label_of[group], witness))

        self.node_groups = {lab: g for g, lab in label_of.items()}
        self.reticulations = sorted(set(reticulations))
        if self.reticulations:
            warnings.warn(
                "Witness(es) {} occur in overlapping (non-laminar) groups and "
                "were kept at a single position in the tree; the overlap was "
                "not resolved upstream. Use anomaly_method='two_by_two' or "
                "lower the cutoff to resolve it, or to_dot() / "
                "to_extended_newick() to keep it as a reticulation.".format(
                    ", ".join(sorted({w for _, w in self.reticulations}))),
                stacklevel=2,
            )

        return {root_witness: root_children}

    @staticmethod
    def _to_newick_string(tree: dict[str, dict]) -> str:
        """Serialize a nested-dict tree to Newick format.

        Args:
            tree: Nested dict as produced by :meth:`_build_tree_structure`.

        Returns:
            Newick string ending with ``;``, e.g. ``(A,(C,D)B)ω;``.
        """
        def serialize(children: dict[str, dict], label: str) -> str:
            if not children:
                return label
            inner = ",".join(
                serialize(grandchildren, child)
                for child, grandchildren in children.items()
            )
            return f"({inner}){label}"

        root_label, root_children = next(iter(tree.items()))
        return serialize(root_children, root_label) + ";"

    def to_newick(self, root_witness: str | None = None) -> str:
        """Return the stemma serialized in Newick format.

        Calls :meth:`build_tree` first if it has not been called yet.

        Args:
            root_witness: Siglum to place at the root.  Defaults to the
                reference witness defined in ``list_apparatus``.

        Returns:
            Newick string ending with ``;``, e.g. ``(A,(C,D)B)ω;``.
        """
        if self.tree is None:
            self.build_tree(root_witness)
        assert self.tree is not None
        return self._to_newick_string(self.tree)

    @staticmethod
    def _render_tree_with_groups(
        tree: dict[str, dict],
        witnesses: frozenset[str],
        frequencies: dict[frozenset[str], int] | None = None,
    ) -> str:
        """Render a nested-dict stemma showing each internal node's group and,
        optionally, the number of loci supporting it.

        Each node is annotated with the witnesses appearing anywhere in its
        subtree — an internal witness-ancestor counted in its own group. The
        ``witnesses`` set distinguishes real sigla from Greek hyperarchetype
        labels, which are excluded from the groups. When ``frequencies`` is
        given, a node whose group is an error group also shows ``×N`` (its
        locus count); nodes whose group is absent from the mapping — notably
        the root, whose group includes the reference — show no count. Leaves
        are left bare.

        Args:
            tree: Nested dict as produced by :meth:`_build_tree_structure`.
            witnesses: Full set of witness sigla (anything else is a
                hyperarchetype label).
            frequencies: Optional mapping from an error group to its number of
                supporting loci, e.g. ``{g: info["frequency"]}`` from
                ``group_by_frequency``.

        Returns:
            A multi-line box-drawing string with ``label  (w1, w2, …)  ×N``
            annotations on every node that has children.
        """
        frequencies = frequencies or {}

        def nat(s: str) -> list:
            return [int(t) if t.isdigit() else t
                    for t in re.split(r"(\d+)", s)]

        def group_of(label: str, children: dict[str, dict]) -> set[str]:
            g = {label} if label in witnesses else set()
            for child, grandchildren in children.items():
                g |= group_of(child, grandchildren)
            return g

        def render(
            label: str,
            children: dict[str, dict],
            prefix: str = "",
            is_last: bool = True,
            is_root: bool = True,
        ) -> list[str]:
            if children:
                group = group_of(label, children)
                sigla = ", ".join(sorted(group, key=nat))
                count = frequencies.get(frozenset(group))
                suffix = f"  ×{count}" if count is not None else ""
                text = f"{label}  ({sigla}){suffix}"
            else:
                text = label
            if is_root:
                lines, cpfx = [text], ""
            else:
                lines = [prefix + ("└── " if is_last else "├── ") + text]
                cpfx = prefix + ("    " if is_last else "│   ")
            items = list(children.items())
            for i, (child, grandchildren) in enumerate(items):
                lines.extend(
                    render(child, grandchildren, cpfx,
                           i == len(items) - 1, False)
                )
            return lines

        root_label = next(iter(tree))
        return "\n".join(render(root_label, tree[root_label]))

    def display_groups(self, root_witness: str | None = None) -> str:
        """Return a box-drawing stemma annotating each node with its group and
        the number of loci supporting it.

        Like :meth:`display_tree`, but every internal node is labelled with the
        set of witnesses it covers and the locus count ``×N`` of the
        corresponding error group. Calls :meth:`build_tree` first if the tree
        has not been built yet.

        Args:
            root_witness: Siglum to place at the root. Defaults to the
                reference witness defined in ``list_apparatus``.

        Returns:
            Multi-line string suitable for printing to a terminal.
        """
        if self.tree is None:
            self.build_tree(root_witness)
        assert self.tree is not None
        witnesses = frozenset(self.members | {self.list_apparatus.reference})
        frequencies = {
            frozenset(g): info["frequency"]
            for g, info in self.list_apparatus.group_by_frequency(
                min_frequency=self.min_frequency).items()
        }
        return self._render_tree_with_groups(self.tree, witnesses, frequencies)

    def build_tree(self,
                   root_witness: str | None = None,
                   resolve_anomaly_threshold: float = 2,
                   anomaly_method: str = "raw_cutoff") -> dict[str, dict]:
        """Build a nested-dict stemma from the variant apparatus.

        Args:
            root_witness (str): Siglum to place at the root. Defaults to the
                reference witness defined in ``list_apparatus``.
            resolve_anomaly_threshold (int): Threshold below which the
                weaker side of an overlapping pair is treated as polygenesis
                and dismissed. Defaults to 2 (a single-site agreement in a
                conflict is dismissed); set to 1 to disable.
            anomaly_method (str): The method selected for resolving the
                anomalies. Can be either:
                    - raw_cutoff (all conflicting groups supported by fewer
                    than `resolve_anomaly_threshold` independent sites are
                    removed);
                    - two_by_two: Overlaps are resolved by a two by
                        two comparison and keeping the group with the higher
                        number of independent sites;
                    - probabilistic_cutoff: The resolve_anomaly threshold is
                        understood as a probability of contamination,
                        computed as the (number of variants supporting each
                        conflicting group) / (total variants
                                            in the two groups).
                        This cut-off was described by Froger himself.

        Returns:
            A nested dict ``{root: {child: {…}, …}}`` representing the stemma.
        """
        if self.children_of is None:
            self.compute_edges()
        assert self.children_of is not None
        if root_witness is None:
            root_witness = self.list_apparatus.reference

        if not self.children_of:
            # No covering relations (no error groups beyond the universe):
            # nothing to build but the root itself.
            self.tree = {root_witness: {}}
            return self.tree

        self._resolve(resolve_anomaly_threshold, anomaly_method)

        # Recomputed *after* anomaly resolution: _rebuild_without() replaces
        # self.children_of, so the root group must be read from the current
        # lattice, not the pre-resolution one.
        assert self.children_of is not None
        if not self.children_of:
            self.tree = {root_witness: {}}
            return self.tree

        root_group = max(self.children_of, key=len)

        self.tree = self._build_tree_structure(
            self.children_of,
            root_group,
            root_witness,
            reserved=frozenset(self.members | {self.list_apparatus.reference}),
            support={
                frozenset(g): info["frequency"]
                for g, info in self.list_apparatus.group_by_frequency(
                    min_frequency=self.min_frequency).items()
            },
        )

        return self.tree

    def _resolve(self,
                 resolve_anomaly_threshold: float,
                 anomaly_method: str) -> None:
        """Dispatch to the requested anomaly-resolution strategy.

        Shared by :meth:`build_tree` and :meth:`build_graph` so that both see
        the same lattice. Groups surviving resolution may still overlap; it is
        then up to the caller to prune them (tree) or keep them as
        reticulations (graph).
        """
        if anomaly_method == "raw_cutoff":
            self._resolve_anomaly_cutoff(
                min_loci=int(resolve_anomaly_threshold))
        elif anomaly_method == "two_by_two":
            self._resolve_anomaly_two_by_two()
        elif anomaly_method == "probabilistic_cutoff":
            if not 0 < resolve_anomaly_threshold < 1:
                raise ValueError(
                    "anomaly_method='probabilistic_cutoff' compares a support "
                    "ratio in ]0, 1[ against resolve_anomaly_threshold; got "
                    f"{resolve_anomaly_threshold!r}. Pass a probability "
                    "(e.g. 0.5), not a locus count."
                )
            self._resolve_anomaly_probabilistic(
                resolve_anomaly_threshold=resolve_anomaly_threshold)
        elif anomaly_method == "none":
            # Keep every group: overlaps are meant to survive as reticulations.
            pass
        else:
            raise ValueError(
                f"Unknown anomaly_method {anomaly_method!r}; expected one of "
                "'raw_cutoff', 'two_by_two', 'probabilistic_cutoff', 'none'."
            )

    @property
    def has_unresolved_overlaps(self) -> bool:
        """True when the built tree had to drop a witness position.

        Set by :meth:`build_tree`. When it holds, at least one witness descends
        from two ancestors: the stemma is not a tree, :meth:`to_newick` silently
        omits one of the two positions, and :meth:`to_dot` or
        :meth:`to_extended_newick` should be used instead.
        """
        return bool(self.reticulations)

    def _node_groups(self) -> dict[str, frozenset[str]]:
        """Map each node label to the error group it stands for.

        Recorded while building the tree, so a node keeps its full group even
        when a reticulate witness was pruned out of its children.
        """
        witnesses = self.members | {self.list_apparatus.reference}
        groups = dict(self.node_groups)
        for witness in witnesses:
            groups.setdefault(witness, frozenset({witness}))
        return groups

    def to_dot(self,
               root_witness: str | None = None,
               show_support: bool = True) -> str:
        """Return the stemma serialized in Graphviz DOT format.

        Unlike :meth:`to_newick`, this keeps the reticulations: a witness left
        in two overlapping error groups gets an edge from each, the weaker one
        drawn dashed and red as horizontal transmission. Calls
        :meth:`build_tree` first if it has not been called yet.

        Args:
            root_witness: Siglum to place at the root. Defaults to the
                reference witness defined in ``list_apparatus``.
            show_support: Append ``×N`` (locus count) to internal nodes.

        Returns:
            A DOT source string, ready for ``dot -Tsvg``.
        """
        if self.tree is None:
            self.build_tree(root_witness)
        assert self.tree is not None

        witnesses = self.members | {self.list_apparatus.reference}
        frequencies = {
            frozenset(g): info["frequency"]
            for g, info in self.list_apparatus.group_by_frequency(
                min_frequency=self.min_frequency).items()
        }
        groups = self._node_groups()
        hybrids = {w for _, w in self.reticulations}

        lines = ["digraph stemma {",
                 "    rankdir=TB;",
                 '    node [fontname="serif"];']

        for label, group in groups.items():
            count = frequencies.get(group) if show_support else None
            text = f"{label} ×{count}" if count is not None else label
            if label in witnesses:
                attrs = f'label="{text}", shape=box, style=rounded'
            else:
                attrs = f'label="{text}", shape=circle, width=0.3'
            if label in hybrids:
                attrs += ", color=red, penwidth=2"
            lines.append(f'    "{label}" [{attrs}];')

        def walk(label: str, children: dict[str, dict]) -> None:
            for child, grandchildren in children.items():
                lines.append(f'    "{label}" -> "{child}";')
                walk(child, grandchildren)

        root_label = next(iter(self.tree))
        walk(root_label, self.tree[root_label])

        for parent, witness in self.reticulations:
            lines.append(f'    "{parent}" -> "{witness}" '
                         "[style=dashed, color=red, constraint=false];")

        lines.append("}")
        return "\n".join(lines)

    def to_extended_newick(self, root_witness: str | None = None) -> str:
        """Return the stemma serialized in extended Newick format.

        Same as :meth:`to_newick` when the stemma is a tree. Otherwise each
        reticulate witness carries a hybrid tag ``#Hk``: it appears with its
        subtree below its genealogical parent, and as a bare ``#Hk`` leaf below
        each contaminating parent. Calls :meth:`build_tree` first if it has not
        been called yet.

        Args:
            root_witness: Siglum to place at the root. Defaults to the
                reference witness defined in ``list_apparatus``.

        Returns:
            An eNewick string ending with ``;``, e.g.
            ``((A,B#H1)β,(#H1,C)γ)ref;``.
        """
        if self.tree is None:
            self.build_tree(root_witness)
        assert self.tree is not None

        tags = {
            witness: i
            for i, witness in enumerate(
                sorted({w for _, w in self.reticulations}), start=1)
        }
        extra: dict[str, list[str]] = {}
        for parent, witness in self.reticulations:
            extra.setdefault(parent, []).append(witness)

        def serialize(children: dict[str, dict], label: str) -> str:
            parts = [
                serialize(grandchildren, child)
                for child, grandchildren in children.items()
            ]
            parts += [f"#H{tags[w]}" for w in sorted(extra.get(label, []))]
            name = f"{label}#H{tags[label]}" if label in tags else label
            return f"({','.join(parts)}){name}" if parts else name

        root_label, root_children = next(iter(self.tree.items()))
        return serialize(root_children, root_label) + ";"

    def display_graph(self,
                      root_witness: str | None = None,
                      show_support: bool = True) -> str:
        """Print and return the stemma with its reticulations spelled out.

        The tree is drawn as usual, hybrid witnesses are marked ``⚯``, and the
        contamination edges dropped from the tree are listed underneath.

        Args:
            root_witness: Siglum to place at the root. Defaults to the
                reference witness defined in ``list_apparatus``.
            show_support: Append ``×N`` (locus count) to internal nodes.

        Returns:
            Multi-line string suitable for printing to a terminal.
        """
        if self.tree is None:
            self.build_tree(root_witness)
        assert self.tree is not None

        if show_support:
            body = self.display_groups(root_witness)
        else:
            body = self._render_tree(self.tree)

        hybrids = sorted({w for _, w in self.reticulations})
        lines = [body]
        if self.reticulations:
            lines.append("")
            lines.append("Contamination (horizontal transmission):")
            for parent, witness in self.reticulations:
                lines.append(f"  {parent} ⇢ {witness}")
        graph_str = "\n".join(lines)

        title = ("[bold yellow]Stemma (reticulate)[/bold yellow]"
                 if self.reticulations else "[bold green]Stemma[/bold green]")
        self.console.print(Panel(
            graph_str,
            title=title,
            border_style="yellow" if self.reticulations else "green",
            box=box.ROUNDED,
        ))
        if hybrids:
            self.console.print(
                f"[yellow]⚯ Reticulate witness(es): {', '.join(hybrids)} — "
                "not expressible in plain Newick; use to_extended_newick() "
                "or to_dot().[/yellow]\n"
            )
        return graph_str

    def _resolve_anomaly_probabilistic(
            self,
            resolve_anomaly_threshold: float):
        """Resolve anomalies by computing how well each group
        is supported as a ratio between the number of variant loci
        and the total number of variant readings in the overlap."""
        if self.children_of is None:
            self.compute_edges()

        dismissed = set()

        if not self.contamination_overlap:
            self.console.print(
                "[green]✓ No overlapping variant groups detected.[/green]")
            return

        if self.verbose:
            # Display contamination info in a nice table
            self.console.print(
                "\n[bold yellow]⚠ Please examine Contamination/Anomalies[/bold yellow]")

        freqs = self.list_apparatus.group_by_frequency(
            min_frequency=self.min_frequency)

        for overlap in self.contamination_overlap:
            # For each incompatible pair, count supporting sites
            locus_counts = {
                group: freqs[group]["frequency"]
                for group in overlap
            }

            if self.verbose:
                # Create a nice display for overlaps
                self.console.print("\n[cyan]Overlapping groups:[/cyan]")
                for group in overlap:
                    info = freqs[group]
                    group_str = f"  • {', '.join(sorted(group))}"
                    count_str = f" ({locus_counts[group]} loci)"
                    loci_str = f"\n    Loci: {', '.join(info['loci'])}"
                    readings_str = f"\n    Readings: {'\n'.join(info['readings'])}"
                    self.console.print(f"{group_str}{count_str}")
                    self.console.print(f"{loci_str}")
                    self.console.print(f"{readings_str}")

            # Check total number of variant per group
            total = 0
            for count in locus_counts.values():
                total += count

            if not total:
                continue

            # Compute probability of apparition per group
            locus_apparition = {}
            for group, count in locus_counts.items():
                locus_apparition[group] = count / total

            for group, probability in locus_apparition.items():
                if probability < resolve_anomaly_threshold:
                    # Dismiss group that are below the threshold
                    dismissed.add(group)
                    if self.verbose:
                        self.console.print(
                            "[red]✗ Dismissing group "
                            f"{', '.join(sorted(group))} "
                            f"(support ratio {probability:.2f} "
                            f"< {resolve_anomaly_threshold})[/red]"
                        )

        # Rebuild stemma without dismissed groups
        if dismissed:
            self.console.print(
                f"\n[yellow]Rebuilding stemma after dismissing {len(dismissed)} group(s)...[/yellow]")
            self._rebuild_without(dismissed)
            self.console.print("[green]✓ Stemma rebuilt successfully[/green]")

        return {
            "dismissed": dismissed,
        }

    def _resolve_anomaly_two_by_two(self):
        """Resolve anomalies by comparing the number of *independent* sites
        in two competing groups and dismissing the weaker.
        """
        if self.children_of is None:
            self.compute_edges()

        dismissed = set()

        if not self.contamination_overlap:
            self.console.print(
                "[green]✓ No overlapping variant groups detected.[/green]")
            return

        if self.verbose:
            # Display contamination info in a nice table
            self.console.print(
                "\n[bold yellow]⚠"
                "Please examine Contamination/Anomalies[/bold yellow]")

        freqs = self.list_apparatus.group_by_frequency(
            min_frequency=self.min_frequency)

        for overlap in self.contamination_overlap:
            # For each incompatible pair, count independent supporting sites
            # (contiguous loci collapse to a single transmission event).
            locus_counts = {
                group: len(freqs[group]["loci"])
                for group in overlap
            }
            if self.verbose:
                # Create a nice display for this overlap
                self.console.print("\n[cyan]Overlapping groups:[/cyan]")
                for group in overlap:
                    info = freqs[group]
                    group_str = f"  • {', '.join(sorted(group))}"
                    count_str = (f" ({info['frequency']} loci, "
                                 f"{locus_counts[group]} independent)")
                    loci_str = f"\n    Loci: {', '.join(info['loci'])}"
                    readings_str = f"\n    Readings: {'\n'.join(info['readings'])}"
                    self.console.print(f"{group_str}{count_str}")
                    self.console.print(f"{loci_str}")
                    self.console.print(f"{readings_str}")

            weakest = min(locus_counts, key=lambda g: locus_counts[g])
            dismissed.add(weakest)
            if self.verbose:
                self.console.print(
                    f"[red]✗ Dismissing weaker group {', '.join(sorted(weakest))} "
                    f"({locus_counts[weakest]} independent site(s))[/red]"
                )

        # Rebuild stemma without dismissed groups
        if dismissed:
            if self.verbose:
                self.console.print(
                    f"\n[yellow]Rebuilding stemma after dismissing {len(dismissed)} group(s)...[/yellow]")
            self._rebuild_without(dismissed)
            self.console.print("[green]✓ Stemma rebuilt successfully[/green]")

        return {
            "dismissed": dismissed,
        }

    def get_contamination_loci(self):
        """Get all readings that led to contamination."""
        contaminated_loci = []
        if self.contamination_overlap:
            for overlap in self.contamination_overlap:
                for group in overlap:
                    info = self.list_apparatus.group_by_frequency(
                        min_frequency=self.min_frequency)[group]
                    contaminated_loci.extend(info['loci'])
            return set(contaminated_loci)
        else:
            print("No contamination to check!")

    def _resolve_anomaly_cutoff(self, min_loci: int = 2):
        """Resolve anomalies by dismissing incompatible groups supported
        by fewer than `min_loci` *independent* sites as probable polygenesis.

        Args:
            min_loci (int, default 2): Minimum number of independent variation
                sites (contiguous loci collapsed) required for the weaker side
                of a conflict to be retained as genuine contamination evidence.
                Set to 1 to disable (never dismiss).
        """
        if self.children_of is None:
            self.compute_edges()

        dismissed = set()

        if not self.contamination_overlap:
            self.console.print(
                "[green]✓ No overlapping variant groups detected.[/green]")
            return

        if self.verbose:
            # Display contamination info in a nice table
            self.console.print(
                "\n[bold yellow]⚠ Please examine Contamination/Anomalies[/bold yellow]")

        freqs = self.list_apparatus.group_by_frequency(
            min_frequency=self.min_frequency)

        for overlap in self.contamination_overlap:
            # For each incompatible pair, count supporting sites
            locus_counts = {
                group: freqs[group]["frequency"]
                for group in overlap
            }

            if self.verbose:
                # Create a nice display for overlaps
                self.console.print("\n[cyan]Overlapping groups:[/cyan]")
                for group in overlap:
                    info = freqs[group]
                    group_str = f"  • {', '.join(sorted(group))}"
                    count_str = f" ({locus_counts[group]} loci)"
                    loci_str = f"\n    Loci: {', '.join(info['loci'])}"
                    readings_str = f"\n    Readings: {'\n'.join(info['readings'])}"
                    self.console.print(f"{group_str}{count_str}")
                    self.console.print(f"{loci_str}")
                    self.console.print(f"{readings_str}")

            # Check for weak groups
            for group in locus_counts:
                if locus_counts[group] < min_loci:
                    if self.verbose:
                        self.console.print(
                            "[red]✗ Dismissing weak group "
                            f"{', '.join(sorted(group))} "
                            f"(only {locus_counts[group]} site(s) "
                            f"< {min_loci})[/red]"
                        )
                    dismissed.add(group)

        # Rebuild stemma without dismissed groups
        if dismissed:
            self.console.print(
                f"\n[yellow]Rebuilding stemma after dismissing {len(dismissed)} group(s)...[/yellow]")
            self._rebuild_without(dismissed)
            self.console.print("[green]✓ Stemma rebuilt successfully[/green]")

        return {
            "dismissed": dismissed,
        }

    @staticmethod
    def _render_tree(tree: dict[str, dict]) -> str:
        """Render a nested-dict tree as a Unicode ASCII-art string.

        Args:
            tree: Nested dict as produced by :meth:`_build_tree_structure`.

        Returns:
            A multi-line string with the root at the top and branches drawn
            using ``└──`` / ``├──`` / ``│`` box-drawing characters.
        """
        def render(
            label: str,
            children: dict[str, dict],
            prefix: str = "",
            is_last: bool = True,
            is_root: bool = True,
        ) -> list[str]:
            lines: list[str] = []
            if is_root:
                lines.append(label)
                cpfx = ""
            else:
                lines.append(prefix + ("└── " if is_last else "├── ") + label)
                cpfx = prefix + ("    " if is_last else "│   ")
            items = list(children.items())
            for i, (child, grandchildren) in enumerate(items):
                lines.extend(
                    render(child, grandchildren, cpfx,
                           i == len(items) - 1, False)
                )
            return lines

        root_label = next(iter(tree))
        return "\n".join(render(root_label, tree[root_label]))

    def display_tree(self, root_witness: str | None = None, show_groups: bool = False) -> str:
        """Return a beautifully formatted Unicode ASCII-art representation of the stemma.

        Calls :meth:`build_tree` first if it has not been called yet.

        Args:
            root_witness: Siglum to place at the root. Defaults to the
                reference witness defined in ``list_apparatus``.
            show_groups: If True, shows the witness groups and support counts.

        Returns:
            Multi-line string suitable for printing to a terminal.
        """
        if self.tree is None:
            self.build_tree(root_witness)
        assert self.tree is not None

        # Print header
        self.console.print(
            "\n[bold cyan]╔═══════════════════════════════════════════════╗[/bold cyan]")
        self.console.print(
            "[bold cyan]║           S T E M M A   C O D I C U M          ║[/bold cyan]")
        self.console.print(
            "[bold cyan]╚═══════════════════════════════════════════════╝[/bold cyan]\n")

        # Show contamination warnings if any
        if self.contaminated_witnesses and self.contamination_overlap:
            self.console.print(
                "[bold yellow]⚠ Please examine Contamination/Anomalies[/bold yellow]")
            for overlaps in self.contamination_overlap:
                for group in overlaps:
                    # Display group info in a nice format
                    info = self.list_apparatus.group_by_frequency(
                        min_frequency=self.min_frequency).get(group)
                    if info:
                        group_str = f"  • Group: {', '.join(sorted(group))}"
                        loci_str = f"    Loci: {', '.join(info['loci'])}"
                        readings_str = f"    Readings: {'\n'.join(info['readings'])}"
                        freq_str = f"    Frequency: {info['frequency']}"
                        self.console.print(f"[yellow]{group_str}[/yellow]")
                        self.console.print(f"[dim]{loci_str}[/dim]")
                        self.console.print(f"[dim]{readings_str}[/dim]")
                        self.console.print(f"[dim]{freq_str}[/dim]")

            self.console.print(
                f"\n[red]✗ Contaminated witnesses: {', '.join(sorted(self.contaminated_witnesses))}[/red]\n")

        # Show the tree
        if show_groups:
            tree_str = self.display_groups(root_witness)
        else:
            tree_str = self._render_tree(self.tree)

        # Wrap in a panel for better visibility
        panel = Panel(
            tree_str,
            title="[bold green]Stemma[/bold green]",
            border_style="green",
            box=box.ROUNDED
        )
        self.console.print(panel)

        # Show statistics
        self._display_statistics()

        return tree_str

    def _display_statistics(self):
        """Display statistics about the stemma."""
        if self.tree is None:
            return

        # Count nodes
        def count_nodes(tree_dict):
            count = 0
            for label, children in tree_dict.items():
                count += 1
                count += count_nodes(children)
            return count

        total_nodes = count_nodes(self.tree)

        stats_table = Table(
            title="[bold]Stemma Statistics[/bold]", box=box.MINIMAL)
        stats_table.add_column("Metric", style="cyan")
        stats_table.add_column("Value", style="green")

        stats_table.add_row("Total Witnesses", str(len(self.members)))
        stats_table.add_row("Total Nodes in Tree", str(total_nodes))
        stats_table.add_row("Reference Witness", self.list_apparatus.reference)
        if self.edges:
            stats_table.add_row("Covering Relations", str(len(self.edges)))
        if self.variant_frequencies:
            stats_table.add_row("Variant Groups", str(
                len(self.variant_frequencies)))
        if self.dismissed:
            stats_table.add_row("Dismissed Groups", str(len(self.dismissed)))
            stats_table.add_row("Retained Variant Groups", str(
                len(self.variant_frequencies) - len(self.dismissed)))
        if self.contaminated_witnesses:
            stats_table.add_row("Contaminated Witnesses",
                                f"{len(self.contaminated_witnesses)} ({', '.join(sorted(self.contaminated_witnesses))})")

        self.console.print(stats_table)
        self.console.print()

    def display_variant_groups(self, minimum_group_size=1) -> None:
        """Print a beautifully formatted colour table for each variant group vs. the reference.

        For every error group the table shows the group's readings alongside
        the corresponding reference readings, rendered with Rich.
        """
        self.console.print(
            "\n[bold cyan]╔═══════════════════════════════════════════════╗[/bold cyan]")
        self.console.print(
            "[bold cyan]║       V A R I A N T   G R O U P S            ║[/bold cyan]")
        self.console.print(
            "[bold cyan]╚═══════════════════════════════════════════════╝[/bold cyan]\n")

        for idx, (group, group_variants) in enumerate(self.variant_frequencies.items(), 1):
            if len(group_variants) > minimum_group_size and len(group) > minimum_group_size:
                reference_group = self.members - group
                reference_group = reference_group.union(
                    set([self.list_apparatus.reference]))

                # Create a nice header for each group
                self.console.print(
                    f"\n[bold cyan]Group #{idx}[/bold cyan] [yellow]{', '.join(sorted(group))}[/yellow]")
                self.console.print(
                    f"[dim]vs. Reference: {', '.join(sorted(reference_group))}[/dim]")
                self.console.print(
                    f"[dim]Frequency: {group_variants['frequency']} loci[/dim]\n")

                # Create a table for readings
                table = Table(
                    box=box.ROUNDED,
                    show_header=True,
                    header_style="bold magenta"
                )
                table.add_column("#", justify="right", style="dim", width=4)
                table.add_column("Group Reading", style="cyan", no_wrap=False)
                table.add_column("Reference Reading", style="green", no_wrap=False)
                table.add_column("Locus", style="dim", no_wrap=False)

                for i, (group_reading, reference_reading, locus) in enumerate(
                        zip(group_variants["readings"],
                            group_variants["reference"],
                            group_variants["loci"]), start=1):
                    # Color code differences
                    if group_reading != reference_reading:
                        group_display = f"[bold yellow]{group_reading}[/bold yellow]"
                        ref_display = f"[bold]{reference_reading}[/bold]"
                    else:
                        group_display = f"[dim]{group_reading}[/dim]"
                        ref_display = f"[dim]{reference_reading}[/dim]"

                    table.add_row(str(i), group_display, ref_display, locus)

                self.console.print(table)

                # Add a separator between groups
                if idx < len(self.variant_frequencies):
                    self.console.print("\n[dim]─" * 80 + "[/dim]")
                self.console.print()

    def display_summary(self) -> None:
        """Display a comprehensive summary of the stemma analysis."""
        self.console.print(
            "\n[bold cyan]╔═══════════════════════════════════════════════╗[/bold cyan]")
        self.console.print(
            "[bold cyan]║      S T E M M A   S U M M A R Y             ║[/bold cyan]")
        self.console.print(
            "[bold cyan]╚═══════════════════════════════════════════════╝[/bold cyan]\n")

        # Basic info
        info_table = Table(title="[bold]Overview[/bold]", box=box.ROUNDED)
        info_table.add_column("Property", style="cyan")
        info_table.add_column("Value", style="white")

        info_table.add_row("Total Witnesses", str(len(self.members)))
        info_table.add_row("Reference Witness", self.list_apparatus.reference)
        info_table.add_row("Variant Groups", str(
            len(self.variant_frequencies)))

        if self.dismissed:
            info_table.add_row("Dismissed Groups", str(len(self.dismissed)))

        if self.edges:
            info_table.add_row("Covering Relations", str(len(self.edges)))

        if self.contaminated_witnesses:
            info_table.add_row("⚠ Contamination Detected",
                            f"[red]{len(self.contaminated_witnesses)} witnesses[/red]")
        else:
            info_table.add_row("Contamination", "[green]None detected[/green]")

        self.console.print(info_table)

        # Show variant groups count by frequency
        if self.variant_frequencies:
            self.console.print("\n[bold]Variant Groups by Frequency[/bold]")
            freq_counts = Counter(info["frequency"]
                                for info in self.variant_frequencies.values())

            freq_table = Table(box=box.MINIMAL)
            freq_table.add_column("Frequency", style="cyan", justify="right")
            freq_table.add_column("Number of Groups",
                                style="green", justify="center")

            for freq in sorted(freq_counts.keys()):
                freq_table.add_row(str(freq), str(freq_counts[freq]))

            self.console.print(freq_table)
        self.console.print()

### BACKUP WITHOUT DAG
# """Given the list apparatus, build an unoriented stemma using the group method."""
# from collections import Counter
# import re
# import warnings
# from rich.console import Console
# from rich.table import Table
# from rich.panel import Panel
# from pyfroger.list_apparatus import ListApparatus
# from rich import box


# class bcolors:
#     HEADER = '\033[95m'
#     OKBLUE = '\033[94m'
#     OKCYAN = '\033[96m'
#     OKGREEN = '\033[92m'
#     WARNING = '\033[93m'
#     FAIL = '\033[91m'
#     ENDC = '\033[0m'
#     BOLD = '\033[1m'
#     UNDERLINE = '\033[4m'

#     # Additional styles
#     INFO = '\033[94m'
#     SUCCESS = '\033[92m'
#     ERROR = '\033[91m'
#     HIGHLIGHT = '\033[103m\033[30m'  # Yellow background, black text


# class UnorientedStemma:
#     """Build an unoriented stemma using the group method
#     ('Méthode par les groupes').

#     Attributes:
#         list_apparatus: The variant groups and their different information.
#         variant_frequencies: Variant groups keyed by frequency.
#         members: Set of witness sigla available to build the tree.
#         edges: Covering-relation edges once computed, or None.
#         children_of: Direct parent→children mapping once computed, or None.
#         tree: Nested-dict stemma once built, or None.
#         contamination_overlap: Pairs of groups that partially overlap
#             (horizontal transmission), or None.
#         contaminated_witnesses: Frozenset of witnesses involved in
#             contamination, or None.
#     """

#     def __init__(self,
#                  list_apparatus: ListApparatus,
#                  min_frequency: int = 1,
#                  verbose: bool = False) -> None:
#         """Initialise the stemma builder from a variant apparatus.

#         Args:
#             list_apparatus (ListApparatus): The variant groups and their different information.
#             min_frequency (int): The minimum frequency of apparition for the group to be
#                 used for building the stemma.

#         """
#         self.min_frequency = min_frequency
#         self.list_apparatus: ListApparatus = list_apparatus
#         self.variant_frequencies =\
#             self.list_apparatus.group_by_frequency(min_frequency=min_frequency)
#         self.members: set[str] = list_apparatus.witnesses
#         self.edges: list[tuple[frozenset[str], frozenset[str]]] | None = None
#         self.children_of: dict[frozenset[str],
#                                list[frozenset[str]]] | None = None
#         self.tree: dict[str, dict] | None = None
#         self.contamination_overlap: (
#             list[tuple[frozenset[str], ...]] | None
#         ) = None
#         self.contaminated_witnesses: frozenset[str] | None = None
#         self.dismissed: set[frozenset[str]] = set()
#         self.verbose = verbose
#         # Setup Rich console
#         self.console = Console()

#     @staticmethod
#     def _compute_covering_relations(
#         groups: frozenset[frozenset[str]],
#     ) -> tuple[
#         list[tuple[frozenset[str], frozenset[str]]],
#         dict[frozenset[str], list[frozenset[str]]],
#         list[tuple[frozenset[str], ...]],
#     ]:
#         """Compute covering relations between variant groups.

#         Two groups are related by a *covering relation* when one strictly
#         contains the other with no intermediate group between them in the
#         inclusion order. Groups that partially overlap signal potential
#         horizontal transmission (contamination).

#         Args:
#             groups: All variant groups expressed as a frozenset of frozensets
#                 of witness sigla.

#         Returns:
#             A three-element tuple:

#             - **edges**: Sorted list of ``(parent, child)`` covering pairs.
#             - **children_of**: Mapping from each parent group to its list of
#               direct child groups, sorted lexicographically.
#             - **contamination_overlap**: List of ``(g1, g2)`` pairs whose
#               intersection is non-empty but neither group contains the other.
#         """
#         def covers(c: frozenset[str], p: frozenset[str]) -> bool:
#             return c < p and not any(c < u < p for u in groups)

#         contamination_overlap: set[frozenset[frozenset[str]]] = set()
#         edges: list[tuple[frozenset[str], frozenset[str]]] = []
#         children_of: dict[frozenset[str], list[frozenset[str]]] = {}

#         for c in groups:
#             for p in groups:
#                 if covers(c, p):
#                     edges.append((p, c))
#                     children_of.setdefault(p, []).append(c)
#                 else:
#                     shared = c & p
#                     if shared and shared != c and shared != p:
#                         contamination_overlap.add(frozenset((c, p)))

#         return (
#             sorted(edges),
#             {k: sorted(v, key=lambda s: sorted(s))
#              for k, v in children_of.items()},
#             [
#                 tuple(sorted(pair, key=lambda s: sorted(s)))
#                 for pair in contamination_overlap
#             ],
#         )

#     def _lattice_groups(
#         self,
#         dismissed: frozenset[frozenset[str]] = frozenset(),
#     ) -> frozenset[frozenset[str]]:
#         """Witness-set groups used to build the lattice.

#         Keeps only error groups attested at >= self.min_frequency loci
#         (problematic low-attestation groups are dropped up front), minus any
#         explicitly dismissed groups, plus the universe (all witnesses + ref).
#         """
#         universe = frozenset(self.members | {self.list_apparatus.reference})
#         kept = frozenset(
#             frozenset(g)
#             for g in self.list_apparatus.group_by_frequency(
#                 min_frequency=self.min_frequency)
#             if frozenset(g) not in dismissed
#         )
#         return kept | {universe}

#     def compute_edges(self) -> list[tuple[frozenset[str], frozenset[str]]]:
#         """Compute covering relations between error groups and store results."""
#         groups = self._lattice_groups()
#         self.edges, self.children_of, self.contamination_overlap = (
#             self._compute_covering_relations(groups)
#         )
#         self.contaminated_witnesses = frozenset(
#             w for g1, g2 in self.contamination_overlap for w in (g1 & g2)
#         )
#         return self.edges

#     def _rebuild_without(self, dismissed: set[frozenset[str]]) -> None:
#         """Recompute covering relations after excluding dismissed groups."""
#         groups = self._lattice_groups(dismissed=frozenset(dismissed))
#         self.edges, self.children_of, self.contamination_overlap = (
#             self._compute_covering_relations(groups)
#         )
#         self.contaminated_witnesses = frozenset(
#             w for g1, g2 in self.contamination_overlap for w in (g1 & g2)
#         )
#         self.dismissed = set(dismissed)
#         self.tree = None

#     @staticmethod
#     def _build_tree_structure(
#         children_of: dict[frozenset[str], list[frozenset[str]]],
#         root_group: frozenset[str],
#         root_witness: str,
#         reserved: frozenset[str] = frozenset(),
#     ) -> dict[str, dict]:
#         """Recursively build a nested-dict stemma from covering relations.

#         Each key in the returned dict is a node label (witness siglum or
#         Greek-letter hypothetical ancestor); its value is a dict of the same
#         shape representing its children.

#         In the group method each error group corresponds to a set of witnesses
#         that share a common ancestor.  Witnesses in a group that are *not*
#         subsumed by any child group sit directly under the node for that group:

#         - If exactly one such witness exists it becomes the node label (the
#           identifiable ancestor of the group at that level).
#         - If zero or more than one exist, no single witness can play that role
#           and a Greek-letter hyperarchetype (α, β, γ, …) is introduced instead.

#         Args:
#             children_of: Mapping from each parent group to its list of direct
#                 child groups.
#             root_group: The largest group, used as the starting point for
#                 recursion.
#             root_witness: Siglum placed at the root of the returned tree
#                 (in our case the selected reference witness).

#         Returns:
#             A nested dict ``{root_witness: {child_label: {…}, …}}``.
#         """
#         alphabet = "αβγδεζηθικλμνξοπρστυφχψω"
#         idx = 0

#         def next_label() -> str:
#             """Return a fresh hyperarchetype label that collides neither with
#             a witness siglum nor with a previously issued label."""
#             nonlocal idx
#             while True:
#                 letter = alphabet[idx % len(alphabet)]
#                 cycle = idx // len(alphabet)
#                 label = letter if cycle == 0 else f"{letter}{cycle + 1}"
#                 idx += 1
#                 if label not in reserved:
#                     return label

#         def build(group: frozenset[str]) -> dict[str, dict]:
#             if len(group) == 1:
#                 return {next(iter(group)): {}}

#             subs = children_of.get(group, [])
#             covered = frozenset().union(*subs) if subs else frozenset()
#             complement = group - covered
#             merged: dict[str, dict] = {}

#             if len(complement) != 1:
#                 for w in sorted(complement):
#                     merged[w] = {}

#             for sg in subs:
#                 merged.update(build(sg))

#             if len(complement) == 1:
#                 return {next(iter(complement)): merged}

#             return {next_label(): merged}

#         def _prune_dups(children: dict[str, dict]) -> set[str]:
#             """Drop leaf witnesses occurring under more than one sibling,
#             keeping the first occurrence, and report the duplicated sigla.

#             Two error groups that partially overlap are not laminar and cannot
#             both be clean clades; anomaly resolution is meant to remove such
#             overlaps upstream. If one slips through (e.g. a strong-vs-strong
#             conflict under ``raw_cutoff``), the shared witness would otherwise
#             be emitted as a leaf under *each* group — the duplicate-leaf bug.
#             We keep it once rather than fabricate a duplicate, and signal the
#             caller so the overlap can be resolved deliberately.
#             """
#             seen: set[str] = set()
#             pruned: set[str] = set()

#             def walk(node: dict[str, dict]) -> None:
#                 for label in list(node):
#                     sub = node[label]
#                     if sub:
#                         walk(sub)
#                     elif label in seen:
#                         del node[label]
#                         pruned.add(label)
#                     else:
#                         seen.add(label)

#             walk(children)
#             return pruned

#         subs = children_of.get(root_group, [])
#         covered = frozenset().union(*subs) if subs else frozenset()
#         root_children: dict[str, dict] = {
#             w: {} for w in sorted((root_group - covered) - {root_witness})
#         }

#         for sg in subs:
#             root_children.update(build(sg))

#         duplicated = _prune_dups(root_children)
#         if duplicated:
#             warnings.warn(
#                 "Witness(es) {} occur in overlapping (non-laminar) groups and "
#                 "were kept at a single position; the overlap was not resolved "
#                 "upstream. Use anomaly_method='two_by_two' or lower the cutoff "
#                 "to resolve it.".format(", ".join(sorted(duplicated))),
#                 stacklevel=2,
#             )

#         return {root_witness: root_children}

#     @staticmethod
#     def _to_newick_string(tree: dict[str, dict]) -> str:
#         """Serialize a nested-dict tree to Newick format.

#         Args:
#             tree: Nested dict as produced by :meth:`_build_tree_structure`.

#         Returns:
#             Newick string ending with ``;``, e.g. ``(A,(C,D)B)ω;``.
#         """
#         def serialize(children: dict[str, dict], label: str) -> str:
#             if not children:
#                 return label
#             inner = ",".join(
#                 serialize(grandchildren, child)
#                 for child, grandchildren in children.items()
#             )
#             return f"({inner}){label}"

#         root_label, root_children = next(iter(tree.items()))
#         return serialize(root_children, root_label) + ";"

#     def to_newick(self, root_witness: str | None = None) -> str:
#         """Return the stemma serialized in Newick format.

#         Calls :meth:`build_tree` first if it has not been called yet.

#         Args:
#             root_witness: Siglum to place at the root.  Defaults to the
#                 reference witness defined in ``list_apparatus``.

#         Returns:
#             Newick string ending with ``;``, e.g. ``(A,(C,D)B)ω;``.
#         """
#         if self.tree is None:
#             self.build_tree(root_witness)
#         assert self.tree is not None
#         return self._to_newick_string(self.tree)

#     @staticmethod
#     def _render_tree_with_groups(
#         tree: dict[str, dict],
#         witnesses: frozenset[str],
#         frequencies: dict[frozenset[str], int] | None = None,
#     ) -> str:
#         """Render a nested-dict stemma showing each internal node's group and,
#         optionally, the number of loci supporting it.

#         Each node is annotated with the witnesses appearing anywhere in its
#         subtree — an internal witness-ancestor counted in its own group. The
#         ``witnesses`` set distinguishes real sigla from Greek hyperarchetype
#         labels, which are excluded from the groups. When ``frequencies`` is
#         given, a node whose group is an error group also shows ``×N`` (its
#         locus count); nodes whose group is absent from the mapping — notably
#         the root, whose group includes the reference — show no count. Leaves
#         are left bare.

#         Args:
#             tree: Nested dict as produced by :meth:`_build_tree_structure`.
#             witnesses: Full set of witness sigla (anything else is a
#                 hyperarchetype label).
#             frequencies: Optional mapping from an error group to its number of
#                 supporting loci, e.g. ``{g: info["frequency"]}`` from
#                 ``group_by_frequency``.

#         Returns:
#             A multi-line box-drawing string with ``label  (w1, w2, …)  ×N``
#             annotations on every node that has children.
#         """
#         frequencies = frequencies or {}

#         def nat(s: str) -> list:
#             return [int(t) if t.isdigit() else t
#                     for t in re.split(r"(\d+)", s)]

#         def group_of(label: str, children: dict[str, dict]) -> set[str]:
#             g = {label} if label in witnesses else set()
#             for child, grandchildren in children.items():
#                 g |= group_of(child, grandchildren)
#             return g

#         def render(
#             label: str,
#             children: dict[str, dict],
#             prefix: str = "",
#             is_last: bool = True,
#             is_root: bool = True,
#         ) -> list[str]:
#             if children:
#                 group = group_of(label, children)
#                 sigla = ", ".join(sorted(group, key=nat))
#                 count = frequencies.get(frozenset(group))
#                 suffix = f"  ×{count}" if count is not None else ""
#                 text = f"{label}  ({sigla}){suffix}"
#             else:
#                 text = label
#             if is_root:
#                 lines, cpfx = [text], ""
#             else:
#                 lines = [prefix + ("└── " if is_last else "├── ") + text]
#                 cpfx = prefix + ("    " if is_last else "│   ")
#             items = list(children.items())
#             for i, (child, grandchildren) in enumerate(items):
#                 lines.extend(
#                     render(child, grandchildren, cpfx,
#                            i == len(items) - 1, False)
#                 )
#             return lines

#         root_label = next(iter(tree))
#         return "\n".join(render(root_label, tree[root_label]))

#     def display_groups(self, root_witness: str | None = None) -> str:
#         """Return a box-drawing stemma annotating each node with its group and
#         the number of loci supporting it.

#         Like :meth:`display_tree`, but every internal node is labelled with the
#         set of witnesses it covers and the locus count ``×N`` of the
#         corresponding error group. Calls :meth:`build_tree` first if the tree
#         has not been built yet.

#         Args:
#             root_witness: Siglum to place at the root. Defaults to the
#                 reference witness defined in ``list_apparatus``.

#         Returns:
#             Multi-line string suitable for printing to a terminal.
#         """
#         if self.tree is None:
#             self.build_tree(root_witness)
#         assert self.tree is not None
#         witnesses = frozenset(self.members | {self.list_apparatus.reference})
#         frequencies = {
#             frozenset(g): info["frequency"]
#             for g, info in self.list_apparatus.group_by_frequency(
#                 min_frequency=self.min_frequency).items()
#         }
#         return self._render_tree_with_groups(self.tree, witnesses, frequencies)

#     def build_tree(self,
#                    root_witness: str | None = None,
#                    resolve_anomaly_threshold: float = 2,
#                    anomaly_method: str = "raw_cutoff") -> dict[str, dict]:
#         """Build a nested-dict stemma from the variant apparatus.

#         Args:
#             root_witness (str): Siglum to place at the root. Defaults to the
#                 reference witness defined in ``list_apparatus``.
#             resolve_anomaly_threshold (int): Threshold below which the
#                 weaker side of an overlapping pair is treated as polygenesis
#                 and dismissed. Defaults to 2 (a single-site agreement in a
#                 conflict is dismissed); set to 1 to disable.
#             anomaly_method (str): The method selected for resolving the
#                 anomalies. Can be either:
#                     - raw_cutoff (all conflicting groups supported by fewer
#                     than `resolve_anomaly_threshold` independent sites are
#                     removed);
#                     - two_by_two: Overlaps are resolved by a two by
#                         two comparison and keeping the group with the higher
#                         number of independent sites;
#                     - probabilistic_cutoff: The resolve_anomaly threshold is
#                         understood as a probability of contamination,
#                         computed as the (number of variants supporting each
#                         conflicting group) / (total variants
#                                             in the two groups).
#                         This cut-off was described by Froger himself.

#         Returns:
#             A nested dict ``{root: {child: {…}, …}}`` representing the stemma.
#         """
#         if self.children_of is None:
#             self.compute_edges()
#         assert self.children_of is not None
#         if root_witness is None:
#             root_witness = self.list_apparatus.reference

#         if not self.children_of:
#             # No covering relations (no error groups beyond the universe):
#             # nothing to build but the root itself.
#             self.tree = {root_witness: {}}
#             return self.tree

#         if anomaly_method == "raw_cutoff":
#             self._resolve_anomaly_cutoff(
#                 min_loci=int(resolve_anomaly_threshold))
#         elif anomaly_method == "two_by_two":
#             self._resolve_anomaly_two_by_two()
#         elif anomaly_method == "probabilistic_cutoff":
#             if not 0 < resolve_anomaly_threshold < 1:
#                 raise ValueError(
#                     "anomaly_method='probabilistic_cutoff' compares a support "
#                     "ratio in ]0, 1[ against resolve_anomaly_threshold; got "
#                     f"{resolve_anomaly_threshold!r}. Pass a value between 0"
#                     "and 1"
#                 )
#             self._resolve_anomaly_probabilistic(
#                 resolve_anomaly_threshold=resolve_anomaly_threshold)
#         else:
#             raise ValueError(
#                 f"Unknown anomaly_method {anomaly_method!r}; expected one of "
#                 "'raw_cutoff', 'two_by_two', 'probabilistic_cutoff'."
#             )

#         # Recomputed *after* anomaly resolution: _rebuild_without() replaces
#         # self.children_of, so the root group must be read from the current
#         # lattice, not the pre-resolution one.
#         assert self.children_of is not None
#         if not self.children_of:
#             self.tree = {root_witness: {}}
#             return self.tree

#         root_group = max(self.children_of, key=len)

#         self.tree = self._build_tree_structure(
#             self.children_of,
#             root_group,
#             root_witness,
#             reserved=frozenset(self.members | {self.list_apparatus.reference}),
#         )

#         return self.tree

#     def _resolve_anomaly_probabilistic(
#             self,
#             resolve_anomaly_threshold: float):
#         """Resolve anomalies by computing how well each group
#         is supported as a ratio between the number of variant loci
#         and the total number of variant readings in the overlap."""
#         if self.children_of is None:
#             self.compute_edges()

#         dismissed = set()

#         if not self.contamination_overlap:
#             self.console.print(
#                 "[green]✓ No overlapping variant groups detected.[/green]")
#             return

#         if self.verbose:
#             # Display contamination info in a nice table
#             self.console.print(
#                 "\n[bold yellow]⚠ Please examine Contamination/Anomalies[/bold yellow]")

#         freqs = self.list_apparatus.group_by_frequency(
#             min_frequency=self.min_frequency)

#         for overlap in self.contamination_overlap:
#             # For each incompatible pair, count supporting sites
#             locus_counts = {
#                 group: freqs[group]["frequency"]
#                 for group in overlap
#             }

#             if self.verbose:
#                 # Create a nice display for overlaps
#                 self.console.print("\n[cyan]Overlapping groups:[/cyan]")
#                 for group in overlap:
#                     info = freqs[group]
#                     group_str = f"  • {', '.join(sorted(group))}"
#                     count_str = f" ({locus_counts[group]} loci)"
#                     loci_str = f"\n    Loci: {', '.join(info['loci'])}"
#                     readings_str = f"\n    Readings: {'\n'.join(info['readings'])}"
#                     self.console.print(f"{group_str}{count_str}")
#                     self.console.print(f"{loci_str}")
#                     self.console.print(f"{readings_str}")

#             # Check total number of variant per group
#             total = 0
#             for count in locus_counts.values():
#                 total += count

#             if not total:
#                 continue

#             # Compute probability of apparition per group
#             locus_apparition = {}
#             for group, count in locus_counts.items():
#                 locus_apparition[group] = count / total

#             for group, probability in locus_apparition.items():
#                 if probability < resolve_anomaly_threshold:
#                     # Dismiss group that are below the threshold
#                     dismissed.add(group)
#                     if self.verbose:
#                         self.console.print(
#                             "[red]✗ Dismissing group "
#                             f"{', '.join(sorted(group))} "
#                             f"(support ratio {probability:.2f} "
#                             f"< {resolve_anomaly_threshold})[/red]"
#                         )

#         # Rebuild stemma without dismissed groups
#         if dismissed:
#             self.console.print(
#                 f"\n[yellow]Rebuilding stemma after dismissing {len(dismissed)} group(s)...[/yellow]")
#             self._rebuild_without(dismissed)
#             self.console.print("[green]✓ Stemma rebuilt successfully[/green]")

#         return {
#             "dismissed": dismissed,
#         }

#     def _resolve_anomaly_two_by_two(self):
#         """Resolve anomalies by comparing the number of sites
#         in two competing groups and dismissing the weaker.
#         """
#         if self.children_of is None:
#             self.compute_edges()

#         dismissed = set()

#         if not self.contamination_overlap:
#             self.console.print(
#                 "[green]✓ No overlapping variant groups detected.[/green]")
#             return

#         if self.verbose:
#             # Display contamination info in a nice table
#             self.console.print(
#                 "\n[bold yellow]⚠"
#                 "Please examine Contamination/Anomalies[/bold yellow]")

#         freqs = self.list_apparatus.group_by_frequency(
#             min_frequency=self.min_frequency)

#         for overlap in self.contamination_overlap:
#             # For each incompatible pair, count independent supporting sites
#             # (contiguous loci collapse to a single transmission event).
#             locus_counts = {
#                 group: len(freqs[group]["loci"])
#                 for group in overlap
#             }
#             if self.verbose:
#                 # Create a nice display for this overlap
#                 self.console.print("\n[cyan]Overlapping groups:[/cyan]")
#                 for group in overlap:
#                     info = freqs[group]
#                     group_str = f"  • {', '.join(sorted(group))}"
#                     count_str = (f" ({info['frequency']} loci, "
#                                  f"{locus_counts[group]} independent)")
#                     loci_str = f"\n    Loci: {', '.join(info['loci'])}"
#                     readings_str = f"\n    Readings: {'\n'.join(info['readings'])}"
#                     self.console.print(f"{group_str}{count_str}")
#                     self.console.print(f"{loci_str}")
#                     self.console.print(f"{readings_str}")

#             weakest = min(locus_counts, key=lambda g: locus_counts[g])
#             dismissed.add(weakest)
#             if self.verbose:
#                 self.console.print(
#                     f"[red]✗ Dismissing weaker group {', '.join(sorted(weakest))} "
#                     f"({locus_counts[weakest]} independent site(s))[/red]"
#                 )

#         # Rebuild stemma without dismissed groups
#         if dismissed:
#             if self.verbose:
#                 self.console.print(
#                     f"\n[yellow]Rebuilding stemma after dismissing"
#                     f"{len(dismissed)} group(s)...[/yellow]")
#             self._rebuild_without(dismissed)
#             self.console.print("[green]✓ Stemma rebuilt successfully[/green]")

#         return {
#             "dismissed": dismissed,
#         }

#     def get_contamination_loci(self):
#         """Get all readings that led to contamination."""
#         contaminated_loci = []
#         if self.contamination_overlap:
#             for overlap in self.contamination_overlap:
#                 for group in overlap:
#                     info = self.list_apparatus.group_by_frequency(
#                         min_frequency=self.min_frequency)[group]
#                     contaminated_loci.extend(info['loci'])
#             return set(contaminated_loci)
#         else:
#             print("No contamination to check!")

#     def _resolve_anomaly_cutoff(self, min_loci: int = 2):
#         """Resolve anomalies by dismissing incompatible groups supported
#         by fewer than `min_loci` *independent* sites as probable polygenesis.

#         Args:
#             min_loci (int, default 2): Minimum number of independent variation
#                 sites (contiguous loci collapsed) required for the weaker side
#                 of a conflict to be retained as genuine contamination evidence.
#                 Set to 1 to disable (never dismiss).
#         """
#         if self.children_of is None:
#             self.compute_edges()

#         dismissed = set()

#         if not self.contamination_overlap:
#             self.console.print(
#                 "[green]✓ No overlapping variant groups detected.[/green]")
#             return

#         if self.verbose:
#             # Display contamination info in a nice table
#             self.console.print(
#                 "\n[bold yellow]⚠ Please examine Contamination/Anomalies[/bold yellow]")

#         freqs = self.list_apparatus.group_by_frequency(
#             min_frequency=self.min_frequency)

#         for overlap in self.contamination_overlap:
#             # For each incompatible pair, count supporting sites
#             locus_counts = {
#                 group: freqs[group]["frequency"]
#                 for group in overlap
#             }

#             if self.verbose:
#                 # Create a nice display for overlaps
#                 self.console.print("\n[cyan]Overlapping groups:[/cyan]")
#                 for group in overlap:
#                     info = freqs[group]
#                     group_str = f"  • {', '.join(sorted(group))}"
#                     count_str = f" ({locus_counts[group]} loci)"
#                     loci_str = f"\n    Loci: {', '.join(info['loci'])}"
#                     readings_str = f"\n    Readings: {'\n'.join(info['readings'])}"
#                     self.console.print(f"{group_str}{count_str}")
#                     self.console.print(f"{loci_str}")
#                     self.console.print(f"{readings_str}")

#             # Check for weak groups
#             for group in locus_counts:
#                 if locus_counts[group] < min_loci:
#                     if self.verbose:
#                         self.console.print(
#                             "[red]✗ Dismissing weak group "
#                             f"{', '.join(sorted(group))} "
#                             f"(only {locus_counts[group]} site(s) "
#                             f"< {min_loci})[/red]"
#                         )
#                     dismissed.add(group)

#         # Rebuild stemma without dismissed groups
#         if dismissed:
#             self.console.print(
#                 f"\n[yellow]Rebuilding stemma after dismissing {len(dismissed)} group(s)...[/yellow]")
#             self._rebuild_without(dismissed)
#             self.console.print("[green]✓ Stemma rebuilt successfully[/green]")

#         return {
#             "dismissed": dismissed,
#         }

#     @staticmethod
#     def _render_tree(tree: dict[str, dict]) -> str:
#         """Render a nested-dict tree as a Unicode ASCII-art string.

#         Args:
#             tree: Nested dict as produced by :meth:`_build_tree_structure`.

#         Returns:
#             A multi-line string with the root at the top and branches drawn
#             using ``└──`` / ``├──`` / ``│`` box-drawing characters.
#         """
#         def render(
#             label: str,
#             children: dict[str, dict],
#             prefix: str = "",
#             is_last: bool = True,
#             is_root: bool = True,
#         ) -> list[str]:
#             lines: list[str] = []
#             if is_root:
#                 lines.append(label)
#                 cpfx = ""
#             else:
#                 lines.append(prefix + ("└── " if is_last else "├── ") + label)
#                 cpfx = prefix + ("    " if is_last else "│   ")
#             items = list(children.items())
#             for i, (child, grandchildren) in enumerate(items):
#                 lines.extend(
#                     render(child, grandchildren, cpfx,
#                            i == len(items) - 1, False)
#                 )
#             return lines

#         root_label = next(iter(tree))
#         return "\n".join(render(root_label, tree[root_label]))

#     def display_tree(self, root_witness: str | None = None, show_groups: bool = False) -> str:
#         """Return a beautifully formatted Unicode ASCII-art representation of the stemma.

#         Calls :meth:`build_tree` first if it has not been called yet.

#         Args:
#             root_witness: Siglum to place at the root. Defaults to the
#                 reference witness defined in ``list_apparatus``.
#             show_groups: If True, shows the witness groups and support counts.

#         Returns:
#             Multi-line string suitable for printing to a terminal.
#         """
#         if self.tree is None:
#             self.build_tree(root_witness)
#         assert self.tree is not None

#         # Print header
#         self.console.print(
#             "\n[bold cyan]╔═══════════════════════════════════════════════╗[/bold cyan]")
#         self.console.print(
#             "[bold cyan]║           S T E M M A   C O D I C U M          ║[/bold cyan]")
#         self.console.print(
#             "[bold cyan]╚═══════════════════════════════════════════════╝[/bold cyan]\n")

#         # Show contamination warnings if any
#         if self.contaminated_witnesses and self.contamination_overlap:
#             self.console.print(
#                 "[bold yellow]⚠ Please examine Contamination/Anomalies[/bold yellow]")
#             for overlaps in self.contamination_overlap:
#                 for group in overlaps:
#                     # Display group info in a nice format
#                     info = self.list_apparatus.group_by_frequency(
#                         min_frequency=self.min_frequency).get(group)
#                     if info:
#                         group_str = f"  • Group: {', '.join(sorted(group))}"
#                         loci_str = f"    Loci: {', '.join(info['loci'])}"
#                         readings_str = f"    Readings: {'\n'.join(info['readings'])}"
#                         freq_str = f"    Frequency: {info['frequency']}"
#                         self.console.print(f"[yellow]{group_str}[/yellow]")
#                         self.console.print(f"[dim]{loci_str}[/dim]")
#                         self.console.print(f"[dim]{readings_str}[/dim]")
#                         self.console.print(f"[dim]{freq_str}[/dim]")

#             self.console.print(
#                 f"\n[red]✗ Contaminated witnesses: {', '.join(sorted(self.contaminated_witnesses))}[/red]\n")

#         # Show the tree
#         if show_groups:
#             tree_str = self.display_groups(root_witness)
#         else:
#             tree_str = self._render_tree(self.tree)

#         # Wrap in a panel for better visibility
#         panel = Panel(
#             tree_str,
#             title="[bold green]Stemma[/bold green]",
#             border_style="green",
#             box=box.ROUNDED
#         )
#         self.console.print(panel)

#         # Show statistics
#         self._display_statistics()

#         return tree_str

#     def _display_statistics(self):
#         """Display statistics about the stemma."""
#         if self.tree is None:
#             return

#         # Count nodes
#         def count_nodes(tree_dict):
#             count = 0
#             for label, children in tree_dict.items():
#                 count += 1
#                 count += count_nodes(children)
#             return count

#         total_nodes = count_nodes(self.tree)

#         stats_table = Table(
#             title="[bold]Stemma Statistics[/bold]", box=box.MINIMAL)
#         stats_table.add_column("Metric", style="cyan")
#         stats_table.add_column("Value", style="green")

#         stats_table.add_row("Total Witnesses", str(len(self.members)))
#         stats_table.add_row("Total Nodes in Tree", str(total_nodes))
#         stats_table.add_row("Reference Witness", self.list_apparatus.reference)
#         if self.edges:
#             stats_table.add_row("Covering Relations", str(len(self.edges)))
#         if self.variant_frequencies:
#             stats_table.add_row("Variant Groups", str(
#                 len(self.variant_frequencies)))
#         if self.dismissed:
#             stats_table.add_row("Dismissed Groups", str(len(self.dismissed)))
#             stats_table.add_row("Retained Variant Groups", str(
#                 len(self.variant_frequencies) - len(self.dismissed)))
#         if self.contaminated_witnesses:
#             stats_table.add_row("Contaminated Witnesses",
#                                 f"{len(self.contaminated_witnesses)} ({', '.join(sorted(self.contaminated_witnesses))})")

#         self.console.print(stats_table)
#         self.console.print()

#     def display_variant_groups(self) -> None:
#         """Print a beautifully formatted colour table for each variant group vs. the reference.

#         For every error group the table shows the group's readings alongside
#         the corresponding reference readings, rendered with Rich.
#         """
#         self.console.print(
#             "\n[bold cyan]╔═══════════════════════════════════════════════╗[/bold cyan]")
#         self.console.print(
#             "[bold cyan]║       V A R I A N T   G R O U P S            ║[/bold cyan]")
#         self.console.print(
#             "[bold cyan]╚═══════════════════════════════════════════════╝[/bold cyan]\n")

#         for idx, (group, group_variants) in enumerate(self.variant_frequencies.items(), 1):
#             reference_group = self.members - group
#             reference_group = reference_group.union(
#                 set([self.list_apparatus.reference]))

#             # Create a nice header for each group
#             self.console.print(
#                 f"\n[bold cyan]Group #{idx}[/bold cyan] [yellow]{', '.join(sorted(group))}[/yellow]")
#             self.console.print(
#                 f"[dim]vs. Reference: {', '.join(sorted(reference_group))}[/dim]")
#             self.console.print(
#                 f"[dim]Frequency: {group_variants['frequency']} loci[/dim]\n")

#             # Create a table for readings
#             table = Table(
#                 box=box.ROUNDED,
#                 show_header=True,
#                 header_style="bold magenta"
#             )
#             table.add_column("#", justify="right", style="dim", width=4)
#             table.add_column("Group Reading", style="cyan", no_wrap=False)
#             table.add_column("Reference Reading", style="green", no_wrap=False)
#             table.add_column("Locus", style="dim", no_wrap=False)

#             for i, (group_reading, reference_reading, locus) in enumerate(
#                     zip(group_variants["readings"],
#                         group_variants["reference"],
#                         group_variants["loci"]), start=1):
#                 # Color code differences
#                 if group_reading != reference_reading:
#                     group_display = f"[bold yellow]{group_reading}[/bold yellow]"
#                     ref_display = f"[bold]{reference_reading}[/bold]"
#                 else:
#                     group_display = f"[dim]{group_reading}[/dim]"
#                     ref_display = f"[dim]{reference_reading}[/dim]"

#                 table.add_row(str(i), group_display, ref_display, locus)

#             self.console.print(table)

#             # Add a separator between groups
#             if idx < len(self.variant_frequencies):
#                 self.console.print("\n[dim]─" * 80 + "[/dim]")

#         self.console.print()

#     def display_summary(self) -> None:
#         """Display a comprehensive summary of the stemma analysis."""
#         self.console.print(
#             "\n[bold cyan]╔═══════════════════════════════════════════════╗[/bold cyan]")
#         self.console.print(
#             "[bold cyan]║      S T E M M A   S U M M A R Y             ║[/bold cyan]")
#         self.console.print(
#             "[bold cyan]╚═══════════════════════════════════════════════╝[/bold cyan]\n")

#         # Basic info
#         info_table = Table(title="[bold]Overview[/bold]", box=box.ROUNDED)
#         info_table.add_column("Property", style="cyan")
#         info_table.add_column("Value", style="white")

#         info_table.add_row("Total Witnesses", str(len(self.members)))
#         info_table.add_row("Reference Witness", self.list_apparatus.reference)
#         info_table.add_row("Variant Groups", str(
#             len(self.variant_frequencies)))

#         if self.dismissed:
#             info_table.add_row("Dismissed Groups", str(len(self.dismissed)))

#         if self.edges:
#             info_table.add_row("Covering Relations", str(len(self.edges)))

#         if self.contaminated_witnesses:
#             info_table.add_row("⚠ Contamination Detected",
#                                f"[red]{len(self.contaminated_witnesses)} witnesses[/red]")
#         else:
#             info_table.add_row("Contamination", "[green]None detected[/green]")

#         self.console.print(info_table)

#         # Show variant groups count by frequency
#         if self.variant_frequencies:
#             self.console.print("\n[bold]Variant Groups by Frequency[/bold]")
#             freq_counts = Counter(info["frequency"]
#                                   for info in self.variant_frequencies.values())

#             freq_table = Table(box=box.MINIMAL)
#             freq_table.add_column("Frequency", style="cyan", justify="right")
#             freq_table.add_column("Number of Groups",
#                                   style="green", justify="center")

#             for freq in sorted(freq_counts.keys()):
#                 freq_table.add_row(str(freq), str(freq_counts[freq]))

#             self.console.print(freq_table)

#         self.console.print()
