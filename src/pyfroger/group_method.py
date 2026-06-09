"""Given the list apparatus, build an unoriented stemma using the group method."""
from rich.console import Console
from rich.table import Table
from pyfroger.list_apparatus import ListApparatus


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

    def __init__(self, list_apparatus: ListApparatus) -> None:
        """Initialise the stemma builder from a variant apparatus.

        Args:
            list_apparatus: The variant groups and their different information.
        """
        self.list_apparatus: ListApparatus = list_apparatus
        self.variant_frequencies = self.list_apparatus.group_by_frequency()
        self.members: set[str] = list_apparatus.witnesses
        self.edges: list[tuple[frozenset[str], frozenset[str]]] | None = None
        self.children_of: dict[frozenset[str], list[frozenset[str]]] | None = None
        self.tree: dict[str, dict] | None = None
        self.contamination_overlap: (
            list[tuple[frozenset[str], ...]] | None
        ) = None
        self.contaminated_witnesses: frozenset[str] | None = None

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

    def compute_edges(self) -> list[tuple[frozenset[str], frozenset[str]]]:
        """Compute covering relations between error groups and store results.

        Populates ``self.edges``, ``self.children_of``,
        ``self.contamination_overlap``, and ``self.contaminated_witnesses``.

        Returns:
            Sorted list of ``(parent, child)`` covering-relation edges.
        """
        groups = frozenset(
            [frozenset(g.witnesses)
             for g in self.list_apparatus.error_groups()]
            # + [frozenset(self.list_apparatus.witnesses)]
            + [frozenset(self.list_apparatus.witnesses |
                         {self.list_apparatus.reference})]
        )

        self.edges, self.children_of, self.contamination_overlap = (
            self._compute_covering_relations(groups)
        )
        self.contaminated_witnesses = frozenset(
            w for g1, g2 in self.contamination_overlap for w in (g1 & g2)
        )

        return self.edges

    @staticmethod
    def _build_tree_structure(
        children_of: dict[frozenset[str], list[frozenset[str]]],
        root_group: frozenset[str],
        root_witness: str,
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
        hyp = list("αβγδεζηθικλμνξοπρστ")
        idx = 0

        def build(group: frozenset[str]) -> dict[str, dict]:
            nonlocal idx

            if len(group) == 1:
                return {next(iter(group)): {}}

            subs = children_of.get(group, [])
            covered = frozenset().union(*subs) if subs else frozenset()
            complement = group - covered
            merged: dict[str, dict] = {}

            if len(complement) != 1:
                for w in sorted(complement):
                    merged[w] = {}

            for sg in subs:
                merged.update(build(sg))

            if len(complement) == 1:
                return {next(iter(complement)): merged}

            label = hyp[idx]
            idx += 1
            return {label: merged}

        subs = children_of.get(root_group, [])
        covered = frozenset().union(*subs) if subs else frozenset()
        root_children: dict[str, dict] = {
            w: {} for w in sorted((root_group - covered) - {root_witness})
        }

        for sg in subs:
            root_children.update(build(sg))

        return {root_witness: root_children}

    def build_tree(self, root_witness: str | None = None) -> dict[str, dict]:
        """Build a nested-dict stemma from the variant apparatus.

        Calls :meth:`compute_edges` first if it has not been called yet.

        Args:
            root_witness: Siglum to place at the root. Defaults to the
                reference witness defined in ``list_apparatus``.

        Returns:
            A nested dict ``{root: {child: {…}, …}}`` representing the stemma.
        """
        if self.children_of is None:
            self.compute_edges()
        assert self.children_of is not None
        if root_witness is None:
            root_witness = self.list_apparatus.reference

        root_group = max(self.children_of, key=len)
        self.tree = self._build_tree_structure(
            self.children_of, root_group, root_witness)
        return self.tree

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
                    render(child, grandchildren, cpfx, i == len(items) - 1, False)
                )
            return lines

        root_label = next(iter(tree))
        return "\n".join(render(root_label, tree[root_label]))

    def display_tree(self, root_witness: str | None = None) -> str:
        """Return a Unicode ASCII-art representation of the stemma.

        Calls :meth:`build_tree` first if it has not been called yet.

        Args:
            root_witness: Siglum to place at the root.  Defaults to the
                reference witness defined in ``list_apparatus``.

        Returns:
            Multi-line string suitable for printing to a terminal.
        """
        if self.tree is None:
            self.build_tree(root_witness)
        assert self.tree is not None
        return self._render_tree(self.tree)

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

    def detect_contamination(self) -> None:
        """Print a warning listing witnesses susceptible to contamination.

        Calls :meth:`compute_edges` first if it has not been called yet.
        Does nothing when no contamination is detected.
        """
        if self.children_of is None:
            self.compute_edges()
        if self.contaminated_witnesses:
            print(f"! Witnesses susceptible to contamination:\
                  {' '.join(self.contaminated_witnesses)}")

    def display_variant_groups(self) -> None:
        """Print a colour table for each variant group vs. the reference.

        For every error group the table shows the group's readings alongside
        the corresponding reference readings, rendered with Rich.
        """
        console = Console()

        for group, group_variants in self.variant_frequencies.items():
            reference_group = self.members - group

            table = Table(
                title=(
                    f"[bold cyan]Groupe[/] ({', '.join(sorted(group))})"
                    f"   vs   "
                    f"[bold magenta]Référence[/] ({', '.join(sorted(reference_group))})"
                ),
                title_justify="left",
                header_style="bold",
            )
            table.add_column("#", justify="right", style="dim")
            table.add_column("Leçon du groupe", style="cyan")
            table.add_column("Leçon de référence", style="magenta")

            for i, (group_reading, reference_reading) in enumerate(
                    zip(group_variants["readings"],
                        group_variants["reference"]), start=1):
                table.add_row(str(i), group_reading, reference_reading)

            console.print(table)
            console.print()
