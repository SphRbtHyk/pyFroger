# """Implementation of the Froger's group method."""
# import sys


# class StemmaPerGroup:
#     """Implementation of the group method as defined by Jacques Froger.

#     The algorithm workflow is as follows:
#         - Order the manuscripts by the size of their group.
#         - Compute the recursive inclusion of each group within another group.
#         - Compute anomalies.
#         - Display the inclusion as a tree.
#     """

#     def __init__(self,
#                  variant_group) -> None:
#         """Initialize the group method.

#         Args:
#             variant_group (ListApparatus): A ListApparatus object
#                 that contains the relationshup between the variant readings
#                 and the manuscripts bearing them.
#         """
#         self.variant_group = variant_group
#         self.variant_frequencies = {k: v["frequency"] for k, v in
#                                     variant_group.group_by_frequency().items()}

#     def compute_levels(self) -> dict[str, list[set[str]]]:
#         """Order the manuscripts by their respective levels in the tree,
#         i.e. the number of manuscripts that compose them.
#         See table p. 253."""
#         levels_ = {str(len(self.variant_group.witnesses)):
#                    [set(self.variant_group.witnesses)]}

#         for witnesses in self.variant_frequencies.keys():
#             if str(len(witnesses)) in levels_:
#                 levels_[str(len(witnesses))].append(set(witnesses))
#             else:
#                 levels_[str(len(witnesses))] = [set(witnesses)]
#         return levels_

#     def compute_group_inclusion(self):
#         """Find out the relationship between two group of manuscripts
#         between two adjacent levels.
#         Outputs a structure such as:
#         {level_1: [{child: parent}, {child: parent}], level_2: ...}

#         ! Question de la gestion des anomalies? 
#         Question du display des variants partagés?
#         """
#         tree_structure = {}
#         levels = self.compute_levels()
#         sorted_levels = sorted([int(level) for level in levels.keys()])
#         for level in sorted_levels[:-1]:  # Last level is necessarily all MS
#             tree_structure[str(level)] = []
#             # For each manuscript of the level
#             # Check until inclusion and store as ancestor
#             for ms_group in levels[str(level)]:

#                 for ix in range(level+1, max(sorted_levels)+1):
#                     if str(ix) in levels:
#                         for candidate_group in levels[str(ix)]:
#                             if ms_group.issubset(candidate_group):
#                                 tree_structure[str(level)].append(
#                                     [{"child": ms_group,
#                                       "parent": candidate_group}])

#         return tree_structure

#     def display_group_inclusion(self, color=None):
#         """Display the inclusion of the groups from the structure."""
#         tree_structure = self.compute_group_inclusion()
#         if color is None:
#             color = sys.stdout.isatty()
#         paint = (lambda c, s: f"\033[{c}m{s}\033[0m") if color else (lambda c, s: s)
#         HEAD = lambda s: paint("1;36", s)     # bold cyan
#         DIM = lambda s: paint("2", s)         # dim
#         RED = lambda s: paint("1;31", s)      # bold red
#         ORANGE = lambda s: paint("33", s)     # yellow/orange
#         GREEN = lambda s: paint("1;32", s)    # bold green
#         SUB = "\u2282"
    
#         # 1. every distinct group
#         groups = set()
#         for level in tree_structure.values():
#             for alternatives in level:
#                 for pair in alternatives:
#                     groups.add(frozenset(pair["child"]))
#                     groups.add(frozenset(pair["parent"]))
    
#         # 2. covering relation (transitive reduction)
#         def covers(c, p):
#             return c < p and not any(c < u < p for u in groups)
    
#         immediate = {}
#         for c in groups:
#             for p in groups:
#                 if covers(c, p):
#                     immediate.setdefault(c, set()).add(p)
    
#         # 3. contamination
#         contaminated_nodes, contaminated_ms = set(), set()
#         for c, ps in immediate.items():
#             ps = sorted(ps, key=lambda s: sorted(s))
#             for i in range(len(ps)):
#                 for j in range(i + 1, len(ps)):
#                     a, b = ps[i], ps[j]
#                     if not (a <= b or b <= a):
#                         contaminated_nodes.add(c)
#                         contaminated_ms |= (a & b)
    
#         def fmt(s, hot=False):
#             if color and hot:
#                 inner = ", ".join(RED(m) if m in contaminated_ms else m
#                                 for m in sorted(s))
#             else:
#                 inner = ", ".join(sorted(s))
#             return "{" + inner + "}"
    
#         # 4. print by level = child cardinality
#         by_level = {}
#         for c, ps in immediate.items():
#             for p in ps:
#                 by_level.setdefault(len(c), []).append((c, p))
    
#         for level in sorted(by_level):
#             print(HEAD(f"── Level {level} ──"))
#             pairs = sorted(by_level[level],
#                         key=lambda cp: (sorted(cp[0]), sorted(cp[1])))
#             for i, (c, p) in enumerate(pairs, 1):
#                 hot = c in contaminated_nodes
#                 line = f"  [{i}]  {fmt(c, hot)}  {SUB}  {fmt(p, hot)}"
#                 if hot:
#                     line += "   " + RED("\u26a0 contamination")
#                 print(line)
#             print()
    
#         # 5. Print summary
#         if contaminated_ms:
#             print(RED("\u26a0 Contamination detected"))
#             for c in sorted(contaminated_nodes, key=lambda s: sorted(s)):
#                 ps = sorted(immediate[c], key=lambda s: sorted(s))
#                 joined = ("  " + ORANGE("\u27f7") + "  ").join(ORANGE(fmt(p)) for p in ps)
#                 print("   " + fmt(c, True)
#                     + DIM(" descends from incompatible groups: ") + joined)
#             print("   " + DIM("witness(es) to inspect: ")
#                 + RED(", ".join(sorted(contaminated_ms))))
#         else:
#             print(GREEN("\u2713 No contamination \u2014 clean tree."))

#     def resolve_witnesses(self):
#         """Position the different witnesses in the tree given the inclusion."""

#     def detect_anomalies(self):
#         """Detect contamination and flag out the different readings
#         that led to this anormality."""

#     def reconstruct_archetype(self) -> None:
#         """Reconstruct archetype?"""
"""Implementation of the Froger's group method."""
import sys
from dataclasses import dataclass, field
from collections import defaultdict
from list_apparatus import ListApparatus


def _group_label(s) -> str:
    """Canonical, human-readable label for a group of witnesses."""
    return "{" + ", ".join(sorted(s)) + "}"


@dataclass
class Stemma:
    """Derived stemma structure, the single source of truth for display,
    anomaly detection and (later) archetype reconstruction."""
    members: dict          # label -> frozenset of witnesses
    edges: list            # (parent_label, child_label), covering relation
    parents: dict          # child_label  -> set(parent_labels)
    children: dict         # parent_label -> set(child_labels)
    depth: dict            # label -> int (0 = root/archetype, grows downward)
    roots: list            # labels with no parent
    convergence: set       # labels with >1 parent (structural contamination)
    conflicts: list        # (labelA, labelB, shared frozenset) incompatible pairs
    contaminated: set      # witnesses to inspect
    positioned: dict = field(default_factory=dict)  # label -> witnesses sitting at node
    lost: set = field(default_factory=set)          # labels = lost/hypothetical MS
    support: dict = field(default_factory=dict)     # label -> agreements (loci count)


class StemmaPerGroup:
    """Implementation of the group method as defined by Jacques Froger.

    The algorithm workflow is as follows:
        - Order the manuscripts by the size of their group.
        - Compute the recursive inclusion of each group within another group.
        - Compute anomalies.
        - Display the inclusion as a tree.
    """

    def __init__(self,
                 variant_group) -> None:
        """Initialize the group method.

        Args:
            variant_group (ListApparatus): A ListApparatus object
                that contains the relationshup between the variant readings
                and the manuscripts bearing them.
        """
        self.variant_group = variant_group
        self.variant_frequencies = {k: v["frequency"] for k, v in
                                    variant_group.group_by_frequency().items()}
        self._stemma = None     # lazy cache for build_stemma()

    def compute_levels(self) -> dict[str, list[set[str]]]:
        """Order the manuscripts by their respective levels in the tree,
        i.e. the number of manuscripts that compose them.
        See table p. 253."""
        levels_ = {str(len(self.variant_group.witnesses)):
                   [set(self.variant_group.witnesses)]}
        levels_["1"] = [set([wit]) for wit in self.variant_group.witnesses]
        for witnesses in self.variant_frequencies.keys():
            if str(len(witnesses)) in levels_:
                levels_[str(len(witnesses))].append(set(witnesses))
            else:
                levels_[str(len(witnesses))] = [set(witnesses)]
        return levels_

    def compute_group_inclusion(self):
        """Find out the relationship between two group of manuscripts
        between two adjacent levels.
        Outputs a structure such as:
        {level_1: [{child: parent}, {child: parent}], level_2: ...}
        """
        tree_structure = {}
        levels = self.compute_levels()
        sorted_levels = sorted([int(level) for level in levels.keys()])
        for level in sorted_levels[:-1]:  # Last level is necessarily all MS
            tree_structure[str(level)] = []
            for ms_group in levels[str(level)]:
                for ix in range(level+1, max(sorted_levels)+1):
                    if str(ix) in levels:
                        for candidate_group in levels[str(ix)]:
                            if ms_group.issubset(candidate_group):
                                tree_structure[str(level)].append(
                                    [{"child": ms_group,
                                      "parent": candidate_group}])
        return tree_structure

    def build_stemma(self, rebuild: bool = False) -> Stemma:
        """Build the stemma from the groups: covering relation (transitive
        reduction), descent depth, and contamination. Structure is DERIVED
        from the distinct groups, so it is correct whether the inclusion data
        lists immediate parents or every ancestor.
        """
        if self._stemma is not None and not rebuild:
            return self._stemma

        # 1. distinct groups = every set produced by compute_levels()
        groups = {frozenset(g)
                  for sets in self.compute_levels().values() for g in sets}
        members = {_group_label(g): g for g in groups}

        # 2. covering relation = transitive reduction
        def covers(c, p):
            return c < p and not any(c < u < p for u in groups)

        edges = []
        parents, children = defaultdict(set), defaultdict(set)
        for c in groups:
            for p in groups:
                if covers(c, p):
                    lc, lp = _group_label(c), _group_label(p)
                    edges.append((lp, lc))
                    parents[lc].add(lp)
                    children[lp].add(lc)
        edges.sort()

        # 3. depth = longest path from a root (memoised)
        depth = {}

        def dep(lab):
            if lab not in depth:
                ps = parents.get(lab)
                depth[lab] = 0 if not ps else 1 + max(dep(x) for x in ps)
            return depth[lab]

        for lab in members:
            dep(lab)
        roots = sorted(lab for lab in members if not parents.get(lab))

        # 4. structural contamination: a node with more than one parent
        #    (in a transitive reduction these parents are always incompatible)
        convergence = {lab for lab in members if len(parents.get(lab, ())) > 1}

        # 5. witness-level contamination: compatibility criterion
        labs = list(members)
        conflicts, contaminated = [], set()
        for i in range(len(labs)):
            for j in range(i + 1, len(labs)):
                A, B = members[labs[i]], members[labs[j]]
                inter = A & B
                if inter and not (A <= B or B <= A):
                    conflicts.append((labs[i], labs[j], inter))
                    contaminated |= inter

        # 6. position witnesses: a node's own witness(es) are what is left
        #    after subtracting every child group. Empty => lost/hypothetical MS.
        positioned, lost = {}, set()
        for lab in members:
            covered = set()
            for ch in children.get(lab, ()):
                covered |= members[ch]
            remainder = members[lab] - covered
            positioned[lab] = frozenset(remainder)
            if not remainder:
                lost.add(lab)

        # 7. agreement support: how many loci back each group (agreements
        #    "against the others"). This is the statistic that tells strong
        #    contamination from a coincidental (polygenetic) single-locus clash.
        support = {}
        for g, freq in self.variant_frequencies.items():
            support[_group_label(g)] = freq

        self._stemma = Stemma(members=members, edges=edges,
                              parents=dict(parents), children=dict(children),
                              depth=depth, roots=roots, convergence=convergence,
                              conflicts=conflicts, contaminated=contaminated,
                              positioned=positioned, lost=lost, support=support)
        return self._stemma

    def display_group_inclusion(self, color=None):
        """Display the simplified inclusion tree and flag contamination."""
        stemma = self.build_stemma()
        if color is None:
            color = sys.stdout.isatty()
        paint = (lambda c, s: f"\033[{c}m{s}\033[0m") if color else (lambda c, s: s)
        HEAD = lambda s: paint("1;36", s)     # bold cyan
        DIM = lambda s: paint("2", s)         # dim
        RED = lambda s: paint("1;31", s)      # bold red
        ORANGE = lambda s: paint("33", s)     # yellow/orange
        GREEN = lambda s: paint("1;32", s)    # bold green
        SUB = "\u2282"

        def fmt(lbl, hot=False):
            s = stemma.members[lbl]
            if color and hot:
                inner = ", ".join(RED(m) if m in stemma.contaminated else m
                                  for m in sorted(s))
            else:
                inner = ", ".join(sorted(s))
            return "{" + inner + "}"

        # group covering edges by level = child cardinality
        by_level = {}
        for parent_lbl, child_lbl in stemma.edges:
            by_level.setdefault(len(stemma.members[child_lbl]), []).append(
                (child_lbl, parent_lbl))

        for level in sorted(by_level):
            print(HEAD(f"── Level {level} ──"))
            pairs = sorted(by_level[level],
                           key=lambda cp: (sorted(stemma.members[cp[0]]),
                                           sorted(stemma.members[cp[1]])))
            for i, (c, p) in enumerate(pairs, 1):
                hot = c in stemma.convergence
                line = f"  [{i}]  {fmt(c, hot)}  {SUB}  {fmt(p, hot)}"
                if hot:
                    line += "   " + RED("\u26a0 contamination")
                print(line)
            print()

        # summary
        if stemma.contaminated:
            print(RED("\u26a0 Contamination detected"))
            for c in sorted(stemma.convergence,
                            key=lambda lbl: sorted(stemma.members[lbl])):
                ps = sorted(stemma.parents[c],
                            key=lambda lbl: sorted(stemma.members[lbl]))
                joined = ("  " + ORANGE("\u27f7") + "  ").join(
                    ORANGE(fmt(p)) for p in ps)
                print("   " + fmt(c, True)
                      + DIM(" descends from incompatible groups: ") + joined)
            print("   " + DIM("witness(es) to inspect: ")
                  + RED(", ".join(sorted(stemma.contaminated))))
        else:
            print(GREEN("\u2713 No contamination \u2014 clean tree."))

    def resolve_witnesses(self, verbose: bool = True, color=None):
        """Position the witnesses in the tree.

        Each node's own witness(es) = its group minus the union of its child
        groups (the readings that single it out). If that difference is empty
        the node is a LOST / hypothetical manuscript -- an inferred ancestor
        with no surviving witness sitting exactly at that place.

        Returns {label: {"witnesses": [...], "lost": bool}} and, if verbose,
        prints the stemma as an indented descent tree.
        """
        stemma = self.build_stemma()
        result = {lab: {"witnesses": sorted(stemma.positioned[lab]),
                        "lost": lab in stemma.lost}
                  for lab in stemma.members}

        if verbose:
            if color is None:
                color = sys.stdout.isatty()
            paint = (lambda c, s: f"\033[{c}m{s}\033[0m") if color \
                else (lambda c, s: s)
            DIM = lambda s: paint("2", s)
            RED = lambda s: paint("1;31", s)
            GREEN = lambda s: paint("1;32", s)
            MAG = lambda s: paint("35", s)

            def show(lab, indent):
                pad = "  " * indent
                if lab in stemma.lost:
                    kind = "archetype" if lab in stemma.roots else "hypothetical"
                    tag = MAG(f"\u2736 lost ({kind})")
                else:
                    tag = GREEN("\u2192 " + ", ".join(sorted(stemma.positioned[lab])))
                warn = "  " + RED("\u26a0 contamination") \
                    if lab in stemma.convergence else ""
                print(f"{pad}{lab}   {tag}{warn}")
                for ch in sorted(stemma.children.get(lab, ()),
                                 key=lambda l: sorted(stemma.members[l])):
                    show(ch, indent + 1)

            for r in stemma.roots:
                show(r, 0)

        return result

    def detect_anomalies(self):
        """Detect contamination and flag out the readings that led to it.

        Returns a dict the rest of the pipeline can act on programmatically.
        """
        stemma = self.build_stemma()
        return {
            "contaminated_witnesses": sorted(stemma.contaminated),
            "convergence_nodes": sorted(stemma.convergence),
            "conflicts": [(a, b, sorted(shared))
                          for a, b, shared in stemma.conflicts],
        }

    def _readings_by_group(self):
        """Best-effort {frozenset(witnesses): [reading, ...]} from the
        apparatus. Looks for a readings field in group_by_frequency()'s value."""
        out = {}
        for group, info in self.variant_group.group_by_frequency().items():
            rs = []
            if isinstance(info, dict):
                for key in ("readings", "reading", "variants", "variant", "lemmas"):
                    if key in info:
                        rs = info[key]
                        break
            if not isinstance(rs, (list, tuple, set)):
                rs = [rs]
            out[frozenset(group)] = list(rs)
        return out

    def plot_stemma(self, readings=None, save_path=None, figsize=(12.5, 8)):
        """Draw the stemma with matplotlib.

        Shows the covering edges, the shared readings at each node, and the
        contamination / lost-manuscript highlights. Pure matplotlib.

        readings : optional {frozenset(witnesses): [reading, ...]}. If omitted,
                   pulled from variant_group via _readings_by_group().
        save_path: save a PNG there; if None, plt.show().
        """
        import random
        import matplotlib.pyplot as plt
        from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

        stemma = self.build_stemma()
        if readings is None:
            readings = self._readings_by_group()
        reading_label = {_group_label(g): ", ".join(map(str, rs))
                         for g, rs in readings.items() if rs}
        conflict_nodes = {g for a, b, _ in stemma.conflicts for g in (a, b)}

        def support_label(lab):
            """Locus-level annotation: number of agreements (loci) backing the
            group; the bare reading value is not shown, only how well-supported
            the grouping is."""
            n = stemma.support.get(lab)
            if not n:
                return None
            return f"{n} loci" if n != 1 else "1 locus"

        FLAG_RED, WARN, GREY = '#c0392b', '#e67e22', '#95a5a6'
        GAP_X, GAP_Y = 3.4, 1.7
        edges = stemma.edges
        depth = stemma.depth
        maxd = max(depth.values()) if depth else 0
        tiers = {}
        for n, d in depth.items():
            tiers.setdefault(d, []).append(n)
        layer_keys = sorted(tiers)

        def y_of(d):
            return (maxd - d) * GAP_Y

        def coords(order):
            pos = {}
            for d in layer_keys:
                row = order[d]
                for i, n in enumerate(row):
                    pos[n] = ((i - (len(row) - 1) / 2) * GAP_X, y_of(d))
            return pos

        def cost(order):
            pos = coords(order)

            def ccw(a, b, c):
                return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])

            cross = 0
            for i in range(len(edges)):
                for j in range(i + 1, len(edges)):
                    a, b = edges[i]
                    c, d = edges[j]
                    if len({a, b, c, d}) < 4:
                        continue
                    p1, p2, p3, p4 = pos[a], pos[b], pos[c], pos[d]
                    if (ccw(p3, p4, p1) > 0) != (ccw(p3, p4, p2) > 0) and \
                       (ccw(p1, p2, p3) > 0) != (ccw(p1, p2, p4) > 0):
                        cross += 1
            through = 0
            for a, b in edges:
                (xa, ya), (xb, yb) = pos[a], pos[b]
                lo, hi = sorted((ya, yb))
                for d in layer_keys:
                    yl = y_of(d)
                    if lo < yl < hi:
                        ex = xa + (xb - xa) * ((yl - ya) / (yb - ya))
                        for m in tiers[d]:
                            if m not in (a, b) and abs(ex - pos[m][0]) < GAP_X * 0.45:
                                through += 1
            return cross * 10 + through

        best, bcost = None, float('inf')
        for seed in range(40):
            rng = random.Random(seed)
            order = {d: list(tiers[d]) for d in layer_keys}
            for d in layer_keys:
                rng.shuffle(order[d])
            improved = True
            while improved:
                improved = False
                for d in layer_keys:
                    row = order[d]
                    for i in range(len(row) - 1):
                        base = cost(order)
                        row[i], row[i + 1] = row[i + 1], row[i]
                        if cost(order) < base:
                            improved = True
                        else:
                            row[i], row[i + 1] = row[i + 1], row[i]
            c = cost(order)
            if c < bcost:
                bcost, best = c, {d: list(order[d]) for d in layer_keys}

        nbrs = {}
        for p, c in edges:
            nbrs.setdefault(p, set()).add(c)
            nbrs.setdefault(c, set()).add(p)
        x = {}
        for d in layer_keys:
            row = best[d]
            for i, n in enumerate(row):
                x[n] = (i - (len(row) - 1) / 2) * GAP_X
        for _ in range(150):
            for d in layer_keys:
                row = best[d]
                for n in row:
                    nb = nbrs.get(n)
                    if nb:
                        x[n] = sum(x[m] for m in nb) / len(nb)
                for i in range(1, len(row)):
                    ov = (x[row[i - 1]] + GAP_X) - x[row[i]]
                    if ov > 0:
                        x[row[i - 1]] -= ov / 2
                        x[row[i]] += ov / 2
        mx = sum(x.values()) / len(x)
        pos = {n: (x[n] - mx, y_of(depth[n])) for n in stemma.members}

        xs = [p[0] for p in pos.values()]
        ys = [p[1] for p in pos.values()]
        fig, ax = plt.subplots(figsize=figsize)
        ax.set_xlim(min(xs) - 2.6, max(xs) + 2.6)
        ax.set_ylim(min(ys) - 1.1, max(ys) + 1.0)
        ax.axis('off')
        fig.canvas.draw()
        inv = ax.transData.inverted()
        rend = fig.canvas.get_renderer()

        def wh(t):
            bb = t.get_window_extent(renderer=rend)
            (x0, y0), (x1, y1) = inv.transform([(bb.x0, bb.y0), (bb.x1, bb.y1)])
            return abs(x1 - x0), abs(y1 - y0)

        def runs(lab):
            mem = sorted(stemma.members[lab])
            r = [("{", 'black')]
            for k, m in enumerate(mem):
                r.append((m, FLAG_RED if m in stemma.contaminated else 'black'))
                if k < len(mem) - 1:
                    r.append((", ", 'black'))
            r.append(("}", 'black'))
            return r

        boxes = {}
        for n, (px, py) in pos.items():
            rtexts, rws, lab_h = [], [], 0
            for s, col in runs(n):
                t = ax.text(0, 0, s, ha='left', va='center', fontsize=10,
                            fontweight='bold', color=col, zorder=5)
                w, h = wh(t)
                rtexts.append(t); rws.append(w); lab_h = max(lab_h, h)
            lab_w = sum(rws)

            sub_t = None
            rl = support_label(n)
            if rl:
                sub_t = ax.text(0, 0, rl, ha='center', va='center', fontsize=8,
                                style='italic', color='#566573', zorder=5)
                sub_w, sub_h = wh(sub_t)

            if sub_t:
                gap = 0.07
                total = lab_h + gap + sub_h
                y_main, y_sub = py + total/2 - lab_h/2, py - total/2 + sub_h/2
                box_w, box_h = max(lab_w, sub_w), total
            else:
                y_main, box_w, box_h = py, lab_w, lab_h

            cx = px - lab_w / 2
            for t, w in zip(rtexts, rws):
                t.set_position((cx, y_main)); cx += w
            if sub_t:
                sub_t.set_position((px, y_sub))

            w, h = box_w + 0.46, box_h + 0.34
            boxes[n] = (w, h)

            lost = n in stemma.lost
            fc = '#e8edf0' if lost else '#5dade2'
            ec, lw, ls = '#1b2631', 1.4, 'solid'
            if lost:
                ec, ls = GREY, (0, (4, 3))
            if n in stemma.convergence:
                ec, lw, ls = FLAG_RED, 2.8, 'solid'
            if n in conflict_nodes:
                ec, lw, ls = WARN, 3.0, 'solid'
            ax.add_patch(FancyBboxPatch((px - w/2, py - h/2), w, h,
                         boxstyle='round,pad=0.02,rounding_size=0.10',
                         linewidth=lw, edgecolor=ec, facecolor=fc,
                         linestyle=ls, zorder=2))

            marks = []
            if n in conflict_nodes:
                marks.append((WARN, "\u26a0"))
            if lost:
                marks.append((GREY, "\u2736"))
            cxm = px + w/2 - 0.05
            for col, sym in marks:
                ax.text(cxm, py + h/2 + 0.04, sym, ha='right', va='bottom',
                        fontsize=12, color=col, zorder=6)
                cxm -= 0.42

        for p, c in edges:
            (x1, y1), (x2, y2) = pos[p], pos[c]
            ax.add_patch(FancyArrowPatch(
                (x1, y1 - boxes[p][1] / 2), (x2, y2 + boxes[c][1] / 2),
                arrowstyle='-|>', mutation_scale=15, color='#34495e',
                lw=1.5, shrinkA=2, shrinkB=2, zorder=1))

        def _loci(n):
            return f"{n} locus" if n == 1 else f"{n} loci"

        if stemma.contaminated:
            lines = ["\u26a0 Contamination \u2014 taux d'accords contre les autres:"]
            for m in sorted(stemma.contaminated):
                grps = sorted({g for a, b, sh in stemma.conflicts if m in sh
                               for g in (a, b)})
                counts = [stemma.support.get(g, 0) for g in grps]
                strong = all(c >= 2 for c in counts)
                verdict = "forte" if strong else "faible (coïncidence/polygénèse possible)"
                detail = "   vs   ".join(
                    f"{g} ({_loci(stemma.support.get(g, 0))})" for g in grps)
                lines.append(f"\u2022 {m}:  {detail}   \u2192 {verdict}")
            ax.text(0.01, 0.99, "\n".join(lines), transform=ax.transAxes,
                    va='top', ha='left', fontsize=9, color='#7b241c',
                    bbox=dict(boxstyle='round,pad=0.5', fc='#fdebd0',
                              ec=WARN, lw=1.4))

        ax.set_title("Stemma \u00b7 appui en loci par n\u0153ud \u00b7 rouge = t\u00e9moin "
                     "contaminant \u00b7 orange = groupes en conflit \u00b7 "
                     "gris pointill\u00e9 \u2736 = perdu / hypoth\u00e9tique", fontsize=9)

        if save_path:
            fig.savefig(save_path, dpi=150, bbox_inches='tight')
            plt.close(fig)
            return save_path
        plt.show()
        return fig, ax

    def reconstruct_archetype(self) -> None:
        """Reconstruct archetype?"""



if __name__ == "__main__":
    import unittest
    import os

    groups = StemmaPerGroup(ListApparatus.from_csv("test_alignment_2.csv"))
                
    print(ListApparatus.from_csv("test_alignment_2.csv"))
    print(groups.display_group_inclusion())

    # def test_build_stemma(self):
    #     """Tests that building the stemma behaves as expected."""
    print(groups.plot_stemma())