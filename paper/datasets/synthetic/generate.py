"""
stemma_benchmark.py
===================

Generator of *artificial textual traditions* with controllable parameters,
for benchmarking stemma-reconstruction methods (Froger / compatibility vs.
NJ / parsimony / Bayesian).

Outputs (mutually consistent):
  * CSV     -> wide collation table with a ``correct`` (reference) column,
              for pyFroger's errors-only apparatus logic.
  * FASTA   -> aligned character matrix (1 witness = 1 sequence,
              1 locus = 1 column), for phylogenetic tools.
  * Newick  -> the *rooted* true stemma, extant witnesses at their real
              position, lost ancestors as Greek letters; for RF, internal-
              ancestor identification, and Roos orientation scoring.

Transmission model
------------------
A rooted genealogy is grown; the root is the (always lost) archetype. Each
non-root node copies its parent and introduces errors locus by locus. The
per-locus alphabet has ``n_states`` states (0 = correct); a finite alphabet
makes two independent errors coincide with probability ~1/(n_states-1) =
*polygenesis*. Optional contamination, reversion (homoplasy) and lacunae.

Two extra controls:
  * ``n_generations`` : exact number of copying generations from the
    archetype (metric/transmission depth). A guaranteed spine fixes the
    depth; it is protected from leaf-pruning. If None, the tree is grown by
    uniform random attachment (depth emergent).
  * hot loci : a fraction ``hot_loci_fraction`` of sites are more prone to
    change (``hot_loci_factor`` x the error rate) and, optionally, draw from
    a smaller alphabet (``hot_loci_states``) so polygenesis concentrates on
    them. Which loci are hot is recorded as ground truth.

No third-party dependency (stdlib only).
"""

from __future__ import annotations

import csv
import json
import math
import string
from collections import defaultdict, deque
from dataclasses import dataclass, field
from random import Random
from typing import Optional


@dataclass
class _Node:
    id: int
    parent: Optional[int] = None
    children: list[int] = field(default_factory=list)
    extant: bool = True
    name: str = ""
    gen: int = 0
    text: list[int] = field(default_factory=list)
    errors_on_edge: int = 0


_GREEK = list("αβγδεζηθικλμνξοπρστυφχψω")


def _greek_label(i: int) -> str:
    if i < len(_GREEK):
        return _GREEK[i]
    return f"{_GREEK[i % len(_GREEK)]}{i // len(_GREEK)}"


@dataclass
class StemmaBenchmark:
    # --- size ---
    n_witnesses: int = 12
    n_loci: int = 200

    # --- corruption ---
    error_rate: float = 0.05
    n_states: int = 4
    allow_reversion: bool = False

    # --- topology / depth / loss ---
    loss_rate: float = 0.30
    max_children: Optional[int] = 3
    n_generations: Optional[int] = None   # exact copying-generation depth
    extant_archetype: bool = True         # archetype survives as a real witness (the root)

    # --- hot loci (site rate heterogeneity) ---
    hot_loci_fraction: float = 0.0        # share of sites that are "weak points"
    hot_loci_factor: float = 5.0          # error-rate multiplier on hot loci
    hot_loci_states: Optional[int] = None  # smaller alphabet on hot loci (more polygenesis)

    # --- horizontal transmission ---
    contamination_rate: float = 0.0
    contamination_strength: float = 0.30

    # --- missing data ---
    lacuna_rate: float = 0.0

    # --- reproducibility ---
    seed: Optional[int] = None

    def __post_init__(self) -> None:
        if self.n_states < 2:
            raise ValueError("n_states must be >= 2")
        if not 0.0 <= self.loss_rate < 1.0:
            raise ValueError("loss_rate must be in [0, 1)")
        if self.n_generations is not None and self.n_generations < 1:
            raise ValueError("n_generations must be >= 1")
        if (self.extant_archetype and self.n_generations is not None
                and self.n_witnesses < 2):
            raise ValueError("extant_archetype with n_generations needs "
                             "n_witnesses >= 2 (root + >=1 deep survivor)")
        if not 0.0 <= self.hot_loci_fraction <= 1.0:
            raise ValueError("hot_loci_fraction must be in [0, 1]")
        if self.hot_loci_states is not None and self.hot_loci_states < 2:
            raise ValueError("hot_loci_states must be >= 2")
        self._rng = Random(self.seed)
        self.nodes: dict[int, _Node] = {}
        self.root: int = 0
        self.witnesses: list[str] = []
        self.name2id: dict[str, int] = {}
        self.contaminations: list[dict] = []
        self.lacunae: dict[int, set[int]] = {}
        self.hot_loci: set[int] = set()
        self._locus_rate: list[float] = []
        self._locus_states: list[int] = []
        self._births: dict[tuple[int, int], list[int]] = defaultdict(list)
        self.generated = False

    # ===================================================================== #
    # generation
    # ===================================================================== #
    def generate(self) -> "StemmaBenchmark":
        self._build_tree()
        self._prune_lost_leaves()
        self._name_nodes()
        self._setup_loci()
        self._simulate_text()
        self._apply_contamination()
        self._apply_lacunae()
        self.generated = True
        return self

    # ---- tree -------------------------------------------------------------
    def _n_lost_target(self) -> int:
        base = (round(self.n_witnesses * self.loss_rate / (1 - self.loss_rate))
                if self.loss_rate > 0 else 0)
        floor = 0 if self.extant_archetype else 1  # lost archetype counts as >=1
        return max(floor, base)

    def _build_tree(self) -> None:
        if self.n_generations is None:
            self._build_random()
        else:
            self._build_generational()

    def _build_random(self) -> None:
        n_total = self.n_witnesses + self._n_lost_target()
        self.nodes = {0: _Node(id=0, parent=None, gen=0)}
        for nid in range(1, n_total):
            if self.max_children is None:
                pool = list(self.nodes)
            else:
                pool = [k for k, v in self.nodes.items()
                        if len(v.children) < self.max_children]
            parent = self._rng.choice(pool)
            self.nodes[parent].children.append(nid)
            self.nodes[nid] = _Node(id=nid, parent=parent,
                                    gen=self.nodes[parent].gen + 1)
        self._assign_loss_random(n_total)

    def _assign_loss_random(self, n_total: int) -> None:
        non_root = [k for k in self.nodes if k != 0]
        self._rng.shuffle(non_root)
        self.nodes[0].extant = self.extant_archetype
        n_keep = self.n_witnesses - 1 if self.extant_archetype else self.n_witnesses
        keep = set(non_root[:n_keep])
        for k in non_root:
            self.nodes[k].extant = k in keep

    def _build_generational(self) -> None:
        g = self.n_generations
        n_total = max(self.n_witnesses + self._n_lost_target(), g + 1)
        # feasibility w.r.t. branching
        if self.max_children is not None:
            if self.max_children == 1:
                cap = g + 1
            else:
                cap = (self.max_children ** (g + 1) - 1) // (self.max_children - 1)
            if n_total > cap:
                raise ValueError(
                    f"infeasible: need {n_total} nodes but a depth-{g} tree with "
                    f"max_children={self.max_children} holds at most {cap}. "
                    f"Increase max_children, increase loss_rate, or decrease n_generations."
                )
        n_extant_non_root = (self.n_witnesses - 1 if self.extant_archetype
                             else self.n_witnesses)
        if n_extant_non_root > n_total - 1:
            raise ValueError("n_witnesses too large for this configuration")

        # guaranteed spine -> fixes depth exactly at g
        self.nodes = {0: _Node(id=0, parent=None, gen=0)}
        prev, nid = 0, 1
        for depth in range(1, g + 1):
            self.nodes[nid] = _Node(id=nid, parent=prev, gen=depth)
            self.nodes[prev].children.append(nid)
            prev, nid = nid, nid + 1
        # fill remaining nodes (never beyond generation g)
        while len(self.nodes) < n_total:
            eligible = [k for k, v in self.nodes.items()
                        if v.gen < g and (self.max_children is None
                                          or len(v.children) < self.max_children)]
            if not eligible:
                raise ValueError("ran out of branching capacity; increase max_children")
            p = self._rng.choice(eligible)
            self.nodes[nid] = _Node(id=nid, parent=p, gen=self.nodes[p].gen + 1)
            self.nodes[p].children.append(nid)
            nid += 1
        self._assign_loss_generational(g)

    def _assign_loss_generational(self, g: int) -> None:
        self.nodes[0].extant = self.extant_archetype
        non_root = [k for k in self.nodes if k != 0]
        # force one extant witness at the deepest generation so the spine
        # survives pruning and the depth is preserved
        deepest = [k for k in non_root if self.nodes[k].gen == g]
        forced = self._rng.choice(deepest)
        rest = [k for k in non_root if k != forced]
        self._rng.shuffle(rest)
        n_more = (self.n_witnesses - 1 if self.extant_archetype
                  else self.n_witnesses) - 1
        keep = {forced} | set(rest[: max(0, n_more)])
        for k in non_root:
            self.nodes[k].extant = k in keep

    def _prune_lost_leaves(self) -> None:
        changed = True
        while changed:
            changed = False
            for nid, node in list(self.nodes.items()):
                if nid == self.root:
                    continue
                if not node.extant and not node.children:
                    self.nodes[node.parent].children.remove(nid)
                    del self.nodes[nid]
                    changed = True

    def _name_nodes(self) -> None:
        w_idx, g_idx = 1, 0
        for nid in self._bfs(self.root):
            node = self.nodes[nid]
            if node.extant:
                node.name = f"w{w_idx}"
                w_idx += 1
            else:
                node.name = _greek_label(g_idx)
                g_idx += 1
        self.witnesses = sorted((n.name for n in self.nodes.values() if n.extant),
                                key=lambda s: int(s[1:]))
        self.name2id = {n.name: n.id for n in self.nodes.values()}

    def _bfs(self, start: int) -> list[int]:
        order, q = [], deque([start])
        while q:
            nid = q.popleft()
            order.append(nid)
            q.extend(self.nodes[nid].children)
        return order

    # ---- loci (hot-loci rate heterogeneity) -------------------------------
    def _setup_loci(self) -> None:
        n_hot = round(self.hot_loci_fraction * self.n_loci)
        self.hot_loci = set(self._rng.sample(range(self.n_loci), n_hot)) if n_hot else set()
        self._locus_rate = [
            min(1.0, self.error_rate * self.hot_loci_factor) if l in self.hot_loci
            else self.error_rate
            for l in range(self.n_loci)
        ]
        self._locus_states = [
            self.hot_loci_states if (l in self.hot_loci and self.hot_loci_states)
            else self.n_states
            for l in range(self.n_loci)
        ]

    # ---- text simulation --------------------------------------------------
    def _simulate_text(self) -> None:
        self.nodes[self.root].text = [0] * self.n_loci
        for nid in self._bfs(self.root):
            node = self.nodes[nid]
            if node.parent is None:
                continue
            text = list(self.nodes[node.parent].text)
            errors = 0
            for locus in range(self.n_loci):
                if self._rng.random() < self._locus_rate[locus]:
                    new = self._mutate(text[locus], locus)
                    if new != text[locus]:
                        text[locus] = new
                        errors += 1
                        if new != 0:
                            self._births[(locus, new)].append(nid)
            node.text = text
            node.errors_on_edge = errors

    def _mutate(self, current: int, locus: int) -> int:
        states = self._locus_states[locus]
        lo = 0 if self.allow_reversion else 1
        choices = [s for s in range(lo, states) if s != current]
        return self._rng.choice(choices) if choices else current

    # ---- contamination ----------------------------------------------------
    def _contam_spans(self) -> list[tuple[int, int]]:
        # contamination_strength now = expected fraction copied as ONE contiguous
        # passage (an exemplar switch), not a per-locus probability.
        length = max(1, min(self.n_loci, round(self.contamination_strength * self.n_loci)))
        start = self._rng.randint(0, self.n_loci - length)
        return [(start, start + length)]      # return >1 span for multi-switch scribes

    def _apply_contamination(self) -> None:
        if self.contamination_rate <= 0:
            return
        # snapshot: every event reads the *genealogical* (pre-contamination) text of
        # its source, so a source used here can't be retro-altered by a later event.
        snapshot = {nid: list(n.text) for nid, n in self.nodes.items()}
        for nid, node in self.nodes.items():
            if not node.extant or nid == self.root:
                continue
            if self._rng.random() >= self.contamination_rate:
                continue
            forbidden = self._ancestors(nid) | self._descendants(nid) | {nid}
            pool = [k for k in self.nodes if k not in forbidden]
            if not pool:
                continue
            source = self._rng.choice(pool)
            src = snapshot[source]
            spans, n_changed = [], 0
            for start, end in self._contam_spans():
                for l in range(start, end):
                    if node.text[l] != src[l]:
                        n_changed += 1
                    node.text[l] = src[l]
                spans.append([start, end])
            self.contaminations.append(
                {"witness": node.name, "source": self.nodes[source].name,
                "spans": spans, "span_len": sum(e - s for s, e in spans),
                "n_changed": n_changed})

    def _ancestors(self, nid: int) -> set[int]:
        out, cur = set(), self.nodes[nid].parent
        while cur is not None:
            out.add(cur)
            cur = self.nodes[cur].parent
        return out

    def _descendants(self, nid: int) -> set[int]:
        out, q = set(), deque(self.nodes[nid].children)
        while q:
            c = q.popleft()
            out.add(c)
            q.extend(self.nodes[c].children)
        return out

    def _is_ancestor(self, a: int, b: int) -> bool:
        return a in self._ancestors(b)

    # ---- lacunae ----------------------------------------------------------
    def _apply_lacunae(self) -> None:
        if self.lacuna_rate <= 0:
            return
        for nid, node in self.nodes.items():
            if not node.extant or nid == self.root:  # reference stays complete
                continue
            miss = {l for l in range(self.n_loci) if self._rng.random() < self.lacuna_rate}
            if miss:
                self.lacunae[nid] = miss

    # ===================================================================== #
    # output
    # ===================================================================== #
    @staticmethod
    def _state_token(state: int) -> str:
        return string.ascii_lowercase[state] if state < 26 else f"s{state}"

    def _cell(self, nid: int, locus: int, missing: str = "?") -> str:
        if locus in self.lacunae.get(nid, ()):
            return missing
        return self._state_token(self.nodes[nid].text[locus])

    def to_csv(self, path: str) -> str:
        self._require()
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["locus"] + self.witnesses)
            for locus in range(self.n_loci):
                row = [locus + 1]
                row += [self._cell(self.name2id[n], locus) for n in self.witnesses]
                w.writerow(row)
        return path

    def to_fasta(self, path: str, alphabet: Optional[str] = None, gap: str = "-") -> str:
        self._require()
        alpha = alphabet or (string.digits + string.ascii_uppercase + string.ascii_lowercase)
        max_states = max(self.n_states, self.hot_loci_states or 0)
        if max_states > len(alpha):
            raise ValueError("state count exceeds FASTA alphabet size; pass alphabet=")
        with open(path, "w", encoding="utf-8") as fh:
            for name in self.witnesses:
                nid = self.name2id[name]
                seq = (gap if locus in self.lacunae.get(nid, ())
                       else alpha[self.nodes[nid].text[locus]]
                       for locus in range(self.n_loci))
                fh.write(f">{name}\n{''.join(seq)}\n")
        return path

    def to_newick(self, path: Optional[str] = None, *,
                  branch_lengths: bool = False, grafted: bool = False) -> str:
        self._require()
        nwk = self._emit(self.root, branch_lengths, grafted) + ";"
        if path:
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(nwk + "\n")
        return nwk

    def _kept(self, nid: int) -> bool:
        node = self.nodes[nid]
        return node.extant or nid == self.root or len(node.children) >= 2

    def _kept_children(self, nid: int):
        result, stack = [], [(c, self.nodes[c].errors_on_edge)
                             for c in self.nodes[nid].children]
        while stack:
            cid, acc = stack.pop()
            if self._kept(cid):
                result.append((cid, acc))
            else:
                for gc in self.nodes[cid].children:
                    stack.append((gc, acc + self.nodes[gc].errors_on_edge))
        return result

    def _emit(self, nid: int, bl: bool, grafted: bool) -> str:
        node = self.nodes[nid]
        kids = self._kept_children(nid)
        if not kids:
            return node.name
        parts = [self._emit(c, bl, grafted) + (f":{ln}" if bl else "") for c, ln in kids]
        if grafted and node.extant:
            parts.append(node.name + (":0" if bl else ""))
            return "(" + ",".join(parts) + ")"
        return "(" + ",".join(parts) + ")" + node.name

    # ===================================================================== #
    # ground-truth metadata
    # ===================================================================== #
    def _full_depth_map(self) -> dict[int, int]:
        d = {self.root: 0}
        for nid in self._bfs(self.root):
            for c in self.nodes[nid].children:
                d[c] = d[nid] + 1
        return d

    def _topo_height(self) -> int:
        best = [0]
        def rec(nid: int, depth: int) -> None:
            kids = self._kept_children(nid)
            if not kids:
                best[0] = max(best[0], depth)
                return
            for c, _ in kids:
                rec(c, depth + 1)
        rec(self.root, 0)
        return best[0]

    def summary(self) -> dict:
        self._require()
        d = self._full_depth_map()
        lost = [n for n in self.nodes.values() if not n.extant]
        extant_internal = [n.name for n in self.nodes.values() if n.extant and n.children]
        root_witness = self.nodes[self.root].name if self.nodes[self.root].extant else None

        poly_loci, poly_states = set(), 0
        for (locus, state), origins in self._births.items():
            if len(origins) < 2:
                continue
            if any(not (self._is_ancestor(a, b) or self._is_ancestor(b, a))
                   for i, a in enumerate(origins) for b in origins[i + 1:]):
                poly_states += 1
                poly_loci.add(locus)
        poly_hot = len(poly_loci & self.hot_loci)
        poly_cold = len(poly_loci - self.hot_loci)
        n_cold = self.n_loci - len(self.hot_loci)

        return {
            "n_witnesses": len(self.witnesses),
            "n_loci": self.n_loci,
            "n_nodes_total": len(self.nodes),
            "n_lost_ancestors": len(lost),
            "loss_rate_realized": round(len(lost) / len(self.nodes), 3),
            "metric_height": max(d[self.name2id[w]] for w in self.witnesses),
            "metric_height_mean": round(
                sum(d[self.name2id[w]] for w in self.witnesses) / len(self.witnesses), 2),
            "topo_height": self._topo_height(),
            "n_generations_set": self.n_generations,
            "extant_archetype": self.extant_archetype,
            "root_witness": root_witness,  # pass as outgroup to phylogenetic tools
            "extant_internal_ancestors": extant_internal,
            "n_extant_internal_ancestors": len(extant_internal),
            "n_contaminated": len(self.contaminations),
            "contaminations": self.contaminations,
            "n_hot_loci": len(self.hot_loci),
            "hot_loci": sorted(self.hot_loci),
            "polygenetic_loci": len(poly_loci),
            "polygenetic_states": poly_states,
            "polygenesis_in_hot_loci": poly_hot,
            "polygenesis_in_cold_loci": poly_cold,
            "polygenesis_rate_hot": round(poly_hot / len(self.hot_loci), 3) if self.hot_loci else None,
            "polygenesis_rate_cold": round(poly_cold / n_cold, 3) if n_cold else None,
            "missing_cells": sum(len(s) for s in self.lacunae.values()),
            "params": {
                "error_rate": self.error_rate, "n_states": self.n_states,
                "allow_reversion": self.allow_reversion, "loss_rate": self.loss_rate,
                "max_children": self.max_children, "n_generations": self.n_generations,
                "extant_archetype": self.extant_archetype,
                "hot_loci_fraction": self.hot_loci_fraction,
                "hot_loci_factor": self.hot_loci_factor,
                "hot_loci_states": self.hot_loci_states,
                "contamination_rate": self.contamination_rate,
                "contamination_strength": self.contamination_strength,
                "lacuna_rate": self.lacuna_rate, "seed": self.seed,
            },
        }

    def write_all(self, prefix: str, *, grafted_newick: bool = False) -> dict:
        paths = {
            "csv": self.to_csv(prefix + ".csv"),
            "fasta": self.to_fasta(prefix + ".fasta"),
            "newick": self.to_newick(prefix + ".nwk", grafted=grafted_newick),
        }
        with open(prefix + ".meta.json", "w", encoding="utf-8") as fh:
            json.dump(self.summary(), fh, ensure_ascii=False, indent=2)
        paths["meta"] = prefix + ".meta.json"
        return paths

    def _require(self) -> None:
        if not self.generated:
            raise RuntimeError("call .generate() first")


if __name__ == "__main__":
    # Generate 10 traditions
    # No contamination
    for ix in range(10):
        bench = StemmaBenchmark(
            n_witnesses=15, n_loci=400, error_rate=0.1, n_states=4,
            n_generations=6, loss_rate=0.0, max_children=4,
            hot_loci_fraction=0.0, hot_loci_factor=6, hot_loci_states=2,
            contamination_rate=0.00, lacuna_rate=0.00
        ).generate()
        print(bench.to_newick(branch_lengths=False))
        print(json.dumps(bench.summary(), ensure_ascii=False, indent=2))
        bench.write_all(f"tmp/demo_tradition_type_1_{ix}")

    # Loss
    for ix in range(10):
        bench = StemmaBenchmark(
            n_witnesses=15, n_loci=400, error_rate=0.1, n_states=4,
            n_generations=6, loss_rate=0.3, max_children=4,
            hot_loci_fraction=0.0, hot_loci_factor=6, hot_loci_states=2,
            contamination_rate=0.00, lacuna_rate=0.00
        ).generate()
        print(bench.to_newick(branch_lengths=False))
        print(json.dumps(bench.summary(), ensure_ascii=False, indent=2))
        bench.write_all(f"tmp/demo_tradition_type_2_{ix}")

    # Contamination
    for ix in range(10):
        bench = StemmaBenchmark(
            n_witnesses=15, n_loci=400, error_rate=0.1, n_states=4,
            n_generations=6, loss_rate=0.0, max_children=4,
            hot_loci_fraction=0.0, hot_loci_factor=6, hot_loci_states=2,
            contamination_rate=0.05, lacuna_rate=0.00
        ).generate()
        print(bench.to_newick(branch_lengths=False))
        print(json.dumps(bench.summary(), ensure_ascii=False, indent=2))
        bench.write_all(f"tmp/demo_tradition_type_3_{ix}")