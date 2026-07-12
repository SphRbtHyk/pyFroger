#!/usr/bin/env python3
"""
prepare.py — derive the three tool-specific inputs from one NEXUS matrix.

The three programs cannot eat the same file:

  IQ-TREE   multistate morphology, up to 32 states  -> PHYLIP (states as-is)
  MrBayes   datatype=standard, AT MOST 10 STATES    -> recoded NEXUS + MB block
  FastME    cannot compute distances from multistate characters
                                                    -> precomputed PHYLIP matrix

MrBayes recoding
----------------
Characters with <= 10 observed states are renumbered 0..9 in order of first
appearance. Characters with more are handled per --many-states:

  binarise (default)  one 0/1 character per observed reading (presence of that
                      reading), which is the usual stemmatological move and
                      loses no grouping information, only the "these are
                      alternatives" constraint
  drop                the character is discarded

Distances for FastME
--------------------
Pairwise-complete Hamming: for each pair, count disagreements over the loci
where BOTH witnesses are extant, divide by that count. Missing data never
contributes a match or a mismatch. This is the whole ballgame for a tradition
as lacunose as this one.

Usage
-----
    python prepare.py heinrichi.nex --outdir work/
"""

from __future__ import annotations

import argparse
import itertools
import json
import re
import sys
from pathlib import Path


# ── NEXUS reading (standard datatype, non-interleaved or interleaved) ─────

def strip_comments(text: str) -> str:
    out, depth, quote = [], 0, False
    for ch in text:
        if quote:
            out.append(ch)
            if ch == "'":
                quote = False
        elif ch == "'" and depth == 0:
            quote = True
            out.append(ch)
        elif ch == "[":
            depth += 1
        elif ch == "]":
            depth = max(0, depth - 1)
        elif depth == 0:
            out.append(ch)
    return "".join(out)


TOKEN_RE = re.compile(r"'([^']*)'|(\S+)")


def read_nexus(path: Path) -> tuple[list[str], dict[str, list[str]], str, str]:
    text = strip_comments(path.read_text(encoding="utf-8", errors="replace"))

    m = re.search(r"begin\s+(?:data|characters)\s*;(.*?)\bend\s*;",
                  text, re.IGNORECASE | re.DOTALL)
    if not m:
        sys.exit("error: no data/characters block")
    block = m.group(1)

    missing, gap = "?", "-"
    fmt = re.search(r"\bformat\b(.*?);", block, re.IGNORECASE | re.DOTALL)
    if fmt:
        if (s := re.search(r"\bmissing\s*=\s*(\S)", fmt.group(1), re.I)):
            missing = s.group(1)
        if (s := re.search(r"\bgap\s*=\s*(\S)", fmt.group(1), re.I)):
            gap = s.group(1)

    mat = re.search(r"\bmatrix\b(.*?);", block, re.IGNORECASE | re.DOTALL)
    if not mat:
        sys.exit("error: no matrix command")

    data: dict[str, list[str]] = {}
    order: list[str] = []
    for line in mat.group(1).splitlines():
        line = line.strip()
        if not line:
            continue
        toks = [a or b for a, b in TOKEN_RE.findall(line)]
        if len(toks) < 2:
            continue
        label = toks[0].strip("'")
        seq = "".join(toks[1:])
        if label not in data:
            data[label] = []
            order.append(label)
        data[label].extend(list(seq))

    lens = {len(v) for v in data.values()}
    if len(lens) != 1:
        sys.exit(f"error: ragged matrix, row lengths {sorted(lens)}")

    return order, data, missing, gap


# ── Outputs ───────────────────────────────────────────────────────────────

def write_phylip_morph(taxa, data, missing, gap, out: Path) -> None:
    """Sequential PHYLIP with the original multistate symbols, for IQ-TREE."""
    n, m = len(taxa), len(data[taxa[0]])
    width = max(len(t) for t in taxa) + 2
    with out.open("w") as fh:
        fh.write(f"{n} {m}\n")
        for t in taxa:
            seq = "".join("?" if s in (missing, gap) else s for s in data[t])
            fh.write(f"{t.ljust(width)}{seq}\n")


def recode_for_mrbayes(taxa, data, missing, gap, many_states: str):
    """Return (new_data, nchar, report) with at most 10 states per character."""
    m = len(data[taxa[0]])
    new: dict[str, list[str]] = {t: [] for t in taxa}
    kept, binarised, dropped = 0, 0, 0

    for j in range(m):
        col = {t: data[t][j] for t in taxa}
        obs = [s for s in dict.fromkeys(col[t] for t in taxa)
               if s not in (missing, gap)]

        if len(obs) <= 10:
            mapping = {s: str(k) for k, s in enumerate(obs)}
            for t in taxa:
                s = col[t]
                new[t].append("?" if s in (missing, gap) else mapping[s])
            kept += 1
        elif many_states == "binarise":
            for reading in obs:
                for t in taxa:
                    s = col[t]
                    if s in (missing, gap):
                        new[t].append("?")
                    else:
                        new[t].append("1" if s == reading else "0")
            binarised += 1
        else:
            dropped += 1

    nchar = len(new[taxa[0]])
    report = {"characters_in": m, "characters_out": nchar,
              "renumbered": kept, "binarised": binarised, "dropped": dropped}
    return new, nchar, report


def write_mrbayes_nexus(taxa, new, nchar, out: Path, params: dict) -> None:
    width = max(len(t) for t in taxa) + 2
    with out.open("w") as fh:
        fh.write("#NEXUS\n\nbegin data;\n")
        fh.write(f"\tdimensions ntax={len(taxa)} nchar={nchar};\n")
        fh.write('\tformat datatype=standard symbols="0123456789" '
                 'missing=? gap=- interleave=no;\n')
        fh.write("\tmatrix\n")
        for t in taxa:
            fh.write(f"\t\t{t.ljust(width)}{''.join(new[t])}\n")
        fh.write("\t;\nend;\n\n")

        fh.write("begin mrbayes;\n")
        fh.write("\tset autoclose=yes nowarn=yes;\n")
        # coding=variable would be wrong: constant characters ARE in the
        # matrix, so no ascertainment-bias correction.
        fh.write("\tlset coding=all rates=gamma;\n")
        fh.write("\tprset ratepr=variable;\n")
        fh.write(f"\tmcmc ngen={params['ngen']} nruns={params['nruns']} "
                 f"nchains={params['nchains']} "
                 f"samplefreq={params['samplefreq']} "
                 f"printfreq={params['samplefreq']} "
                 f"diagnfreq={10 * params['samplefreq']} "
                 f"file=mrbayes;\n")
        fh.write(f"\tsumt burninfrac={params['burninfrac']} contype=allcompat;\n")
        fh.write(f"\tsump burninfrac={params['burninfrac']};\n")
        fh.write("end;\n")


def pairwise_complete_distances(taxa, data, missing, gap):
    """Hamming with pairwise-complete deletion. Returns (matrix, coverage)."""
    m = len(data[taxa[0]])
    lac = {missing, gap}
    n = len(taxa)
    D = [[0.0] * n for _ in range(n)]
    C = [[0] * n for _ in range(n)]
    warnings = []

    for a, b in itertools.combinations(range(n), 2):
        ta, tb = taxa[a], taxa[b]
        diff = shared = 0
        for j in range(m):
            sa, sb = data[ta][j], data[tb][j]
            if sa in lac or sb in lac:
                continue
            shared += 1
            if sa != sb:
                diff += 1
        if shared == 0:
            warnings.append((ta, tb))
            d = 1.0            # no overlap at all: maximally distant
        else:
            d = diff / shared
        D[a][b] = D[b][a] = d
        C[a][b] = C[b][a] = shared

    return D, C, warnings


def write_phylip_matrix(taxa, D, out: Path) -> None:
    width = max(len(t) for t in taxa) + 2
    with out.open("w") as fh:
        fh.write(f"{len(taxa)}\n")
        for i, t in enumerate(taxa):
            row = " ".join(f"{d:.6f}" for d in D[i])
            fh.write(f"{t.ljust(width)}{row}\n")


# ── main ──────────────────────────────────────────────────────────────────

def main(argv=None) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("nexus", type=Path)
    p.add_argument("--outdir", type=Path, default=Path("work"))
    p.add_argument("--many-states", choices=("binarise", "drop"),
                   default="binarise",
                   help="what to do with characters having >10 states "
                        "(MrBayes cannot represent them)")
    p.add_argument("--ngen", type=int, default=2_000_000)
    p.add_argument("--nruns", type=int, default=2)
    p.add_argument("--nchains", type=int, default=4)
    p.add_argument("--samplefreq", type=int, default=1000)
    p.add_argument("--burninfrac", type=float, default=0.25)
    args = p.parse_args(argv)

    args.outdir.mkdir(parents=True, exist_ok=True)

    taxa, data, missing, gap = read_nexus(args.nexus)
    nchar = len(data[taxa[0]])
    print(f"read {len(taxa)} taxa x {nchar} characters", file=sys.stderr)

    states = {s for row in data.values() for s in row} - {missing, gap}
    print(f"observed states: {len(states)} -> {''.join(sorted(states))}",
          file=sys.stderr)

    # IQ-TREE
    phy = args.outdir / "heinrichi.morph.phy"
    write_phylip_morph(taxa, data, missing, gap, phy)
    print(f"wrote {phy}", file=sys.stderr)

    # MrBayes
    new, mb_nchar, report = recode_for_mrbayes(
        taxa, data, missing, gap, args.many_states)
    mbnex = args.outdir / "heinrichi.mb.nex"
    write_mrbayes_nexus(taxa, new, mb_nchar, mbnex, {
        "ngen": args.ngen, "nruns": args.nruns, "nchains": args.nchains,
        "samplefreq": args.samplefreq, "burninfrac": args.burninfrac,
    })
    print(f"wrote {mbnex}  ({json.dumps(report)})", file=sys.stderr)

    # FastME
    D, C, warn = pairwise_complete_distances(taxa, data, missing, gap)
    dist = args.outdir / "heinrichi.dist"
    write_phylip_matrix(taxa, D, dist)
    print(f"wrote {dist}", file=sys.stderr)

    # Coverage report: which pairs rest on almost no shared evidence?
    thin = []
    for i, a in enumerate(taxa):
        for j in range(i + 1, len(taxa)):
            if C[i][j] < 50:
                thin.append((a, taxa[j], C[i][j]))
    if thin:
        print(f"warning: {len(thin)} pairs share <50 loci; "
              f"their distances are near-meaningless:", file=sys.stderr)
        for a, b, c in sorted(thin, key=lambda x: x[2])[:15]:
            print(f"    {a:>4} {b:>4}  {c} shared loci", file=sys.stderr)
    for a, b in warn:
        print(f"warning: {a} and {b} share NO extant locus", file=sys.stderr)

    (args.outdir / "prepare_report.json").write_text(json.dumps({
        "taxa": taxa, "nchar": nchar, "states": sorted(states),
        "mrbayes": report,
        "thin_pairs": [[a, b, c] for a, b, c in thin],
    }, indent=2))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
