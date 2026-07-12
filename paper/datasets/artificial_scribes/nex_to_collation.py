#!/usr/bin/env python3
"""
nexus_to_collation.py — turn a NEXUS `standard` datatype matrix into a
collation CSV suitable for pyFroger.

Output shape (matches the Roos/Heikkilä-style collation):

    locus,A,Ab,Ac,...
    locus_1,A,C,C,...
    locus_2,B,B,B,...

  * one row per character, named `locus_<i>` (1-based)
  * one column per witness, named with the taxon label from the NEXUS
  * missing data (`?`) and gaps (`-`) become empty cells

Handles: interleaved or non-interleaved matrices, MESQUITE-style comments,
quoted/underscored taxon labels, `format` options in any order, and the
NEXUS defaults (missing=?, gap=-, datatype=standard).

Usage
-----
    python nexus_to_collation.py henry.nex -o collation.csv
    python nexus_to_collation.py henry.nex --drop Da I J
    python nexus_to_collation.py henry.nex --keep-missing '?'
    python nexus_to_collation.py henry.nex --constant-chars drop
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path


# --------------------------------------------------------------------------
# NEXUS parsing
# --------------------------------------------------------------------------

def strip_comments(text: str) -> str:
    """Remove NEXUS [square-bracket] comments, respecting 'quoted' strings."""
    out, depth, in_quote = [], 0, False
    for ch in text:
        if in_quote:
            out.append(ch)
            if ch == "'":
                in_quote = False
        elif ch == "'" and depth == 0:
            in_quote = True
            out.append(ch)
        elif ch == "[":
            depth += 1
        elif ch == "]":
            depth = max(0, depth - 1)
        elif depth == 0:
            out.append(ch)
    return "".join(out)


def find_block(text: str, name: str) -> str | None:
    """Return the body of `begin <name>; ... end;` (case-insensitive)."""
    pat = re.compile(
        rf"begin\s+{name}\s*;(.*?)\bend(?:block)?\s*;",
        re.IGNORECASE | re.DOTALL,
    )
    m = pat.search(text)
    return m.group(1) if m else None


def parse_dimensions(block: str) -> tuple[int | None, int | None]:
    m = re.search(r"\bdimensions\b(.*?);", block, re.IGNORECASE | re.DOTALL)
    if not m:
        return None, None
    body = m.group(1)
    ntax = re.search(r"\bntax\s*=\s*(\d+)", body, re.IGNORECASE)
    nchar = re.search(r"\bnchar\s*=\s*(\d+)", body, re.IGNORECASE)
    return (int(ntax.group(1)) if ntax else None,
            int(nchar.group(1)) if nchar else None)


def parse_format(block: str) -> dict:
    """NEXUS defaults per the standard, overridden by the format command."""
    fmt = {
        "datatype": "standard",
        "missing": "?",
        "gap": "-",
        "symbols": None,
        "interleave": False,
        "matchchar": None,
    }
    m = re.search(r"\bformat\b(.*?);", block, re.IGNORECASE | re.DOTALL)
    if not m:
        return fmt
    body = m.group(1)

    if (s := re.search(r"\bdatatype\s*=\s*(\w+)", body, re.IGNORECASE)):
        fmt["datatype"] = s.group(1).lower()
    if (s := re.search(r"\bmissing\s*=\s*(\S)", body, re.IGNORECASE)):
        fmt["missing"] = s.group(1)
    if (s := re.search(r"\bgap\s*=\s*(\S)", body, re.IGNORECASE)):
        fmt["gap"] = s.group(1)
    if (s := re.search(r"\bmatchchar\s*=\s*(\S)", body, re.IGNORECASE)):
        fmt["matchchar"] = s.group(1)
    if (s := re.search(r'\bsymbols\s*=\s*"([^"]*)"', body, re.IGNORECASE)):
        fmt["symbols"] = "".join(s.group(1).split())
    if re.search(r"\binterleave\b(?!\s*=\s*no)", body, re.IGNORECASE):
        fmt["interleave"] = True
    return fmt


TOKEN_RE = re.compile(r"'([^']*)'|(\S+)")


def parse_matrix(block: str, fmt: dict) -> dict[str, list[str]]:
    """
    Parse the `matrix` command into {taxon: [state, state, ...]}.

    Works for interleaved and non-interleaved files alike: we simply append
    each new chunk to whatever we already have for that taxon.
    """
    m = re.search(r"\bmatrix\b(.*?);", block, re.IGNORECASE | re.DOTALL)
    if not m:
        raise ValueError("no `matrix` command found in the data block")
    body = m.group(1)

    data: dict[str, list[str]] = {}
    order: list[str] = []

    for line in body.splitlines():
        line = line.strip()
        if not line:
            continue
        toks = [a or b for a, b in TOKEN_RE.findall(line)]
        if len(toks) < 2:
            # A label with its sequence on the following line is not legal
            # NEXUS; if we hit that, the file needs manual repair.
            continue
        label = toks[0].strip("'").replace("_", " ").strip()
        seq = "".join(toks[1:])

        states = expand_states(seq, fmt)

        if label not in data:
            data[label] = []
            order.append(label)
        data[label].extend(states)

    # Preserve first-seen taxon order.
    return {k: data[k] for k in order}


def expand_states(seq: str, fmt: dict) -> list[str]:
    """
    Split a raw sequence string into one entry per character.

    Polymorphisms/uncertainties written as (AB) or {AB} are kept as a single
    entry with the braces stripped, e.g. "AB" — pyFroger will see an unusual
    reading rather than silently losing the site. Everything else is one
    character per site.
    """
    states: list[str] = []
    i = 0
    while i < len(seq):
        ch = seq[i]
        if ch in "({":
            close = ")" if ch == "(" else "}"
            j = seq.find(close, i)
            if j == -1:
                raise ValueError(f"unclosed {ch!r} in matrix row")
            states.append(seq[i + 1:j])
            i = j + 1
        else:
            states.append(ch)
            i += 1
    return states


# --------------------------------------------------------------------------
# Conversion
# --------------------------------------------------------------------------

def to_collation(
    data: dict[str, list[str]],
    fmt: dict,
    missing_out: str = "",
    drop: set[str] | None = None,
    constant_chars: str = "keep",
) -> tuple[list[str], list[list[str]]]:
    """Return (header, rows) for the collation CSV."""
    drop = drop or set()
    witnesses = [w for w in data if w not in drop]
    if not witnesses:
        raise ValueError("every witness was dropped")

    lengths = {w: len(data[w]) for w in witnesses}
    if len(set(lengths.values())) != 1:
        raise ValueError(
            "ragged matrix — rows have different lengths:\n  "
            + "\n  ".join(f"{w}: {n}" for w, n in sorted(lengths.items()))
        )
    nchar = next(iter(lengths.values()))

    lacunae = {fmt["missing"], fmt["gap"]}

    rows: list[list[str]] = []
    for i in range(nchar):
        cells = []
        for w in witnesses:
            s = data[w][i]
            cells.append(missing_out if s in lacunae else s)

        if constant_chars == "drop":
            observed = {c for c in cells if c != missing_out}
            if len(observed) <= 1:
                continue

        rows.append([f"locus_{len(rows) + 1}" if constant_chars == "drop"
                     else f"locus_{i + 1}"] + cells)

    return ["locus"] + witnesses, rows


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="Convert a NEXUS standard-datatype matrix into a "
                    "pyFroger collation CSV."
    )
    p.add_argument("nexus", type=Path, help="input .nex file")
    p.add_argument("-o", "--output", type=Path, default=None,
                   help="output CSV (default: <input>.csv)")
    p.add_argument("--drop", nargs="*", default=[], metavar="TAXON",
                   help="witnesses to exclude (e.g. --drop Da I J)")
    p.add_argument("--keep-missing", default="", metavar="STR",
                   help="what to write for missing data "
                        "(default: empty cell; use '?' to keep the symbol)")
    p.add_argument("--constant-chars", choices=("keep", "drop"),
                   default="keep",
                   help="keep or drop characters with <=1 observed state "
                        "(default: keep; loci are renumbered if dropped)")
    args = p.parse_args(argv)

    text = strip_comments(args.nexus.read_text(encoding="utf-8", errors="replace"))

    block = find_block(text, "data")
    if block is None:
        block = find_block(text, "characters")
    if block is None:
        print("error: no `data` or `characters` block found", file=sys.stderr)
        return 1

    ntax, nchar = parse_dimensions(block)
    fmt = parse_format(block)

    if fmt["datatype"] not in ("standard", "protein", "dna", "rna", "nucleotide"):
        print(f"warning: unexpected datatype={fmt['datatype']!r}", file=sys.stderr)
    if fmt["matchchar"]:
        print("error: matchchar= is not supported; expand the matrix first",
              file=sys.stderr)
        return 1

    data = parse_matrix(block, fmt)

    # --- sanity checks against the declared dimensions -------------------
    if ntax is not None and len(data) != ntax:
        print(f"warning: dimensions say ntax={ntax} but {len(data)} taxa "
              f"were parsed", file=sys.stderr)
    if nchar is not None:
        bad = {w: len(s) for w, s in data.items() if len(s) != nchar}
        if bad:
            print(f"warning: dimensions say nchar={nchar}, but these rows "
                  f"differ:", file=sys.stderr)
            for w, n in sorted(bad.items()):
                print(f"    {w}: {n}", file=sys.stderr)

    if fmt["symbols"]:
        allowed = set(fmt["symbols"]) | {fmt["missing"], fmt["gap"]}
        stray = {s for row in data.values() for s in row
                 if len(s) == 1 and s not in allowed}
        if stray:
            print(f"warning: states outside symbols= : {sorted(stray)}",
                  file=sys.stderr)

    unknown_drop = set(args.drop) - set(data)
    if unknown_drop:
        print(f"warning: --drop names not in the matrix: {sorted(unknown_drop)}",
              file=sys.stderr)

    header, rows = to_collation(
        data, fmt,
        missing_out=args.keep_missing,
        drop=set(args.drop),
        constant_chars=args.constant_chars,
    )

    out = args.output or args.nexus.with_suffix(".csv")
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)

    print(f"wrote {out}: {len(rows)} loci x {len(header) - 1} witnesses")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())