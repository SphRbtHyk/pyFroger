#!/usr/bin/env python3
"""
check_nexus.py — does the matrix match what `dimensions` promises?

PAUP reads a standard-datatype matrix as a fixed-length list: ntax rows of
nchar states. Any row that is short makes PAUP overrun into the next line and
swallow its taxon label as data, so the taxon count silently collapses and you
get, several rows later:

    Error(#608): A semicolon (;) was found where a taxon name was expected.

The offending row is NOT the one named in the error. This script finds it.

    python check_nexus.py heinrichi.nex
"""

import re
import sys
from collections import Counter
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 2:
        sys.exit("usage: check_nexus.py <file.nex>")
    p = Path(sys.argv[1])
    text = p.read_text(encoding="utf-8", errors="replace")

    # strip [comments]
    text = re.sub(r"\[[^\]]*\]", "", text)

    ntax = nchar = None
    if (m := re.search(r"\bntax\s*=\s*(\d+)", text, re.I)):
        ntax = int(m.group(1))
    if (m := re.search(r"\bnchar\s*=\s*(\d+)", text, re.I)):
        nchar = int(m.group(1))
    print(f"declared: ntax={ntax}  nchar={nchar}")

    m = re.search(r"\bmatrix\b(.*?)(?<![\w?])\;", text, re.I | re.S)
    if not m:
        sys.exit("no matrix ... ; found — is the terminator missing entirely?")
    body = m.group(1)

    rows = []
    for lineno, line in enumerate(body.splitlines(), 1):
        s = line.strip()
        if not s:
            continue
        parts = s.split()
        if len(parts) == 1:
            print(f"  !! line {lineno}: a lone token {parts[0]!r} — a label "
                  f"whose sequence went missing, or a wrapped line?")
            continue
        label, seq = parts[0], "".join(parts[1:])
        rows.append((lineno, label, seq))

    print(f"parsed  : {len(rows)} rows")
    if ntax is not None and len(rows) != ntax:
        print(f"  !! ROW COUNT MISMATCH: {len(rows)} rows vs ntax={ntax}")

    bad = [(ln, lab, len(s)) for ln, lab, s in rows if nchar and len(s) != nchar]
    if bad:
        print(f"\n  !! {len(bad)} ROW(S) DO NOT MATCH nchar={nchar}:")
        for ln, lab, n in bad:
            print(f"       line {ln:>4}  {lab:>4}  {n} chars  "
                  f"({'short' if n < nchar else 'long'} by {abs(nchar - n)})")
    else:
        print(f"  all rows are exactly {nchar} characters")

    dup = [lab for lab, c in Counter(lab for _, lab, _ in rows).items() if c > 1]
    if dup:
        print(f"\n  !! DUPLICATE LABELS: {dup}")

    # symbols declared vs symbols used
    syms = None
    if (m := re.search(r'\bsymbols\s*=\s*"([^"]*)"', text, re.I)):
        syms = set(m.group(1).replace(" ", ""))
    used = {c for _, _, s in rows for c in s}
    if syms:
        stray = sorted(used - syms - {"?", "-"})
        if stray:
            print(f"\n  !! SYMBOLS NOT DECLARED in symbols=: {stray}")
            for ch in stray:
                where = [(ln, lab) for ln, lab, s in rows if ch in s][:5]
                print(f"       {ch!r} appears in {where}")

    if not bad and not dup and (ntax is None or len(rows) == ntax):
        print("\nmatrix is consistent with its dimensions.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())