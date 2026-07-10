"""Format input matrix to the different requirements.
"""
from __future__ import annotations
import argparse
import json
import sys
from collections import OrderedDict


MORPH_ALPHABET = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
DEFAULT_MISSING = set("?-N.Xx")


def read_fasta(path):
    seqs = OrderedDict()
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].strip(), []
            elif line:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
        raise ValueError(
            f"unequal sequence lengths in {path}: {sorted(lengths)}")
    return seqs


def recode(seqs, missing=DEFAULT_MISSING):
    """Compact per-column states to MORPH_ALPHABET; missing -> '?'.
    Returns (recoded OrderedDict, global_max_states, per_col_states)."""
    names = list(seqs)
    ncol = len(next(iter(seqs.values())))
    rows = {n: [] for n in names}
    per_col = []
    gmax = 0
    for c in range(ncol):
        mapping = {}
        for n in names:
            ch = seqs[n][c]
            if ch in missing:
                continue
            if ch not in mapping:
                if len(mapping) >= len(MORPH_ALPHABET):
                    raise ValueError(
                        f"column {c}: >{len(MORPH_ALPHABET)} states")
                mapping[ch] = MORPH_ALPHABET[len(mapping)]
        per_col.append(len(mapping))
        gmax = max(gmax, len(mapping))
        for n in names:
            ch = seqs[n][c]
            rows[n].append("?" if ch in missing else mapping[ch])
    recoded = OrderedDict((n, "".join(rows[n])) for n in names)
    return recoded, gmax, per_col


def name_map(sigla, width=4):
    fwd = OrderedDict((s, f"t{i:0{width}d}") for i, s in enumerate(sigla))
    rev = {v: k for k, v in fwd.items()}
    return fwd, rev


def distance_matrix(recoded):
    """Normalized Hamming, pairwise deletion of '?'.
    Returns (names, dict[(i,j)]->d)."""
    names = list(recoded)
    seq = {n: recoded[n] for n in names}
    D = {}
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = seq[names[i]], seq[names[j]]
            comp = mism = 0
            for x, y in zip(a, b):
                if x == "?" or y == "?":
                    continue
                comp += 1
                if x != y:
                    mism += 1
            if comp == 0:
                sys.stderr.write(
                    f"WARN: no comparable sites for {names[i]},{names[j]} -> d=1.0\n")
                d = 1.0
            else:
                d = mism / comp
            D[(i, j)] = D[(j, i)] = d
    return names, D


def mean_pairwise(names, D):
    n = len(names)
    if n < 2:
        return 0.0
    return sum(D[(i, j)] for i in range(n) for j in range(i + 1, n)) / (n * (n - 1) / 2)


def write_phylip_matrix(path, order_ids, D):
    n = len(order_ids)
    w = max(10, max(len(x) for x in order_ids) + 2)
    with open(path, "w") as fh:
        fh.write(f" {n}\n")
        for i, ni in enumerate(order_ids):
            cells = " ".join(f"{D.get((i, j), 0.0):.6f}" for j in range(n))
            fh.write(f"{ni.ljust(w)}{cells}\n")


def write_phylip_aln(path, ids_seqs):
    ids = list(ids_seqs)
    n = len(ids)
    m = len(ids_seqs[ids[0]])
    w = max(10, max(len(x) for x in ids) + 2)
    with open(path, "w") as fh:
        fh.write(f" {n} {m}\n")
        for i in ids:
            fh.write(f"{i.ljust(w)}{ids_seqs[i]}\n")


def write_nexus_standard(path, ids_seqs, gmax, outgroup_id, ngen, samplefreq,
                         stopval=0.01, diagnfreq=10000):
    ids = list(ids_seqs)
    n = len(ids)
    m = len(ids_seqs[ids[0]])
    symbols = MORPH_ALPHABET[:max(2, gmax)]
    stop = (f" mcmcdiagn=yes stoprule=yes stopval={stopval} "
            "relburnin=yes burninfrac=0.25") if stopval and stopval > 0 else ""
    lines = [
        "#NEXUS", "begin data;",
        f"  dimensions ntax={n} nchar={m};",
        f'  format datatype=standard missing=? gap=- symbols="{symbols}";',
        "  matrix",
    ]
    lines += [f"  {i}  {ids_seqs[i]}" for i in ids]
    lines += ["  ;", "end;", "begin mrbayes;",
              "  set autoclose=yes nowarn=yes;",
              f"  outgroup {outgroup_id};",
              "  lset coding=informative rates=gamma;",
              f"  mcmc ngen={ngen} samplefreq={samplefreq} nchains=4 nruns=2 "
              f"printfreq=10000 diagnfreq={diagnfreq}{stop};",
              "  sumt relburnin=yes burninfrac=0.25 conformat=simple;",
              "  sump relburnin=yes burninfrac=0.25;",
              "end;"]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def prepare(fasta, outgroup="w1", missing=DEFAULT_MISSING):
    """One-stop: returns everything the orchestrator needs for a single collation."""
    seqs = read_fasta(fasta)
    if outgroup not in seqs:
        raise ValueError(
            f"outgroup {outgroup!r} not among sigla: {list(seqs)[:5]}...")
    recoded, gmax, per_col = recode(seqs, missing)
    fwd, rev = name_map(list(seqs))
    ids_seqs = OrderedDict((fwd[s], recoded[s]) for s in seqs)
    # names == list(recoded) == original order
    names, D = distance_matrix(recoded)
    id_order = [fwd[s] for s in names]
    return {
        "ids_seqs": ids_seqs, "id_order": id_order, "D": D,
        "fwd": fwd, "rev": rev, "gmax": gmax,
        "outgroup_id": fwd[outgroup],
        "n_taxa": len(seqs), "n_chars": len(recoded[next(iter(recoded))]),
        "mean_hamming": mean_pairwise(names, D),
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Emit FastME matrix / IQ-TREE phylip / MrBayes nexus from a collation FASTA.")
    ap.add_argument("fasta")
    ap.add_argument("--outgroup", default="w1")
    ap.add_argument("--prefix", default="out")
    ap.add_argument("--ngen", type=int, default=1_000_000)
    ap.add_argument("--samplefreq", type=int, default=1000)
    a = ap.parse_args()
    p = prepare(a.fasta, a.outgroup)
    write_phylip_matrix(f"{a.prefix}.dist.phy", p["id_order"], p["D"])
    write_phylip_aln(f"{a.prefix}.aln.phy", p["ids_seqs"])
    write_nexus_standard(f"{a.prefix}.nex", p["ids_seqs"], p["gmax"],
                         p["outgroup_id"], a.ngen, a.samplefreq)
    json.dump({"fwd": p["fwd"], "rev": p["rev"], "outgroup_id": p["outgroup_id"]},
              open(f"{a.prefix}.namemap.json", "w"), indent=0)
    sys.stderr.write(
        f"{a.fasta}: {p['n_taxa']} taxa x {p['n_chars']} chars, "
        f"gmax_states={p['gmax']}, mean_hamming={p['mean_hamming']:.4f}\n")
