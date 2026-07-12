import os
import argparse
import warnings
import multiprocessing as mp
from collections import defaultdict
import numpy as np

import dendropy
import pandas as pd

from pyfroger import evaluate_reconstructed_trees


# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
TRUTH_ROOT = "../../datasets/synthetic/generated_trees"
RECON_ROOT = "../experiment_outputs"
OUT_DIR = "./results"
N_REPLICATES = 100
METHODS = ["iqtree", "fastme", "froger","mrbayes"]

# rooted (clade) RF + Roos orientation only make sense for oriented recons.
# Froger's UnorientedStemma output is scored unrooted / orientation-free.
METHOD_CONFIG = {
    "iqtree":  {"rooted": True,  "orientation": True},
    "fastme":  {"rooted": True,  "orientation": True},
    "mrbayes": {"rooted": True,  "orientation": True},
    "froger":  {"rooted": True, "orientation": True},
}

# Define which metrics are better when higher vs lower
BETTER_HIGHER = {
    "clade_P": True, "clade_R": True, "clade_F1": True,
    "anc_P": True, "anc_R": True, "anc_F1": True
}
BETTER_LOWER = {
    "rf": True, "nRF": True,
    "skel%": True
}


def truth_path(scenario, basename):
    return os.path.join(TRUTH_ROOT, scenario, basename + ".nwk")


def recon_path(method, scenario, basename):
    # ADJUST HERE if run_bench.py uses a different layout. This mirrors the
    # {scenario}_{basename} convention implied by the draft.
    return os.path.join(RECON_ROOT, method, f"{scenario}_{basename}.nwk")


def sanitize(s):
    return s.replace(".", "_")


# --------------------------------------------------------------------------
# Experimental design grid  (mirrors generate_trees.py exactly)
# Each entry: parameter names, the value grid, a filename template, and a
# short display-label builder.
# --------------------------------------------------------------------------
DESIGN = {
    "scenario_1": {
        "grid": [{"error_rate": e} for e in (0.05, 0.1, 0.2, 0.3)],
        "fname": lambda p, ix: f"tree_{p['error_rate']}_{ix}",
        "label": lambda p: f"err={p['error_rate']}",
    },
    "scenario_2": {
        "grid": [{"loss_rate": l} for l in (0.1, 0.2, 0.5)],
        "fname": lambda p, ix: f"tree_{p['loss_rate']}_{ix}",
        "label": lambda p: f"loss={p['loss_rate']}",
    },
    "scenario_3": {
        "grid": [{"hot_loci_fraction": h} for h in (0.1, 0.2, 0.5)],
        "fname": lambda p, ix: f"tree_{p['hot_loci_fraction']}_{ix}",
        "label": lambda p: f"hot={p['hot_loci_fraction']}",
    },
    # "scenario_4": {
    #     "grid": [{"contamination_rate": c} for c in (0.1, 0.2, 0.3, 0.7)],
    #     "fname": lambda p, ix: f"tree_{p['contamination_rate']}_{ix}",
    #     "label": lambda p: f"cont={p['contamination_rate']}",
    # },
    # "scenario_5": {
    #     "grid": [{"contamination_rate": c, "hot_loci_fraction": h}
    #              for c in (0.1, 0.2, 0.3, 0.7) for h in (0.1, 0.2, 0.5)],
    #     "fname": lambda p, ix: (f"tree_{p['contamination_rate']}_"
    #                             f"{p['hot_loci_fraction']}_{ix}"),
    #     "label": lambda p: f"cont={p['contamination_rate']},hot={p['hot_loci_fraction']}",
    # },
    # "scenario_6": {
    #     "grid": [{"n_loci": n} for n in (50, 100, 500)],
    #     "fname": lambda p, ix: f"tree_{p['n_loci']}_{ix}",
    #     "label": lambda p: f"loci={p['n_loci']}",
    # },
}


METRIC_SPECS = [
    ("rf",        "robinson_foulds",    "rf_distance",                    "count"),
    ("nRF",       "robinson_foulds",    "normalised_rf",                  "rate"),
    ("clade_P",   "robinson_foulds",    "clade_precision",                "rate"),
    ("clade_R",   "robinson_foulds",    "clade_recall",                   "rate"),
    ("clade_F1",  "robinson_foulds",    "clade_f1",                       "rate"),
    ("skel%",     "roos_sign_distance", "skeleton_sign_distance",         "pct"),
    ("anc_P",     "internal_ancestors", "internal_ancestor_precision",    "rate"),
    ("anc_R",     "internal_ancestors", "internal_ancestor_recall",       "rate"),
    ("anc_F1",    "internal_ancestors", "internal_ancestor_f1",           "rate"),
]
METRIC_COLS = [m[0] for m in METRIC_SPECS]
KIND = {m[0]: m[3] for m in METRIC_SPECS}

DISPLAY_COLS = ["nRF", "clade_P", "clade_R", "clade_F1", "skel%"]

# Emoji indicators for pretty performance


def get_performance_indicator(value, metric):
    """Return emoji and color based on performance."""
    if pd.isna(value):
        return "—"

    if metric in BETTER_LOWER:
        # Lower is better
        if value < 0.2:
            return "🟢"
        elif value < 0.4:
            return "🟡"
        else:
            return "🔴"
    else:
        # Higher is better
        if value >= 0.8:
            return "🟢"
        elif value >= 0.6:
            return "🟡"
        else:
            return "🔴"


def print_pretty_header(text, width=80):
    """Print a formatted header with box drawing."""
    print("\n" + "=" * width)
    print(f"  {text}")
    print("=" * width)


def print_scenario_summary(scenario, method_rankings):
    """Print a summary of method rankings for a scenario."""
    print(f"\n📊 {scenario.upper()} - Method Performance Summary")
    print("-" * 60)

    # Show best method for each cell
    for cell, rankings in method_rankings.items():
        if rankings:
            best_method = rankings[0][0]
            best_score = rankings[0][1]
            print(f"  {cell:20} → Best: {best_method:10} (F1={best_score:.3f})")


def print_ranking_table(disp_df, scenario):
    """Print a ranked table with performance indicators."""
    print(f"\n📈 {scenario.upper()} - Detailed Results (Ranked by clade_F1)")
    print("-" * 80)

    # Create ranking for each cell
    for cell in disp_df.index.get_level_values('cell').unique():
        cell_data = disp_df.xs(cell, level='cell')

        # Sort by clade_F1 (descending)
        if 'clade_F1' in cell_data.columns:
            sorted_data = cell_data.sort_values('clade_F1', ascending=False)

            print(f"\n  🎯 {cell}")
            print(
                f"  {'Rank':<6} {'Method':<12} {'nRF':<12} {'Precision':<12} {'Recall':<12} {'F1 Score':<12} {'Performance'}")
            print(
                f"  {'-'*6} {'-'*12} {'-'*12} {'-'*12} {'-'*12} {'-'*12} {'-'*12}")

            for rank, (idx, row) in enumerate(sorted_data.iterrows(), 1):
                method = idx[0]  # method is first level after cell

                # Parse values (they are strings like "0.123±0.045")
                nrf_val = float(row['nRF'].split(
                    '±')[0]) if '±' in row['nRF'] else float(row['nRF'])
                prec_val = float(row['clade_P'].split(
                    '±')[0]) if '±' in row['clade_P'] else float(row['clade_P'])
                rec_val = float(row['clade_R'].split(
                    '±')[0]) if '±' in row['clade_R'] else float(row['clade_R'])
                f1_val = float(row['clade_F1'].split(
                    '±')[0]) if '±' in row['clade_F1'] else float(row['clade_F1'])
                skel_val = float(row['skel%'].split(
                    '±')[0]) if '±' in row['skel%'] else float(row['skel%'])
                # Determine ranking indicator
                if rank == 1:
                    indicator = "🥇"
                elif rank == 2:
                    indicator = "🥈"
                elif rank == 3:
                    indicator = "🥉"
                else:
                    indicator = f"  {rank}"

                # Get performance emoji
                perf_emoji = get_performance_indicator(f1_val, 'clade_F1')

                print(
                    f"  {indicator:<6} {method:<12} {nrf_val:<12.3f} {prec_val:<12.3f} {rec_val:<12.3f} {f1_val:<12.3f} {skel_val:<12.3f} {perf_emoji}")

            # Show gap to best
            if len(sorted_data) > 1:
                best_f1 = float(sorted_data.iloc[0]['clade_F1'].split('±')[0])
                worst_f1 = float(
                    sorted_data.iloc[-1]['clade_F1'].split('±')[0])
                gap = best_f1 - worst_f1
                print(
                    f"  {'':<6} {'':<12} {'':<12} {'':<12} {'':<12} Gap: {gap:.3f}")


def print_overall_summary(mean_df):
    """Print overall summary across all scenarios and metrics."""
    print("\n" + "=" * 80)
    print("🏆 OVERALL METHOD RANKING (Average across all scenarios, cells & metrics)")
    print("=" * 80)

    # Calculate average score for each method across all metrics
    method_scores = {}

    for method in METHODS:
        scores_by_metric = {}
        total_scores = []

        for metric in METRIC_COLS:
            metric_values = []
            for scenario in DESIGN:
                for cell in cell_order(scenario):
                    try:
                        key = (scenario, cell, method)
                        if key in mean_df.index:
                            val = mean_df.loc[key, metric]
                            if not pd.isna(val):
                                metric_values.append(val)
                    except:
                        pass

            if metric_values:
                # For metrics where lower is better, invert the score for ranking
                if metric in BETTER_LOWER:
                    # Normalize and invert: 1 - normalized_value
                    max_val = max(metric_values) if metric_values else 1
                    min_val = min(metric_values) if metric_values else 0
                    if max_val > min_val:
                        normalized_scores = [
                            (max_val - v) / (max_val - min_val) for v in metric_values]
                        avg_score = np.mean(normalized_scores)
                    else:
                        avg_score = 1.0  # All values equal
                else:
                    # For metrics where higher is better, just normalize
                    max_val = max(metric_values) if metric_values else 1
                    min_val = min(metric_values) if metric_values else 0
                    if max_val > min_val:
                        normalized_scores = [
                            (v - min_val) / (max_val - min_val) for v in metric_values]
                        avg_score = np.mean(normalized_scores)
                    else:
                        avg_score = 1.0  # All values equal

                scores_by_metric[metric] = {
                    'raw_mean': np.mean(metric_values),
                    'normalized': avg_score,
                    'n_values': len(metric_values)
                }
                total_scores.append(avg_score)

        if total_scores:
            method_scores[method] = {
                'overall': np.mean(total_scores),
                'by_metric': scores_by_metric,
                'n_total': len(total_scores)
            }

    # Rank methods by overall score
    ranked_methods = sorted(method_scores.items(),
                            key=lambda x: x[1]['overall'], reverse=True)

    # Print overall ranking
    print(f"\n{'Rank':<6} {'Method':<15} {'Overall Score':<15} {'Performance'}")
    print("-" * 60)
    for rank, (method, scores) in enumerate(ranked_methods, 1):
        if rank == 1:
            medal = "🥇"
        elif rank == 2:
            medal = "🥈"
        elif rank == 3:
            medal = "🥉"
        else:
            medal = f"  {rank}"

        overall = scores['overall']
        perf_emoji = "🟢" if overall >= 0.8 else "🟡" if overall >= 0.6 else "🔴"
        print(f"{medal:<6} {method:<15} {overall:<15.3f} {perf_emoji}")

    # Print detailed breakdown by metric
    if ranked_methods:
        print(
            f"\n🏆 Best overall method: {ranked_methods[0][0]}\
                (Score={ranked_methods[0][1]['overall']:.3f})")
        if len(ranked_methods) > 1:
            print(
                f"   Runner-up: {ranked_methods[1][0]}\
                    (Score={ranked_methods[1][1]['overall']:.3f})")

        # Print detailed metric breakdown for top methods
        print(f"\n📊 Detailed metric breakdown for top methods:")
        print("-" * 80)
        print(f"{'Method':<15} ", end="")
        for metric in METRIC_COLS:
            print(f"{metric:<12} ", end="")
        print()
        print("-" * 80)

        for method, scores in ranked_methods[:3]:  # Show top 3
            print(f"{method:<15} ", end="")
            for metric in METRIC_COLS:
                if metric in scores['by_metric']:
                    raw = scores['by_metric'][metric]['raw_mean']
                    norm = scores['by_metric'][metric]['normalized']
                    print(f"{raw:.3f}({norm:.2f}) ", end="")
                else:
                    print(f"{'N/A':<12} ", end="")
            print()


def witnesses_from_truth(true_newick):
    """Real siglum universe = non-placeholder labels of the TRUE tree only."""
    tt = dendropy.Tree.get(data=true_newick, schema="newick",
                           suppress_internal_node_taxa=False)
    return evaluate_reconstructed_trees.infer_witnesses(tt)


def score_pair(true_newick, recon_newick, method, witnesses):
    cfg = METHOD_CONFIG[method]
    blocks = {
        "robinson_foulds":
        evaluate_reconstructed_trees.clade_based_metrics(
            true_newick, recon_newick, witnesses=witnesses,
            rooted=cfg["rooted"]),
        "roos_sign_distance":
        evaluate_reconstructed_trees.average_sign_distance(
            true_newick, recon_newick, witnesses),
        "internal_ancestors":
        evaluate_reconstructed_trees.internal_ancestor_identification(
            true_newick, recon_newick, witnesses),
    }
    return {col: blocks[block].get(field)
            for col, block, field, _ in METRIC_SPECS}


def _score_task(task):
    """Worker: read one truth cell, score every method against it."""
    rows, miss = [], []
    true_nwk = open(task["tpath"], encoding="utf-8").read()
    witnesses = witnesses_from_truth(true_nwk)
    for method, rpath in task["recons"].items():
        if not os.path.exists(rpath):
            miss.append((task["scenario"], method))
            continue
        recon_nwk = open(rpath, encoding="utf-8").read()
        try:
            row = score_pair(true_nwk, recon_nwk, method, witnesses)
        except Exception as exc:
            warnings.warn(f"score failed {task['tpath']} / {method}: {exc}")
            miss.append((task["scenario"], method + ":error"))
            continue
        row.update(scenario=task["scenario"], cell=task["cell"],
                   method=method, replicate=task["ix"])
        rows.append(row)
    return rows, miss


def _build_tasks(methods, n_rep):
    tasks, missing = [], defaultdict(int)
    for scenario, spec in DESIGN.items():
        for params in spec["grid"]:
            label = spec["label"](params)
            for ix in range(n_rep):
                basename = sanitize(spec["fname"](params, ix))
                tpath = truth_path(scenario, basename)
                if not os.path.exists(tpath):
                    missing[(scenario, "TRUTH")] += 1
                    continue
                recons = {m: recon_path(m, scenario, basename)
                          for m in methods}
                tasks.append(dict(scenario=scenario, cell=label, ix=ix,
                                  tpath=tpath, recons=recons))
    return tasks, missing


def collect_records(methods, n_rep, jobs):
    """Walk the design grid -> one record per (scenario, cell, method, replicate)."""
    tasks, missing = _build_tasks(methods, n_rep)
    all_rows = []
    if jobs > 1:
        with mp.Pool(jobs) as pool:
            for rows, miss in pool.imap_unordered(_score_task, tasks, chunksize=8):
                all_rows.extend(rows)
                for k in miss:
                    missing[k] += 1
    else:
        for task in tasks:
            rows, miss = _score_task(task)
            all_rows.extend(rows)
            for k in miss:
                missing[k] += 1
    return pd.DataFrame.from_records(all_rows), missing


# --------------------------------------------------------------------------
# Aggregation & rendering
# --------------------------------------------------------------------------
def cell_order(scenario):
    return [DESIGN[scenario]["label"](p) for p in DESIGN[scenario]["grid"]]


def aggregate(df):
    """(scenario, cell, method) -> mean, std, n for every metric."""
    g = df.groupby(["scenario", "cell", "method"])
    mean = g[METRIC_COLS].mean()
    std = g[METRIC_COLS].std()
    n = g.size().rename("n")
    return mean, std, n


def _fmt(mean, std, kind):
    if pd.isna(mean):
        return "—"
    if kind == "count":
        return f"{mean:.1f}±{std:.1f}" if not pd.isna(std) else f"{mean:.1f}"
    if kind == "pct":
        return f"{mean:.1f}±{std:.1f}" if not pd.isna(std) else f"{mean:.1f}"
    return f"{mean:.3f}±{std:.3f}" if not pd.isna(std) else f"{mean:.3f}"


def wide_display(mean, std, n, scenario, cols=DISPLAY_COLS):
    """A formatted mean±std table for one scenario: rows = (cell, method)."""
    rows, index = [], []
    for cell in cell_order(scenario):
        for method in METHODS:
            key = (scenario, cell, method)
            if key not in mean.index:
                continue
            index.append((cell, method, int(n.loc[key])))
            rows.append([_fmt(mean.loc[key, c], std.loc[key, c], KIND[c])
                         for c in cols])
    idx = pd.MultiIndex.from_tuples(index, names=["cell", "method", "n"])
    return pd.DataFrame(rows, index=idx, columns=cols)


def numeric_wide(mean, scenario, cols=METRIC_COLS):
    """Machine-friendly mean-only wide table for one scenario."""
    rows, index = [], []
    for cell in cell_order(scenario):
        for method in METHODS:
            key = (scenario, cell, method)
            if key in mean.index:
                index.append((scenario, cell, method))
                rows.append([mean.loc[key, c] for c in cols])
    idx = pd.MultiIndex.from_tuples(
        index, names=["scenario", "cell", "method"])
    return pd.DataFrame(rows, index=idx, columns=cols)


def to_booktabs(disp_df, scenario, caption):
    """LaTeX booktabs body from a formatted display table."""
    def esc(s):
        return str(s).replace("_", r"\_").replace("%", r"\%")

    cols = list(disp_df.columns)
    head = "cell & method & n & " + " & ".join(esc(c) for c in cols) + r" \\"
    lines = [r"\begin{table}[t]", r"\centering",
             r"\caption{" + esc(caption) + "}",
             r"\label{tab:" + scenario.replace("_", "") + "}",
             r"\begin{tabular}{ll r " + "r" * len(cols) + "}",
             r"\toprule", head, r"\midrule"]
    prev_cell = None
    for (cell, method, n), row in disp_df.iterrows():
        cell_txt = esc(cell) if cell != prev_cell else ""
        prev_cell = cell
        vals = " & ".join(str(v).replace("±", r"\,$\pm$\,").replace("%", r"\%")
                          for v in row)
        lines.append(f"{cell_txt} & {esc(method)} & {n} & {vals} " + r"\\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    return "\n".join(lines)


def main():
    global TRUTH_ROOT, RECON_ROOT, METHODS
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--truth-root", default=TRUTH_ROOT)
    ap.add_argument("--recon-root", default=RECON_ROOT)
    ap.add_argument("--out-dir", default=OUT_DIR)
    ap.add_argument("--replicates", type=int, default=N_REPLICATES)
    ap.add_argument("--methods", default=",".join(METHODS))
    ap.add_argument("--jobs", type=int, default=1,
                    help="parallel workers (each pair is independent)")
    a = ap.parse_args()

    TRUTH_ROOT, RECON_ROOT = a.truth_root, a.recon_root
    METHODS = a.methods.split(",")
    os.makedirs(a.out_dir, exist_ok=True)

    print_pretty_header("🔬 STEMMATIC RECONSTRUCTION BENCHMARK")
    print(f"📁 Truth root: {TRUTH_ROOT}")
    print(f"📁 Recon root: {RECON_ROOT}")
    print(f"🔢 Replicates: {a.replicates}")
    print(f"⚙️  Methods: {', '.join(METHODS)}")

    df, missing = collect_records(METHODS, a.replicates, a.jobs)
    if df.empty:
        print("❌ No records scored. Check truth/recon paths and naming.")
        for k, v in sorted(missing.items()):
            print(f"  missing {k}: {v}")
        return

    # long tidy table (everything downstream can be built from this)
    long_path = os.path.join(a.out_dir, "metrics_long.csv")
    df.to_csv(long_path, index=False)

    mean, std, n = aggregate(df)

    # combined machine-friendly wide mean table
    num = pd.concat([numeric_wide(mean, s) for s in DESIGN
                     if any((s, c, m) in mean.index
                            for c in cell_order(s) for m in METHODS)])
    num.round(4).to_csv(os.path.join(a.out_dir, "metrics_mean_wide.csv"))

    # per-scenario eyeball tables
    md_all = ["# Benchmark results (mean ± std over replicates)\n"]

    # Store method rankings for overall summary
    method_rankings_by_scenario = {}

    for scenario in DESIGN:
        disp = wide_display(mean, std, n, scenario)
        if disp.empty:
            continue

        # Print pretty scenario header
        print_pretty_header(f"📊 {scenario.upper()}")

        # Print the standard markdown table (keep for compatibility)
        title = f"## {scenario}"
        print("\n" + title)
        print(disp.to_markdown())

        # Print enhanced ranking table
        print_ranking_table(disp, scenario)

        # Collect rankings for overall summary
        rankings = {}
        for cell in disp.index.get_level_values('cell').unique():
            cell_data = disp.xs(cell, level='cell')
            if 'clade_F1' in cell_data.columns:
                sorted_data = cell_data.sort_values(
                    'clade_F1', ascending=False)
                rankings[cell] = [(idx[0], float(row['clade_F1'].split('±')[0]))
                                  for idx, row in sorted_data.iterrows()]
        method_rankings_by_scenario[scenario] = rankings

        md_all += [title, "", disp.to_markdown(), ""]
        with open(os.path.join(a.out_dir, f"{scenario}.tex"), "w") as fh:
            fh.write(to_booktabs(disp, scenario,
                                 f"Reconstruction accuracy, {scenario}."))

    # Print overall summary with ALL metrics
    print_overall_summary(mean)

    # coverage report
    if missing:
        print("\n⚠️  Coverage Gaps")
        print("-" * 40)
        for k, v in sorted(missing.items()):
            print(f"  {k[0]} / {k[1]}: {v} missing")

    print(f"\n✅ Done! Results written to:")
    print(f"  📄 {long_path}")
    print(f"  📄 metrics_mean_wide.csv")
    print(f"  📄 results.md")
    print(f"  📄 per-scenario .tex files → {a.out_dir}/")


if __name__ == "__main__":
    main()
