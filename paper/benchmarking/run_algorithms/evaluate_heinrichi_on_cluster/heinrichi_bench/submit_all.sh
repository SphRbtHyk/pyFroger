#!/usr/bin/env bash
# Submit the whole chain with dependencies:
#   prepare -> {iqtree, mrbayes, fastme} in parallel -> score
set -euo pipefail
cd "$(dirname "$0")"
source /home/roberths/dev/pyfroger/paper/benchmarking/run_algorithms/evaluate_heinrichi_on_cluster/heinrichi_bench/config.sh
mkdir -p logs

echo "SIF MRBAYES is here" $SIF_MRBAYES

if [[ "$ENGINE" == "singularity" || "$ENGINE" == "apptainer" ]]; then
    for s in "$SIF_IQTREE" "$SIF_MRBAYES" "$SIF_FASTME"; do
        if [[ ! -f "$s" ]]; then
            echo "error: $(basename "$s") is missing." >&2
            echo "       Run ./pull_images.sh on a LOGIN node before submitting." >&2
            exit 1
        fi
    done
fi

for f in "$NEXUS" "$TRUE_GRAPH"; do
    [[ -f "$f" ]] || { echo "error: missing input $f" >&2; exit 1; }
done

PREP=$(sbatch --parsable sbatch/00_prepare.sbatch)
echo "prepare  : $PREP"

PP=$(sbatch --parsable sbatch/50_paup.sbatch)
echo "paup     : $PP"


IQ=$(sbatch --parsable --dependency=afterok:$PREP sbatch/10_iqtree.sbatch)
MB=$(sbatch --parsable --dependency=afterok:$PREP sbatch/20_mrbayes.sbatch)
FM=$(sbatch --parsable --dependency=afterok:$PREP sbatch/30_fastme.sbatch)
echo "iqtree   : $IQ"
echo "mrbayes  : $MB"
echo "fastme   : $FM"

# afterany, not afterok: if one tool fails we still want the others' scores.
SC=$(sbatch --parsable --dependency=afterany:$PP:$IQ:$MB:$FM sbatch/40_score.sbatch)
echo "score    : $SC"
