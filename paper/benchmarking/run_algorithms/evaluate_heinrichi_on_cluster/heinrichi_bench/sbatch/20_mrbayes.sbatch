#!/usr/bin/env bash
#SBATCH --job-name=hein_mb
#SBATCH --output=logs/20_mrbayes.%j.out
#SBATCH --error=logs/20_mrbayes.%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

source /home/roberths/dev/pyfroger/paper/benchmarking/run_algorithms/evaluate_heinrichi_on_cluster/heinrichi_bench/config.sh

DIR="$OUT/mrbayes"
mkdir -p "$DIR"
cp "$WORK/heinrichi.mb.nex" "$DIR/"

# The mrbayes block (mcmc / sumt / sump) is inside the nexus file, written by
# prepare.py, so the file is simply executed. Note that prepare.py had to
# recode the matrix: MrBayes standard datatype allows at most 10 states and
# the Heinrichi matrix has 17.
run_container "$(image_for mrbayes)" "$DIR" \
    mb heinrichi.mb.nex

echo "--- convergence: ASDSF must be < 0.01, otherwise do NOT use the tree ---"
grep -h "Average standard deviation of split frequencies" \
    "$LOGS"/20_mrbayes."$SLURM_JOB_ID".out 2>/dev/null | tail -3 || true

echo "mrbayes done -> $DIR/heinrichi.mb.nex.con.tre"
