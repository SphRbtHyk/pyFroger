#!/usr/bin/env bash
#SBATCH --job-name=hein_iqtree
#SBATCH --output=logs/10_iqtree.%j.out
#SBATCH --error=logs/10_iqtree.%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

source /home/roberths/dev/pyfroger/paper/benchmarking/run_algorithms/evaluate_heinrichi_on_cluster/heinrichi_bench/config.sh

DIR="$OUT/iqtree"
mkdir -p "$DIR"
cp "$WORK/heinrichi.morph.phy" "$DIR/"

# -st MORPH : multistate morphology, Mk model (handles the 17 states).
# No +ASC   : the matrix keeps its constant characters, so an
#             ascertainment-bias correction would be wrong here.
run_container "$(image_for iqtree)" "$DIR" \
    iqtree2 \
        -s heinrichi.morph.phy \
        -st MORPH \
        -m "$IQ_MODEL" \
        -B "$IQ_BB" \
        -T "${SLURM_CPUS_PER_TASK:-4}" \
        -pre heinrichi \
        -redo

echo "iqtree done -> $DIR/heinrichi.treefile"
