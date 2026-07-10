#!/usr/bin/env bash
#SBATCH --job-name=stemma-bench
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#
# Prereqs:
#   ./build_sifs.sh
#   ls /home/roberths/dev/pyfroger/paper/datasets/synthetic/generated_trees/*/*.fasta > files.txt
# Submit (array range derived from the manifest length, %20 = max concurrent):
#   sbatch --array=0-$(( $(wc -l < files.txt) - 1 ))%20 bench.sbatch

set -euo pipefail

# --- environment ------------------------------------------------------------
module load singularity/4.2.2
# Activate the python environment
source /home/roberths/dev/pyfroger/venv/bin/activate

SING_DIR="/scratch/roberths/sing_tmp"
DATADIR=/home/roberths/dev/pyfroger/paper/datasets/synthetic/generated_trees
WORKDIR=/home/roberths/dev/pyfroger/paper/benchmarking/run_algorithms
OUT=/scratch/roberths/stemmabench
cd $WORKDIR


MANIFEST=$WORKDIR/files.txt
mkdir -p logs

BASE=${BASE:-0}
IDX=$(( SLURM_ARRAY_TASK_ID + BASE ))
FASTA=$(sed -n "$((IDX + 1))p" "$MANIFEST")
[ -n "${FASTA:-}" ] || { echo "no file at index $IDX"; exit 0; }

mkdir -p "$OUT/metrics"

echo "Working on $FASTA"

python3 run_bench.py "$FASTA" \
    --out "$OUT" \
    --methods "mrbayes" \
    --runtime singularity \
    --iqtree-image  $SING_DIR/iqtree.sif \
    --fastme-image  $SING_DIR/fastme.sif \
    --mrbayes-image $SING_DIR/mrbayes.sif \
    --iqtree-threads "$SLURM_CPUS_PER_TASK" \
    --iqtree-model MK+ASC \
    --outgroup w1 \
    --ngen 400000 --samplefreq 100 \
    --metrics "$OUT/metrics/${IDX}.csv"
