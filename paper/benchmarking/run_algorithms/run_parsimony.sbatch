#!/usr/bin/env bash
#SBATCH --job-name=stemma-bench
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00


set -euo pipefail

# Activate the python environment
source /home/roberths/dev/pyfroger/venv/bin/activate

DATADIR=/home/roberths/dev/pyfroger/paper/datasets/synthetic/generated_trees
WORKDIR=/home/roberths/dev/pyfroger/paper/benchmarking/run_algorithms
OUTDIR=/scratch/roberths/stemmabench/trees/parsimony
cd $WORKDIR

MANIFEST=$WORKDIR/files.txt
mkdir -p logs

python run_parsimony.py