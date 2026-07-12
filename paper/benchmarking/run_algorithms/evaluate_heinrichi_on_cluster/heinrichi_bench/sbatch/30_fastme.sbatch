#!/usr/bin/env bash
#SBATCH --job-name=hein_fastme
#SBATCH --output=logs/30_fastme.%j.out
#SBATCH --error=logs/30_fastme.%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G

source /home/roberths/dev/pyfroger/paper/benchmarking/run_algorithms/evaluate_heinrichi_on_cluster/heinrichi_bench/config.sh

DIR="$OUT/fastme"
mkdir -p "$DIR"
cp "$WORK/heinrichi.dist" "$DIR/"

# -i    input is a PHYLIP distance matrix. FastME cannot derive distances from
#       multistate characters, so prepare.py computed them with
#       pairwise-complete deletion (missing data contributes neither a match
#       nor a mismatch).
# -m B  balanced minimum evolution
# -n B  BME nearest-neighbour interchange
# -s    SPR postprocessing
run_container "$(image_for fastme)" "$DIR" \
    fastme \
        -i heinrichi.dist \
        -o heinrichi.nwk \
        -I heinrichi.stat \
        -m B -n B -s

echo "fastme done -> $DIR/heinrichi.nwk"
