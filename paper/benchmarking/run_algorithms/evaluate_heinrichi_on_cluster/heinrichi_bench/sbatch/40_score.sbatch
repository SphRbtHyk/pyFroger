#!/usr/bin/env bash
#SBATCH --job-name=hein_score
#SBATCH --output=logs/40_score.%j.out
#SBATCH --error=logs/40_score.%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

source /home/roberths/dev/pyfroger/paper/benchmarking/run_algorithms/evaluate_heinrichi_on_cluster/heinrichi_bench/config.sh
source /home/roberths/dev/pyfroger/venv/bin/activate

python "$ROOT/scripts/score.py" \
    --true "$TRUE_GRAPH" \
    --collapse-eps "$COLLAPSE_EPS" \
    --tree "iqtree=$OUT/iqtree/heinrichi.treefile" \
    --tree "mrbayes=$OUT/mrbayes/heinrichi.mb.nex.con.tre" \
    --tree "fastme=$OUT/fastme/heinrichi.nwk" \
    --tree "paup_strict=$OUT/paup/strict.tre" \
    --tree "paup_majrule=$OUT/paup/majrule.tre" \
    --out "$OUT/scores.tsv"

echo
echo "scores written to $OUT/scores.tsv"
