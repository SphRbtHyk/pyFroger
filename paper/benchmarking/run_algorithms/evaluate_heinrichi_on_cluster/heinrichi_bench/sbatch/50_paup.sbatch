#!/usr/bin/env bash
#SBATCH --job-name=hein_paup
#SBATCH --output=logs/05_paup.%j.out
#SBATCH --error=logs/05_paup.%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

# PAUP* — the reference run.
#
# This is the ONLY job that scores the matrix as it actually is: PAUP's
# standard datatype takes the 17 states directly, so there is no recoding
# (MrBayes caps at 10, PHYLIP's pars at 8) and no distance approximation
# (FastME). It is also the program Roos & Heikkilae used, so this is the run
# that tells you whether the published 74.4 is reproducible at all.
#
# No container: PAUP* is a static binary. Put it at bin/paup (chmod +x) or
# set PAUP_BIN.

export PAUP_BIN=/home/roberths/dev/pyfroger/paper/benchmarking/run_algorithms/paup4a168_centos64

source /home/roberths/dev/pyfroger/paper/benchmarking/run_algorithms/evaluate_heinrichi_on_cluster/heinrichi_bench/config.sh


if [[ ! -x "$PAUP_BIN" ]]; then
    echo "error: $PAUP_BIN is not an executable." >&2
    echo "       scp the PAUP* Linux binary there and chmod +x it." >&2
    exit 1
fi

DIR="$OUT/paup"
mkdir -p "$DIR"
cp "$NEXUS" "$DIR/heinrichi.nex"
cd "$DIR"

# Note on collapse: PAUP's default is collapse=yes, which is exactly the
# behaviour the Roos metric rewards (polytomies produce the ties that the true
# stemma genuinely has). We set it explicitly rather than relying on a default
# that could differ between builds.
cat > run.nex <<'PAUP'
#NEXUS
begin paup;
    set autoclose=yes warntree=no warnreset=no warnroot=no errorstop=no;
    log file=paup.log replace=yes;

    execute heinrichi.nex;

    set criterion=parsimony;
    pset collapse=minbrlen;          /* = PAUP's collapse=yes default */

    /* A real search: 100 random-addition replicates, TBR, keep all MPTs.
       This is the search that Biopython's NNI-from-one-NJ-start was not. */
    hsearch addseq=random nreps=100 swap=tbr multrees=yes;

    savetrees file=mpts.tre format=newick brlens=yes replace=yes;
    contree all / strict=yes treefile=strict.tre format=newick replace=yes;
    contree all / majrule=yes percent=50 treefile=majrule.tre
                  format=newick replace=yes;

    log stop;
    quit;
end;
PAUP

"$PAUP_BIN" -n run.nex

echo
echo "--- tree length / number of MPTs (check against the paper) ---"
grep -iE "Length|tree(s) retained|Number of trees" paup.log | head -20 || true

echo
echo "paup done -> $DIR/mpts.tre (all MPTs), $DIR/strict.tre, $DIR/majrule.tre"
