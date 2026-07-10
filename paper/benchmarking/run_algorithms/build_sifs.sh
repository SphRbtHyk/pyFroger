#!/usr/bin/env bash
# Pull the three evolbioinfo images to local SIF files.
set -euo pipefail

export SINGULARITY_TMPDIR="/scratch/roberths/sing_tmp"; mkdir -p "$SINGULARITY_TMPDIR"

module load singularity/4.2.2

singularity pull --force $SINGULARITY_TMPDIR/iqtree.sif  docker://evolbioinfo/iqtree:v3.1.1
singularity pull --force $SINGULARITY_TMPDIR/mrbayes.sif docker://pegi3s/mrbayes:3.2.7
singularity pull --force $SINGULARITY_TMPDIR/fastme.sif  docker://evolbioinfo/fastme:v2.1.6.4

echo "built: $(ls -1 $SINGULARITY_TMPDIR/*.sif | tr '\n' ' ')"
