#!/usr/bin/env bash
# Build the .sif images ONCE, on a LOGIN node (compute nodes generally have no
# outbound network, and three concurrent jobs converting the same image race
# each other in the cache).
#
#   ./pull_images.sh
#
set -euo pipefail
source "$(dirname "$0")/config.sh"

# Keep the (large) conversion cache off $HOME, which is often quota-limited.
export SINGULARITY_CACHEDIR="${SINGULARITY_CACHEDIR:-$SIF/.cache}"
export APPTAINER_CACHEDIR="$SINGULARITY_CACHEDIR"
export SINGULARITY_TMPDIR="${SINGULARITY_TMPDIR:-${TMPDIR:-/tmp}}"
export APPTAINER_TMPDIR="$SINGULARITY_TMPDIR"
mkdir -p "$SINGULARITY_CACHEDIR"

ENG="${ENGINE}"
[[ "$ENG" == "docker" || "$ENG" == "none" ]] && { echo "ENGINE=$ENG: nothing to pull"; exit 0; }

pull() {
    local sif="$1" uri="$2"
    [[ "$uri" != docker://* && "$uri" != library://* ]] && uri="docker://$uri"
    if [[ -f "$sif" ]]; then
        echo "have    $(basename "$sif")"
    else
        echo "pulling $(basename "$sif")  <- $uri"
        "$ENG" pull "$sif" "$uri"
    fi
}

pull "$SIF_IQTREE"  "$IMG_IQTREE"
pull "$SIF_MRBAYES" "$IMG_MRBAYES"
pull "$SIF_FASTME"  "$IMG_FASTME"

echo
echo "--- verifying the entrypoints ---"
"$ENG" exec "$SIF_IQTREE"  iqtree2 -h 2>&1 | head -1 || echo "!! iqtree2 not on PATH in the image"
"$ENG" exec "$SIF_MRBAYES" mb -v    2>&1 | head -3 || echo "!! mb not on PATH in the image"
"$ENG" exec "$SIF_FASTME"  fastme --version 2>&1 | head -1 || echo "!! fastme not on PATH in the image"

echo
echo "images ready in $SIF"
