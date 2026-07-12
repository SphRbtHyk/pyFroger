#!/usr/bin/env bash
# Shared configuration for the Heinrichi benchmark.
# Source this from every sbatch script.

set -euo pipefail

# ── Paths ────────────────────────────────────────────────────────────────
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Inputs you must provide:
NEXUS="${NEXUS:-$ROOT/data/heinrichi_bench/heinrichi.nex}"              # the 37x1208 matrix
TRUE_GRAPH="${TRUE_GRAPH:-$ROOT/data/heinrichi_bench/heinrichi.graph}"  # the 65-node adjacency file

WORK="${WORK:-$ROOT/work}"          # prepared inputs
OUT="${OUT:-$ROOT/results}"         # tool outputs
LOGS="${LOGS:-$ROOT/logs}"
SIF="${SIF:-$ROOT/../sif}"             # pre-pulled container images

mkdir -p "$WORK" "$OUT" "$LOGS" "$SIF"

# ── Containers ───────────────────────────────────────────────────────────
# singularity on the cluster; docker on a workstation; none = tools on $PATH.
ENGINE="${ENGINE:-singularity}"

# Point these at whatever images the simulated benchmark already uses.
IMG_IQTREE="${IMG_IQTREE:-docker://evolbioinfo/iqtree:v3.1.1}"
IMG_MRBAYES="${IMG_MRBAYES:-docker://pegi3s/mrbayes:3.2.7}"
IMG_FASTME="${IMG_FASTME:-docker://evolbioinfo/fastme:v2.1.6.4}"

# Local .sif filenames, produced by ./pull_images.sh on the LOGIN node.
SIF_IQTREE="$SIF/iqtree.sif"
SIF_MRBAYES="$SIF/mrbayes.sif"
SIF_FASTME="$SIF/fastme.sif"

# run_container <sif-or-image> <workdir-on-host> <command...>
run_container() {
    local img="$1"; shift
    local wd="$1"; shift
    case "$ENGINE" in
        singularity|apptainer)
            if [[ ! -f "$img" ]]; then
                echo "error: $img not found." >&2
                echo "       Run ./pull_images.sh on a LOGIN node first —" >&2
                echo "       compute nodes usually have no outbound network." >&2
                exit 1
            fi
            # --cleanenv  : do not leak the host environment (conda, PYTHONPATH)
            # --no-home   : do not auto-mount $HOME over the container
            # --bind/--pwd: the job directory is the container's /data
            "$ENGINE" exec \
                --cleanenv \
                --no-home \
                --bind "$wd":/data \
                --pwd /data \
                "$img" "$@"
            ;;
        docker)
            docker run --rm \
                -u "$(id -u):$(id -g)" \
                -v "$wd":/data -w /data \
                "${img#docker://}" "$@"
            ;;
        none)
            ( cd "$wd" && "$@" )
            ;;
        *)
            echo "unknown ENGINE=$ENGINE" >&2; exit 1
            ;;
    esac
}

# Pick the right image handle for the active engine.
image_for() {
    case "$ENGINE" in
        singularity|apptainer)
            case "$1" in
                iqtree)  echo "$SIF_IQTREE" ;;
                mrbayes) echo "$SIF_MRBAYES" ;;
                fastme)  echo "$SIF_FASTME" ;;
            esac
            ;;
        *)
            case "$1" in
                iqtree)  echo "$IMG_IQTREE" ;;
                mrbayes) echo "$IMG_MRBAYES" ;;
                fastme)  echo "$IMG_FASTME" ;;
            esac
            ;;
    esac
}

# ── Analysis parameters ──────────────────────────────────────────────────
MB_NGEN="${MB_NGEN:-2000000}"
MB_NRUNS="${MB_NRUNS:-2}"
MB_NCHAINS="${MB_NCHAINS:-4}"
MB_SAMPLEFREQ="${MB_SAMPLEFREQ:-1000}"
MB_BURNINFRAC="${MB_BURNINFRAC:-0.25}"

IQ_MODEL="${IQ_MODEL:-MK+FQ+G4}"   # set to MFP to let ModelFinder choose
IQ_BB="${IQ_BB:-1000}"             # ultrafast bootstrap replicates

# Branches shorter than this are collapsed to polytomies before scoring.
# PAUP*'s default is collapse=yes; matching it matters for the Roos metric.
COLLAPSE_EPS="${COLLAPSE_EPS:-1e-6}"

# Witnesses in the matrix but absent from the true stemma: they take part in
# the search, not in the scoring (Roos et al. build on 37, score on 34).
# score.py handles this by intersecting the label sets; noted here for the record.
EXTRA_TAXA="Da I J"

export ROOT NEXUS TRUE_GRAPH WORK OUT LOGS SIF ENGINE
export IMG_IQTREE IMG_MRBAYES IMG_FASTME
export SIF_IQTREE SIF_MRBAYES SIF_FASTME
export MB_NGEN MB_NRUNS MB_NCHAINS MB_SAMPLEFREQ MB_BURNINFRAC
export IQ_MODEL IQ_BB COLLAPSE_EPS
