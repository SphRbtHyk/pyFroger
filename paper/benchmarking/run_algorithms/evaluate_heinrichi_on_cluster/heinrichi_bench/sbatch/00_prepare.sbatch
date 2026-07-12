#!/usr/bin/env bash
#SBATCH --job-name=hein_prep
#SBATCH --output=logs/00_prepare.%j.out
#SBATCH --error=logs/00_prepare.%j.err
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

source "$(dirname "$0")/../config.sh"
source /home/roberths/dev/pyfroger/venv/bin/activate

python3 "$ROOT/scripts/prepare.py" "$NEXUS" \
    --outdir "$WORK" \
    --many-states binarise \
    --ngen "$MB_NGEN" \
    --nruns "$MB_NRUNS" \
    --nchains "$MB_NCHAINS" \
    --samplefreq "$MB_SAMPLEFREQ" \
    --burninfrac "$MB_BURNINFRAC"

echo "prepare done"
