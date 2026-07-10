source /home/roberths/dev/pyfroger/venv/bin/activate

DATADIR=/home/roberths/dev/pyfroger/paper/datasets/synthetic/generated_trees
WORKDIR=/home/roberths/dev/pyfroger/paper/benchmarking/run_algorithms
OUTDIR=/scratch/roberths/stemmabench/trees/parsimony
cd $WORKDIR

# Run the jobs
sbatch --array=0-999%20 --export=ALL,BASE=0,MANIFEST=$PWD/manifest_pars.txt,DATADIR=$DATADIR,OUTDIR=$OUTDIR run_parsimony.sbatch
sbatch --array=0-999%20 --export=ALL,BASE=1000,MANIFEST=$PWD/manifest_pars.txt,DATADIR=$DATADIR,OUTDIR=$OUTDIR run_parsimony.sbatch
sbatch --array=0-899%20 --export=ALL,BASE=2000,MANIFEST=$PWD/manifest_pars.txt,DATADIR=$DATADIR,OUTDIR=$OUTDIR run_parsimony.sbatch
