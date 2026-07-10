# RUN THE PIPELINE (ARRAY LIMIT IS SET TO 1000)

sbatch --export=ALL,BASE=0 --array=0-999%20 bench.sbatch   # files 0–999
sbatch --export=ALL,BASE=1000 --array=0-999%20 bench.sbatch   # files 1000–1999
sbatch --export=ALL,BASE=2000 --array=0-899%20 bench.sbatch   # files 2000–2899