
```
#!/bin/bash
#SBATCH -J quip_decompress
#SBATCH -p day
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH -t 24:00:00
#SBATCH -o quip_%A_%a.out
#SBATCH -e quip_%A_%a.err

module --force purge
module load StdEnv
module load GCCcore/10.2.0
module load Quip/20171217-GCCcore-10.2.0

BASE_DIR="/gpfs/gibbs/pi/guo/mg2684/ATXN2_mouse"

SAMPLE_DIRS=($(find "$BASE_DIR" -type d -name Unaligned | sort))
NUM=${#SAMPLE_DIRS[@]}

DIR="${SAMPLE_DIRS[$SLURM_ARRAY_TASK_ID]}"

echo "Decompressing in: $DIR"

for qp in "$DIR"/*.qp; do
    out="${qp%.qp}"
    echo "quip -d $qp"
    quip -d "$qp"
done

touch "$DIR/.DECOMPRESS_DONE"
```
