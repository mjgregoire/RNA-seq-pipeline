
```
#!/bin/bash
#SBATCH -J quip_to_gzip
#SBATCH -p day
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH -t 24:00:00
#SBATCH -o quip_to_gzip_%j.out
#SBATCH -e quip_to_gzip_%j.err

set -euo pipefail

# Load the modules needed for Quip
module --force purge
module load StdEnv
module load GCCcore/10.2.0
module load Quip/20171217-GCCcore-10.2.0

BASE="/gpfs/gibbs/pi/guo/mg2684/ATXN2_mouse"

echo "========= Starting Quip → Gzip conversion ========="
echo "Time: $(date)"
echo "Base directory: $BASE"
echo "---------------------------------------------------"

# Find all .fastq.qp files
mapfile -t QP_FILES < <(find "$BASE" -type f -name "*.fastq.qp" | sort)

echo "Found ${#QP_FILES[@]} .fastq.qp files to convert."

for qp in "${QP_FILES[@]}"; do
    out="${qp%.fastq.qp}.fastq.gz"
    echo "Processing:"
    echo "  QP : $qp"
    echo "  OUT: $out"

    # Decompress and gzip in one step
    quip -d "$qp" | gzip > "$out"

    # Optional: remove the old qp file only after success
    # rm "$qp"

    echo "  ✔ Done"
    echo "---------------------------------------------------"
done

echo "========= Finished at $(date) ========="
```
`sbatch quip_to_gzip.sh`


```
#!/bin/bash
#SBATCH -J star_align
#SBATCH -p day
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=60G
#SBATCH -t 24:00:00
#SBATCH -o star_%A_%a.out
#SBATCH -e star_%A_%a.err
#SBATCH --mail-type=ALL

module --force purge
module load StdEnv
module load STAR/2.7.11a-GCC-12.2.0
module load SAMtools/1.21-GCC-12.2.0

BASE_DIR="/gpfs/gibbs/pi/guo/mg2684/ATXN2_mouse"
GENOME_DIR="/gpfs/gibbs/pi/guo/mg2684/reference/STAR_index_GRCh38"
OUT_DIR="/gpfs/gibbs/pi/guo/mg2684/ATXN2_mouse/star_out"
mkdir -p "$OUT_DIR"

# Identify Unaligned dirs
SAMPLE_DIRS=($(find "$BASE_DIR" -type d -name Unaligned | sort))
NUM_SAMPLES=${#SAMPLE_DIRS[@]}

echo "Detected $NUM_SAMPLES sample folders"

if [[ -z "$SLURM_ARRAY_TASK_ID" ]]; then
    echo "Run as:"
    echo "  sbatch --array=0-$((NUM_SAMPLES-1)) $0"
    exit 1
fi

UNALIGNED="${SAMPLE_DIRS[$SLURM_ARRAY_TASK_ID]}"
SAMPLE_NAME=$(basename "$(dirname "$UNALIGNED")")

# Only use fastq.gz files — ignore plain FASTQ and QP files
FASTQS=($(find "$UNALIGNED" -maxdepth 1 -type f -name "*.fastq.gz" | sort))

if (( ${#FASTQS[@]} == 0 )); then
    echo "ERROR: No .fastq.gz files found in $UNALIGNED"
    exit 1
elif (( ${#FASTQS[@]} == 1 )); then
    FASTQ_CMD="--readFilesIn ${FASTQS[0]}"
elif (( ${#FASTQS[@]} == 2 )); then
    FASTQ_CMD="--readFilesIn ${FASTQS[0]} ${FASTQS[1]}"
else
    echo "ERROR: Expected 1 or 2 FASTQs, found ${#FASTQS[@]} in $UNALIGNED"
    exit 1
fi

# Required for gzipped input
READ_FILES_COMMAND="--readFilesCommand zcat"

echo "Running STAR for sample: $SAMPLE_NAME"
echo "FASTQs:"
printf "  %s\n" "${FASTQS[@]}"

STAR \
    --runThreadN 8 \
    --genomeDir "$GENOME_DIR" \
    $FASTQ_CMD \
    $READ_FILES_COMMAND \
    --outFileNamePrefix "$OUT_DIR/${SAMPLE_NAME}_" \
    --outSAMtype BAM SortedByCoordinate \
    --limitBAMsortRAM 55000000000

# Index BAM
samtools index "$OUT_DIR/${SAMPLE_NAME}_Aligned.sortedByCoord.out.bam"
```
