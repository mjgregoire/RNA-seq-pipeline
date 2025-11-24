
```
#!/bin/bash
#SBATCH -J star_quip_align
#SBATCH -p day
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=60G
#SBATCH -t 24:00:00
#SBATCH -o star_%A_%a.out
#SBATCH -e star_%A_%a.err
#SBATCH --mail-type ALL

#############################################
# USER SETTINGS
#############################################

BASE_DIR="/gpfs/gibbs/pi/guo/mg2684/ATXN2_mouse"
GENOME_DIR="/gpfs/gibbs/pi/guo/mg2684/reference/STAR_index_GRCh38"    
OUT_DIR="/gpfs/gibbs/pi/guo/mg2684/ATXN2_mouse/star_out"

mkdir -p "$OUT_DIR"

#############################################
# MODULES
#############################################

module purge
module load StdEnv
module load Quip/20171217-GCCcore-10.2.0
module load STAR/2.7.10a-GCC-12.2.0
module load SAMtools/1.17-GCC-12.2.0

echo "Loaded modules:"
module list

#############################################
# DISCOVER SAMPLE FOLDERS
#############################################

SAMPLE_DIRS=($(find "$BASE_DIR" -type d -name Unaligned | sort))
NUM_SAMPLES=${#SAMPLE_DIRS[@]}

echo "Detected $NUM_SAMPLES samples."

if [[ -z "$SLURM_ARRAY_TASK_ID" ]]; then
    echo "Run using: sbatch --array=0-$((NUM_SAMPLES-1)) run_star_with_quip_array.sh"
    exit 1
fi

SAMPLE_UNALIGNED=${SAMPLE_DIRS[$SLURM_ARRAY_TASK_ID]}
SAMPLE_NAME=$(basename "$(dirname "$SAMPLE_UNALIGNED")")

echo "Processing sample: $SAMPLE_NAME"
echo "Unaligned folder: $SAMPLE_UNALIGNED"

#############################################
# STEP 1 — DECOMPRESS any .qp files with Quip
#############################################

QP_FILES=($(ls "$SAMPLE_UNALIGNED"/*.fastq.qp 2>/dev/null))

if (( ${#QP_FILES[@]} > 0 )); then
    echo "Decompressing ${#QP_FILES[@]} .qp files..."

    for qp in "${QP_FILES[@]}"; do
        quip -d "$qp"
    done
else
    echo "No .qp files found. Using existing FASTQs."
fi

#############################################
# STEP 2 — GET FASTQ FILES
#############################################

FASTQS=($(ls "$SAMPLE_UNALIGNED"/*.fastq 2>/dev/null))

if (( ${#FASTQS[@]} == 0 )); then
    echo "ERROR: No FASTQ files found after decompression."
    exit 1
fi

echo "Found FASTQs:"
printf '%s\n' "${FASTQS[@]}"

# Detect paired vs single
if (( ${#FASTQS[@]} == 1 )); then
    echo "Detected single-end read"
    FASTQ_CMD="--readFilesIn ${FASTQS[0]}"
elif (( ${#FASTQS[@]} == 2 )); then
    echo "Detected paired-end reads"
    FASTQ_CMD="--readFilesIn ${FASTQS[0]} ${FASTQS[1]}"
else
    echo "ERROR: Unexpected number of FASTQs for $SAMPLE_NAME"
    exit 1
fi

#############################################
# STEP 3 — RUN STAR
#############################################

echo "Running STAR for $SAMPLE_NAME..."

STAR \
    --runThreadN 8 \
    --genomeDir "$GENOME_DIR" \
    $FASTQ_CMD \
    --readFilesCommand cat \
    --outFileNamePrefix "$OUT_DIR/${SAMPLE_NAME}_" \
    --outSAMtype BAM SortedByCoordinate

#############################################
# STEP 4 — BAM INDEX
#############################################

samtools index "$OUT_DIR/${SAMPLE_NAME}_Aligned.sortedByCoord.out.bam"

echo "Completed sample: $SAMPLE_NAME"
```
