
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

############################
### MODULES
############################
module --force purge
module load StdEnv
module load STAR/2.7.11a-GCC-12.2.0
module load SAMtools/1.21-GCC-12.2.0

############################
### PATHS
############################
BASE_DIR="/gpfs/gibbs/pi/guo/mg2684/ATXN2_mouse"
GENOME_DIR="/gpfs/gibbs/pi/guo/mg2684/reference/STAR_index_GRCh38"
OUT_DIR="$BASE_DIR/star_out"
TEMP_BASE="$OUT_DIR/temp_concat"
mkdir -p "$OUT_DIR"
mkdir -p "$TEMP_BASE"

############################
### SAMPLE LIST
############################
mapfile -t SAMPLE_DIRS < <(find "$BASE_DIR" -type d -name Unaligned | sort)
NUM_SAMPLES=${#SAMPLE_DIRS[@]}

echo "Detected $NUM_SAMPLES sample directories."

if [[ -z "$SLURM_ARRAY_TASK_ID" ]]; then
    echo "Run with: sbatch --array=0-$((NUM_SAMPLES-1)) $0"
    exit 1
fi

UNALIGNED="${SAMPLE_DIRS[$SLURM_ARRAY_TASK_ID]}"
SAMPLE_NAME=$(basename "$(dirname "$UNALIGNED")")

echo "Processing sample: $SAMPLE_NAME"
echo "Unaligned folder: $UNALIGNED"

############################
### FIND ALL LANES (fastq)
############################
R1s=($(find "$UNALIGNED" -maxdepth 1 -type f -name "*_R1_*.fastq" | sort))
R2s=($(find "$UNALIGNED" -maxdepth 1 -type f -name "*_R2_*.fastq" | sort))

if (( ${#R1s[@]} == 0 )); then
    echo "ERROR: No R1 files found — cannot align."
    exit 1
fi

if (( ${#R1s[@]} != ${#R2s[@]} )); then
    echo "ERROR: Unequal number of R1 and R2 files."
    printf "R1 files: %s\n" "${R1s[@]}"
    printf "R2 files: %s\n" "${R2s[@]}"
    exit 1
fi

echo "Found ${#R1s[@]} lanes for this sample."

############################
### CHECK IF OUTPUT EXISTS
############################
OUT_BAM="$OUT_DIR/${SAMPLE_NAME}_Aligned.sortedByCoord.out.bam"
if [[ -f "$OUT_BAM" ]]; then
    echo "Output BAM exists — skipping STAR for $SAMPLE_NAME."
    exit 0
fi

############################
### CONCATENATE ALL LANES
############################
TMP_DIR="$TEMP_BASE/${SAMPLE_NAME}"
mkdir -p "$TMP_DIR"

TMP_R1="$TMP_DIR/${SAMPLE_NAME}_R1_merged.fastq"
TMP_R2="$TMP_DIR/${SAMPLE_NAME}_R2_merged.fastq"

echo "Concatenating R1 files → $TMP_R1"
cat "${R1s[@]}" > "$TMP_R1"

echo "Concatenating R2 files → $TMP_R2"
cat "${R2s[@]}" > "$TMP_R2"

FASTQ_CMD="--readFilesIn $TMP_R1 $TMP_R2"

############################
### RUN STAR
############################
echo "Running STAR for $SAMPLE_NAME..."

STAR \
    --runThreadN 8 \
    --genomeDir "$GENOME_DIR" \
    $FASTQ_CMD \
    --outFileNamePrefix "$OUT_DIR/${SAMPLE_NAME}_" \
    --outSAMtype BAM SortedByCoordinate \
    --limitBAMsortRAM 55000000000 \
    --outFilterMultimapNmax 1

STAR_EXIT=$?
if (( STAR_EXIT != 0 )); then
    echo "ERROR: STAR failed with exit code $STAR_EXIT"
    exit $STAR_EXIT
fi

############################
### BAM INDEX
############################
echo "Indexing BAM..."
samtools index "$OUT_BAM"

############################
### CLEAN UP TEMP FILES
############################
echo "Cleaning temp files..."
rm -rf "$TMP_DIR"

echo "✔ Finished sample: $SAMPLE_NAME"
exit 0
```
