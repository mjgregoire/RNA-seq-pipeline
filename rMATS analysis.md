# rMATS analysis
in R set up:
```
#LOAD LIBRARIES
suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(pheatmap)
  library(tidyr)
  library(stringr)
  library(GEOquery)
  library(rentrez)
  library(tibble)
  library(DESeq2)
  library(future)
  library(furrr)
  library(ggrepel)
  library(GenomicRanges)
  library(rtracklayer)
  library(gt)
  library(DESeq2)
  library(BiocParallel)
})

# --- MAKE SAMPLESHEET LIST ---
#load metadata to make master file with info on BAMS
meta_clean <- read_csv("/gpfs/gibbs/pi/guo/mg2684/GSE201407/GSE201407_metadata_clean.csv")
colnames(meta_clean)
# path to your bam folder
bam_dir <- "/gpfs/gibbs/pi/guo/mg2684/GSE201407/star_output/bam_files"
# list all BAMs (not .bai)
bam_files <- list.files(bam_dir, pattern = "Aligned.sortedByCoord.out.bam$", full.names = TRUE)
# make a tibble with SRR extracted
bam_df <- tibble(
  bam_path = bam_files,
  SRR = str_extract(basename(bam_files), "SRR[0-9]+")
)
head(bam_df)
#merge with metadata
meta_joined <- meta_clean %>%
  inner_join(bam_df, by = "SRR")
dim(meta_clean)
dim(meta_joined)
head(meta_joined)

samplesheet <- meta_joined %>%
  transmute(
    sample_id = sample_id,
    SRR = SRR,
    bam_path = bam_path,
    disease_state = disease_state,   # or genotype
    timepoint = developmental_stage, # e.g., "day21", "day42"
    sex = sex
  )
samplesheet <- samplesheet %>%
  mutate(
    timepoint = recode(timepoint, "21 days in vitro" = "D21", "42 days in vitro" = "D42"),
    disease_state = recode(disease_state, "Healthy control" = "CTRL", "ALS" = "ALS")
  )

head(samplesheet)

write_csv(samplesheet, "/gpfs/gibbs/pi/guo/mg2684/GSE201407/star_output/samplesheet_master.csv")
read_csv("/gpfs/gibbs/pi/guo/mg2684/GSE201407/star_output/samplesheet_master.csv")


# --- PREP .TXT FILES FOR RMATS ---
# ALS vs CTRL at Day21
samplesheet %>%
  filter(timepoint == "D21", disease_state == "ALS") %>%
  pull(bam_path) %>%
  write_lines("/gpfs/gibbs/pi/guo/mg2684/GSE201407/star_output/bams_ALS_D21.txt")

samplesheet %>%
  filter(timepoint == "D21", disease_state == "CTRL") %>%
  pull(bam_path) %>%
  write_lines("/gpfs/gibbs/pi/guo/mg2684/GSE201407/star_output/bams_CTRL_D21.txt")

# ALS vs CTRL at Day42
samplesheet %>%
  filter(timepoint == "D42", disease_state == "ALS") %>%
  pull(bam_path) %>%
  write_lines("/gpfs/gibbs/pi/guo/mg2684/GSE201407/star_output/bams_ALS_D42.txt")

samplesheet %>%
  filter(timepoint == "D42", disease_state == "CTRL") %>%
  pull(bam_path) %>%
  write_lines("/gpfs/gibbs/pi/guo/mg2684/GSE201407/star_output/bams_CTRL_D42.txt")
```
Now you can make sbatch script to run rmats! --change files for each contrast you want
```
#make the files ascii and unix compatible with rmats
# Clean control list
dos2unix /gpfs/gibbs/pi/guo/mg2684/GSE201407/star_output/bams_CTRL_D21.txt
sed -i '/^$/d' /gpfs/gibbs/pi/guo/mg2684/GSE201407/star_output/bams_CTRL_D21.txt

# Clean ALS list
dos2unix /gpfs/gibbs/pi/guo/mg2684/GSE201407/star_output/bams_ALS_D21.txt
sed -i '/^$/d' /gpfs/gibbs/pi/guo/mg2684/GSE201407/star_output/bams_ALS_D21.txt

paste -sd, bams_CTRL_D21.txt > bams_CTRL_D21_fixed.txt
paste -sd, bams_ALS_D21.txt > bams_ALS_D21_fixed.txt
paste -sd, bams_CTRL_D21.txt > bams_CTRL_D21_fixed.txt
paste -sd, bams_ALS_D42.txt > bams_ALS_D42_fixed.txt
```
neither of those worked with rmats so move bam files into folder and then try making file name in command line: 
```
#use code like this and repeat for all samples
mkdir -p /gpfs/gibbs/pi/guo/mg2684/GSE201407/star_output/bam_files/ALS_D42

for srr in SRR18907546 SRR18907547 SRR18907548 SRR18907549 SRR18907550 SRR18907551; do
  cp /gpfs/gibbs/pi/guo/mg2684/GSE201407/star_output/bam_files/${srr}_Aligned.sortedByCoord.out.bam \
     /gpfs/gibbs/pi/guo/mg2684/GSE201407/star_output/bam_files/ALS_D42/
done
```
```
#!/bin/bash
#SBATCH --job-name=rmats_ALS_vs_CTRL_D21
#SBATCH --output=rmats_%j.out
#SBATCH --error=rmats_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=90G
#SBATCH --mail-type=END,FAIL

# -----------------------
# Load environment
# -----------------------
module purge
module load rMATS-turbo/4.2.0-foss-2022b

# -----------------------
# Define paths
# -----------------------
CTRL_LIST=/gpfs/gibbs/pi/guo/mg2684/GSE201407/star_output/bam_files/CTRL_D21_bams.txt
ALS_LIST=/gpfs/gibbs/pi/guo/mg2684/GSE201407/star_output/bam_files/ALS_D21_bams.txt
GTF=/gpfs/gibbs/pi/guo/mg2684/reference/gencode/gencode.v43.annotation.gtf
OUTDIR=/gpfs/gibbs/pi/guo/mg2684/GSE201407/rMATS_analysis/rmats_ALS_vs_CTRL_D21

# -----------------------
# Make comma-separated BAM lists (rMATS prefers this format)
# -----------------------
B1=$(paste -sd, $CTRL_LIST)
B2=$(paste -sd, $ALS_LIST)

echo "Control BAMs: $B1"
echo "ALS BAMs: $B2"

mkdir -p $OUTDIR/tmp

# -----------------------
# Run rMATS
# -----------------------
python $(which rmats.py) \
  --b1 "$B1" \
  --b2 "$B2" \
  --gtf $GTF \
  --od $OUTDIR \
  --tmp $OUTDIR/tmp \
  --readLength 150 \
  --nthread 12 \
  --libType fr-firststrand \
  --variable-read-length \
  --novelSS
```
for f in $(cat /gpfs/gibbs/pi/guo/mg2684/GSE201407/star_output/bams_CTRL_D21_fixed.txt | tr ',' ' '); do
  if [[ ! -f "$f" ]]; then echo "❌ Missing: $f"; else echo "✅ Found: $f"; fi
done
