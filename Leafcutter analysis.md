# Leafcutter and analysis
Need:
Samtools, regtools, python3, R, leafcutter repo
All should be under rnaseq_tools environment

1. Make per-junction sample files from BAM (shell)
```
#!/bin/bash
#SBATCH --job-name=leafcutter_juncs
#SBATCH --output=leafcutter_juncs_%j.out
#SBATCH --error=leafcutter_juncs_%j.err
#SBATCH --mail-type-=ALL
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

# === 1. Load dependencies ===
module load miniconda
conda activate rnaseq_tools

# === 2. Define paths ===
BAM_DIR="/gpfs/gibbs/pi/guo/mg2684/GSE201407/star_output"          # folder with your .bam files
OUT_DIR="/gpfs/gibbs/pi/guo/mg2684/GSE201407/star_output/leafcutter_jxns"    # where junction files will go
SCRIPT_DIR="/gpfs/gibbs/pi/guo/mg2684/tools/leafcutter/scripts"   # location of LeafCutter scripts

# === 3. Run junction extraction ===
mkdir -p $OUT_DIR

for bam in ${BAM_DIR}/*.sortedByCoord.out.bam; do
    sample=$(basename $bam .sortedByCoord.out.bam)
    echo "Processing ${sample}..."
    samtools index $bam  # makes .bai if missing

    # run the leafcutter script to extract junctions
    python ${SCRIPT_DIR}/leafcutter/scripts/leafcutter/clustering/extract_juncs.py \
        -a ${sample} \
        -o ${OUT_DIR}/${sample}.junc \
        $bam
done
```
Notes:
	•	If your BAMs are already coordinate-sorted and indexed, you can skip the sort step. The script checks and sorts into *.sorted.bam.
	•	regtools options used above are examples; -a (anchor length), -m/-M (min/max intron length). You may tune these based on library prep.

2. Combine per-sample junctions (one combined file)
```
cd "$OUT_DIR"
# Add sample name as 6th column and combine (safer)
for f in *.junc; do
  sample=$(basename "$f" .junc)
  awk -v s="$sample" 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,s}' "$f"
done > combined_juncs_with_sample.txt
# optionally gzip
gzip -c combined_juncs_with_sample.txt > combined_juncs_with_sample.txt.gz
```
3. Create phenotype (groups.txt) file from metadata (R)
```
library(readr)
meta <- read_csv("/gpfs/gibbs/pi/guo/mg2684/GSE201407/GSE201407_metadata_clean.csv", show_col_types = FALSE)
meta %>%
  filter(SRR %in% c("SRR18907541","SRR18907548","SRR18907557", "SRR18907558", "SRR18907555", "SRR18907556", "SRR18907552", "SRR18907540", "SRR18907547", "SRR18907544", "SRR18907543", "SRR18907550", "SRR18907549", "SRR18907551", "SRR18907561", "SRR18907562", "SRR18907563", "SRR18907554", "SRR18907545", "SRR18907560", "SRR18907542", "SRR18907546", "SRR18907559", "SRR18907553")) %>%
  select(SRR, disease_state) %>%
  write_tsv("/gpfs/gibbs/pi/guo/mg2684/GSE201407/leafcutter_groups.txt", col_names = FALSE)
meta %>%
  filter(SRR %in% samples) %>%
  select(SRR, disease_state, sex, developmental_stage) %>%
  write_csv("/gpfs/gibbs/pi/guo/mg2684/GSE201407/leafcutter_covariates.csv")
```
4. Run leafcutter clustering (HPC
```
LEAFCUTTER_DIR="/home/youruser/leafcutter"   # update
cd "$LEAFCUTTER_DIR"

# Use combined file (gzip accepted)
python3 $LEAFCUTTER_DIR/leafcutter_cluster.py \
  -j /gpfs/gibbs/pi/guo/mg2684/GSE201407/leafcutter_input_bam/combined_juncs_with_sample.txt.gz \
  -o /gpfs/gibbs/pi/guo/mg2684/GSE201407/leafcutter_clusters \
  --minclustsize 30    # optional: minimum reads per cluster; tune as needed
```
5. Run leafcutter differential splicing
```
Rscript $LEAFCUTTER_DIR/leafcutter_ds.R \
  -i /gpfs/gibbs/pi/guo/mg2684/GSE201407/leafcutter_clusters/perind_numers.counts.gz \
  -g /gpfs/gibbs/pi/guo/mg2684/GSE201407/leafcutter_groups.txt \
  -o /gpfs/gibbs/pi/guo/mg2684/GSE201407/leafcutter_ds_out \
  --min_samples_per_intron 3 \
  --min_reads_per_intron 10
```
6. Visualize leafcutter results
	•	LeafCutter will output cluster-level significance files you can open in R.
	•	For sashimi-style validation you need BAMs — now that you used BAMs, you can plot reads across the junction region in IGV or use ggsashimi/leafcutter plotting helpers. I can produce a sashimi job script for the top clusters if you want.
