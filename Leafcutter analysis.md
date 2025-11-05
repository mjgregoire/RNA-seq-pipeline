# Leafcutter and analysis
Need:
Samtools, regtools, python3, R, leafcutter repo
All should be under rnaseq_tools environment

1. Make per-junction sample files from BAM (shell)
```
#first add this line to the bam2junc.sh file in the leafcutter/scripts directory
leafCutterDir="/gpfs/gibbs/pi/guo/mg2684/tools/leafcutter"
```
```
#!/bin/bash
#SBATCH --job-name=leafcutter_juncs
#SBATCH --output=leafcutter_juncs_%j.out
#SBATCH --error=leafcutter_juncs_%j.err
#SBATCH --mail-type=ALL
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=day
set -euo pipefail
echo "=== Starting LeafCutter junction extraction job ==="
date
echo "Running on node: $(hostname)"
echo "Working directory: $(pwd)"
# === 1. Clean environment ===
module --force purge
module load miniconda
# === 2. Initialize Conda manually (robust method) ===
# Use the conda init script provided by the module itself
if command -v conda &> /dev/null; then
    echo "✅ Conda already available"
else
    echo "🔧 Initializing Conda manually..."
    source $(module show miniconda 2>&1 | grep 'CONDA_DIR' | awk '{print $3}')/etc/profile.d/conda.sh 2>/dev/null || \
    source $(conda info --base)/etc/profile.d/conda.sh 2>/dev/null || true
fi
# If the above failed, try fallback to your home installation
if ! command -v conda &> /dev/null; then
    echo "⚠️ Conda still not found — trying ~/.bashrc"
    source ~/.bashrc || true
fi
# === 3. Activate environment ===
conda activate rnaseq_tools || { echo "❌ Failed to activate rnaseq_tools"; exit 1; }
# === 4. Confirm environment ===
echo "Conda environment: $(conda info --envs | grep '*' | awk '{print $1}')"
which samtools
samtools --version | head -n 1
which python
which R
# === 5. Define paths ===
BAM_DIR="/gpfs/gibbs/pi/guo/mg2684/GSE201407/star_output"
OUT_DIR="/gpfs/gibbs/pi/guo/mg2684/GSE201407/leafcutter_jxns"
SCRIPT_DIR="/gpfs/gibbs/pi/guo/mg2684/tools/leafcutter/scripts"
# === 6. Run junction extraction ===
mkdir -p "$OUT_DIR"
echo "Starting junction extraction for BAMs in: $BAM_DIR"
echo "Output directory: $OUT_DIR"
echo "Script directory: $SCRIPT_DIR"
echo "----------------------------------------------------"
for bam in ${BAM_DIR}/*.sortedByCoord.out.bam; do
    sample=$(basename "$bam" .sortedByCoord.out.bam)
    if [ ! -s "$bam" ]; then
        echo "⚠️ Skipping empty BAM: $bam"
        continue
    fi
    echo "Processing ${sample}..."
    samtools index "$bam"
    bash "${SCRIPT_DIR}/bam2junc.sh" "$bam" "${OUT_DIR}/${sample}.junc"
    echo "✅ Done: ${sample}"
    echo "----------------------------------------------------"
done
echo "=== All junctions processed successfully ==="
date
```

2. Combine per-sample junctions (one combined file)
```
cd "$OUT_DIR"
# Add sample name as 6th column and combine (safer)
for f in *.junc; do
  sample=$(basename "$f" .junc)
  awk -v s=$sample '{print $1"\t"$2"\t"$3"\t.\t"$5"\t"$6"\t"s}' "$f"
done > combined_juncs_with_sample.txt
#add header
sed -i '1ichrom\tstart\tend\tstrand\tcount\tsample_id' combined_juncs_with_sample.txt
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
4. Run leafcutter clustering (HPC) NOTE leafcutter scripts directly from the GitHub repo are python2, if your HPC is python3, they won't work -- make a new environment with python2:
```
conda create -n leafcutter_py2 python=2.7
conda activate leafcutter_py2
conda install numpy scipy pandas
```
```
#!/bin/bash
#SBATCH --job-name=leafcutter_cluster
#SBATCH --output=leafcutter_cluster_%j.out
#SBATCH --error=leafcutter_cluster_%j.err
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=day
echo "=== Starting LeafCutter clustering job ==="
date
echo "Running on node: $(hostname)"
echo "Working directory: $(pwd)"
# === 1. Clean environment ===
module --force purge
module load miniconda
# === 2. Initialize Conda manually (robust method) ===
if command -v conda &> /dev/null; then
    echo "✅ Conda already available"
else
    echo "🔧 Initializing Conda manually..."
    source $(module show miniconda 2>&1 | grep 'CONDA_DIR' | awk '{print $3}')/etc/profile.d/conda.
sh 2>/dev/null || \
    source $(conda info --base)/etc/profile.d/conda.sh 2>/dev/null || true
fi
# Fallback if above failed
if ! command -v conda &> /dev/null; then
    echo "⚠️ Conda still not found — trying ~/.bashrc"
    source ~/.bashrc || true
fi

# === 3. Activate environment ===
conda activate leafcutter_py2 || { echo "❌ Failed to activate leafcutter_py2"; exit 1; }
# === 4. Confirm environment ===
echo "Conda environment: $(conda info --envs | grep '*' | awk '{print $1}')"
which python
python --version
# === 5. Define files and directories ===
LEAFCUTTER_DIR=/gpfs/gibbs/pi/guo/mg2684/tools/leafcutter
OUT_DIR=/gpfs/gibbs/pi/guo/mg2684/GSE201407/leafcutter_jxns/leafcutter_clusters
mkdir -p $OUT_DIR
JUNC_DIR=/gpfs/gibbs/pi/guo/mg2684/GSE201407/leafcutter_jxns
JUNC_FILE=$JUNC_DIR/juncfiles.txt
ls $JUNC_DIR/*.junc > $JUNC_FILE
if [ ! -s "$JUNC_FILE" ]; then
  echo "❌ No junction files found in $JUNC_DIR"
  exit 1
fi
echo "✅ Found $(wc -l < $JUNC_FILE) junction files to process"
# === 6. Run LeafCutter clustering ===
cd $OUT_DIR
echo "🚀 Running LeafCutter clustering..."
# Force Python 2 explicitly
PY2=~/.conda/envs/leafcutter_py2/bin/python
echo "Using Python interpreter: $PY2"

$PY2 $LEAFCUTTER_DIR/clustering/leafcutter_cluster.py \
  -j $JUNC_FILE \
  -m 50 \
  -o leafcutter_clusters
echo "=== LeafCutter clustering complete ==="
date
```
5. Run leafcutter differential splicing
```
#make new leafcutter_R environment
conda create -n leafcutter_R r-base=4.3.1 -y
conda activate leafcutter_R
conda install -c conda-forge r-devtools r-optparse r-tidyverse r-dplyr r-ggplot2 -y
R
	devtools::install_github("davidaknowles/leafcutter/leafcutter")
q()
```
```
#!/bin/bash
#SBATCH --job-name=leafcutter_ds
#SBATCH --output=leafcutter_ds_%j.out
#SBATCH --error=leafcutter_ds_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=day
#SBATCH --mail-type=ALL

echo "=== Starting LeafCutter differential splicing analysis ==="
date
echo "Running on node: $(hostname)"
echo "Working directory: $(pwd)"

# === 1. Clean environment ===
module --force purge
module load miniconda

# === 2. Initialize Conda ===
if ! command -v conda &> /dev/null; then
    echo "🔧 Initializing Conda manually..."
    source $(conda info --base)/etc/profile.d/conda.sh 2>/dev/null || \
    source ~/.bashrc || true
fi

# === 3. Activate environment ===
echo "Activating conda environment: leafcutter_R"
conda activate leafcutter_R || { echo "❌ Failed to activate leafcutter_R"; exit 1; }

# === 4. Confirm environment ===
echo "✅ Environment activated: $(conda info --envs | grep '*' | awk '{print $1}')"
which Rscript
Rscript --version

# === 5. Define directories ===
LEAFCUTTER_DIR=/gpfs/gibbs/pi/guo/mg2684/tools/leafcutter
CLUSTER_DIR=/gpfs/gibbs/pi/guo/mg2684/GSE201407/leafcutter_clusters
GROUP_FILE=/gpfs/gibbs/pi/guo/mg2684/GSE201407/leafcutter_jxns/leafcutter_groups.txt
COV_FILE=/gpfs/gibbs/pi/guo/mg2684/GSE201407/leafcutter_jxns/leafcutter_covariates.csv
OUT_DIR=/gpfs/gibbs/pi/guo/mg2684/GSE201407/leafcutter_ds_out
OUT_DIR_COV=/gpfs/gibbs/pi/guo/mg2684/GSE201407/leafcutter_ds_out_cov

mkdir -p $OUT_DIR $OUT_DIR_COV

# === 6. Run standard LeafCutter DS analysis ===
echo "🚀 Running standard DS analysis (no covariates)..."
Rscript $LEAFCUTTER_DIR/scripts/leafcutter_ds.R \
  -i $CLUSTER_DIR/leafcutter_clusters_perind_numbers.counts.gz \
  -g $GROUP_FILE \
  -o $OUT_DIR \
  -p 8 \
  --min_samples_per_intron 3 \
  --min_reads_per_intron 10
echo "✅ Standard DS analysis complete"
date

# === 7. Run covariate-adjusted LeafCutter DS analysis ===
echo "🚀 Running covariate-adjusted DS analysis..."
Rscript $LEAFCUTTER_DIR/scripts/leafcutter_ds.R \
  -i $CLUSTER_DIR/leafcutter_clusters_perind_numbers.counts.gz \
  -g $GROUP_FILE \
  -o $OUT_DIR_COV \
  -p 8 \
  --min_samples_per_intron 3 \
  --min_reads_per_intron 10 \
  --covariates $COV_FILE
echo "✅ Covariate-adjusted DS analysis complete"
date

echo "=== All LeafCutter DS analyses finished successfully ==="
```
6. Visualize leafcutter results
	•	LeafCutter will output cluster-level significance files you can open in R.
	•	For sashimi-style validation you need BAMs — now that you used BAMs, you can plot reads across the junction region in IGV or use ggsashimi/leafcutter plotting helpers. I can produce a sashimi job script for the top clusters if you want.
