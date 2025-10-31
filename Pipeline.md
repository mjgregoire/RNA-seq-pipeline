# RNA-seq Pipeline Master Script 
**Purpose: Reproducible workflow for RNA-seq analysis. This project is currently in progress and will be updated regularly from Agust 2025 through and to November 2025**

## Log in to HPC and load conda
Log in to HPC with [OpenOnDemand GUI](https://docs.ycrc.yale.edu/clusters-at-yale/access/ood/) or via SSH on your computer's terminal. 

To log in with SSH you need to create a specific key do this with the following:
```
ssh-keygen
cat ~/YaleSSHkey.pub #I named my key "YaleSSHkey" but you can name it whatever
```
Upload the key to Yale: https://sshkeys.ycrc.yale.edu/ 

Check that you can SSH into the Yale HPC and have the proper settings: 
```
ssh -i ~/YaleSSHkey {your ID}@mccleary.ycrc.yale.edu
nano ~/.ssh/config: Host mccleary
    HostName mccleary.ycrc.yale.edu
    User {your ID}
    IdentityFile ~/YaleSSHkey
chmod 600 ~/.ssh/config
chmod 600 ~/YaleSSHkey
chmod 700 ~/.ssh
```
Now you can use: `ssh -i ~/YaleSSHkey {your ID}@mccleary.ycrc.yale.edu`


To work in the HPC you need to activate an environment with packages will we use for RNA seq qnalysis.
Conda environment: rnaseq_tools (Packages: sra-tools, entrez-direct, fastqc, multiqc, samtools, etc...)

`module load miniconda`
`conda activate rnaseq_tools`

This was installed via: 

`conda install -c bioconda {package}` `conda create -n{name} -c bioconda {packages separated by spaces}` 
ex: `conda install -c conda-forge r-base`

You can check what packages are you the environment with `conda list`.

With the environment loaded you can install things needed in it later with the `conda install -c bioconda {package}` code

## Set up/go to working directory
`mkdir {path to your directory here, make an RNA seq folder and then sub folders for each project}`
`cd {path to your directory}`

### For downloading FASTQ from SRA projects or ENA projects with SRA tools, make accession list first
Get the SRR accession list from SRA Project (SRP). But you *don't want to use the login node* on HPC b/c will take up too much compute power and others can't do stuff, so request some time in a cluster to do manual work: `salloc --mem=1G --time=1:00:00`

Use the following code to get SRR accession list: `esearch -db sra -query {SRA name here} | efetch -format runinfo | cut -d ',' -f 1 | grep SRR > Run_Acc_List.txt`

Or for European files: `wget -qO - {link to ENA e.g.: "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=ERP131847&result=read_run&fields=run_accession"} | tail -n + 2 > Run_Acc_List.txt`

### Run FASTQ download script
This bash script will convert .sra to .fastq files if applicable, and downloads and compresses fastq files.
Create the script in your folder with: `nano download_fastq.slurm`

```
#!/bin/bash
#SBATCH -p day
#SBATCH --job-name=fastq_array
#SBATCH --mem-per-cpu=1G
#SBATCH -t 6:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=logs/fastq_%A_%a.out
#SBATCH --error=logs/fastq_%A_%a.err

module load miniconda
conda activate rnaseq_tools

# --- SETUP ---
WORKDIR="/gpfs/gibbs/pi/guo/mg2684/ERP126666"   # <-- change to your path
ACC_LIST="${WORKDIR}/SRR_Acc_List.txt"          # or Run_Acc_List.txt if you used that name

cd "$WORKDIR"

# --- DETERMINE ARRAY SIZE ---
NUM_JOBS=$(wc -l < "$ACC_LIST")
ARRAY_MAX=$((NUM_JOBS - 1))

# If SLURM_ARRAY_TASK_ID is not set (e.g., running interactively), exit
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "This script is meant to be run with sbatch using --array=0-${ARRAY_MAX}"
    exit 1
fi

# --- GET SRR ID ---
SRR=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$ACC_LIST")

echo "[$(date)] Task $SLURM_ARRAY_TASK_ID: Downloading $SRR"

# --- DOWNLOAD ---
fasterq-dump --split-files --threads 4 --outdir fastq "$SRR"

# --- COMPRESS OUTPUT ---
gzip "${SRR}_1.fastq" "${SRR}_2.fastq"

echo "[$(date)] Finished $SRR"
```

Run the script with the following commmand: `sbatch --array=0-$(( $(wc -l < Run_Acc_List.txt) - 1 )) download_fastq.slurm`
Using sbatch on bash scripts will send them to an available cluster and they will run in the background until completed, the time runs out, or they run into an error. You can at this point run other scripts or leave the command line and your script should still run. 

**Verify output from download**
Check that the bash script is working with: `squeue --me`
Or with the output files using: `tail -f ~{your file path here}/logs/*.out`

Run ls or check file counts: `ls -lh fastq_files | grep '.fastq.gz'`
`wc -l SRR_Acc_List.txt`
`ls *_1.fastq.gz | wc -l`

`ls *_1.fastq.gz | sed 's/_1.fastq.gz//' | sort > downloaded.txt`
`comm -23 <(sort SRR_Acc_List.txt) downloaded.txt > not_downloaded.txt`

## Quality control on FASTQ files with FASTQC
Run fastqc as as bash script. 

Make a folder for the results: `mkdir fastqc_results`

Make the script using: `nano fastqc.sh`

```
#!/bin/bash
#SBATCH -p day
#SBATCH --job-name=fastqc
#SBATCH --mail-type=ALL
#SBATCH --mem=4G
#SBATCH -t 2:00:00 
#SBATCH --array=0-49 #change this to the number of files you have
#SBATCH --output=logs/fastqc_%A_%a.out
#SBATCH --error=logs/fastqc_%A_%a.err
module load miniconda
conda activate rnaseq_tools
cd ~{path to your folders here}
# Get the nth fastq.gz file
FILE=$(ls *.fastq.gz | sort | sed -n "$((SLURM_ARRAY_TASK_ID+1))p")
echo "Running FastQC on $FILE"
fastqc "$FILE" --outdir ~{path to your folders here}`
```

OR: 
```
#!/bin/bash
#SBATCH -p day
#SBATCH --job-name=fastqc
#SBATCH --mail-type=ALL
#SBATCH --mem=4G
#SBATCH -t 2:00:00
#SBATCH --output=logs/fastqc_%A_%a.out
#SBATCH --error=logs/fastqc_%A_%a.err

module load miniconda
eval "$(conda shell.bash hook)"   # ensures conda works in sbatch
conda activate rnaseq_tools

cd /gpfs/gibbs/pi/guo/mg2684/GSE201407/fastqs 
#OR:
# Use the directory you are in when you type sbatch
#cd "$SLURM_SUBMIT_DIR" || exit 1   # exit if cd fails


# Make a list of all fastq.gz files
FILES=( $(ls *.fastq.gz | sort) )
NFILES=${#FILES[@]}

# Safety check: prevent running out of range
if [ "$SLURM_ARRAY_TASK_ID" -ge "$NFILES" ]; then
    echo "Error: SLURM_ARRAY_TASK_ID ($SLURM_ARRAY_TASK_ID) >= number of FASTQ files ($NFILES)"
    exit 1
fi

# Pick the correct file
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

echo "Running FastQC on $FILE"
fastqc "$FILE" --outdir /gpfs/gibbs/pi/guo/mg2684/GSE201407/fastqs/fastqc_results
```

Run fastqc with: `sbatch fastqc.sh`
Then compile all fastqc files into one using: `multiqc .`

To open the multiqc .html report you need to transfer the file to your local downloads folder. To do this you need to go to your local terminal. If you have not already set up an ssh key (see above loging into HPC section to login with SSH).

On your terminal type the following to SSH into the HPC and download the file locally (this downloads to your Downloads folder): 
`scp -i ~/YaleSSHkey {your ID}@transfer-mccleary.ycrc.yale.edu:{path to your accesion here}/fastqs/fastqc_results/multiqc_report.html \ 

~/Downloads/multiqc_report_$(basename $(dirname $(dirname $(dirname /{path to your accession here}/fastqs/fastqc_results/multiqc_report.html)))).html`
The backslash in the code is a line continuation character and tells terminal to continue onto the next line.

## Trimming FASTQ files
To trim poor sequence quality or adapters on FASTQ files make and run the following bash script. 
`nano trim_fastp.sh`
`sbatch trim_fastp.sh`

```
#!/bin/bash
#SBATCH --job-name=fastp_trim
#SBATCH --output=fastp_trim.%A_%a.out
#SBATCH --error=fastp_trim.%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL 

# -----------------------------
# Load conda safely
# -----------------------------
module purge
module load miniconda
eval "$(conda shell.bash hook)"
conda activate rnaseq_tools

# -----------------------------
# Directories
# -----------------------------
IN_DIR=/gpfs/gibbs/pi/guo/mg2684/GSE201407/fastqs #change this to your directory!
OUT_DIR="${IN_DIR}/trimmed_fastq"
mkdir -p "$OUT_DIR"

# -----------------------------
# Trimming loop
# -----------------------------
for R1 in "$IN_DIR"/*_1.fastq.gz; do
    BASE=$(basename "$R1" _1.fastq.gz)
    R2="${IN_DIR}/${BASE}_2.fastq.gz"

    # Skip if R2 doesn't exist
    if [ ! -f "$R2" ]; then
        echo "Skipping $BASE: paired file $R2 not found"
        continue
    fi

    echo "Processing sample: $BASE"

    # Estimate read length from the first read
    READ_LEN=$(zcat "$R1" | head -n 2 | tail -n 1 | wc -c)
    READ_LEN=$((READ_LEN - 1))

    # Set minimum length depending on read length
    if [ "$READ_LEN" -le 50 ]; then
        MIN_LEN=20
    else
        MIN_LEN=50
    fi

    # Extract accession folder (two levels up from fastq_results)
    ACCESSION=$(basename "$(dirname "$(dirname "$IN_DIR")")")

    # Run fastp
    fastp \
      -i "$R1" \
      -I "$R2" \
      -o "$OUT_DIR/${BASE}_R1.trimmed.fastq.gz" \
      -O "$OUT_DIR/${BASE}_R2.trimmed.fastq.gz" \
      --detect_adapter_for_pe \
      --cut_front \
      --cut_tail \
      --cut_window_size 4 \
      --cut_mean_quality 15 \
      --length_required "$MIN_LEN" \
      --thread 4 \
      --html "$OUT_DIR/${BASE}_${ACCESSION}_fastp.html" \
      --json "$OUT_DIR/${BASE}_${ACCESSION}_fastp.json"

done

echo "All trimming done."
```
Check the script is working while running with:
```
squeue --me
tail -f fastp_trim.out
wc -l fastp_trim.out
```
Now you can re-Run FASTQC using the same script as before and check that the sequences have better quality now:

```
mkdir trimmed_fastqc_results
nano trimmed_fastqc.sh

    #!/bin/bash
    #SBATCH -p day
    #SBATCH --job-name=fastqc
    #SBATCH --mail-type=ALL
    #SBATCH --mem=4G
    #SBATCH -t 2:00:00
    #SBATCH --array=0-103
    #SBATCH --output=logs/fastqc_%A_%a.out
    #SBATCH --error=logs/fastqc_%A_%a.err

    module load miniconda
    eval "$(conda shell.bash hook)"
    conda activate rnaseq_tools

    cd /gpfs/gibbs/pi/guo/mg2684/ERP131847/fastqs/trimmed_fastq

    # Get the nth fastq.gz file
    FILE=$(ls *.fastq.gz | sort | sed -n "$((SLURM_ARRAY_TASK_ID+1))p")
    if [ -z "$FILE" ]; then
        echo "No file found for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
        exit 1
    fi

    echo "Running FastQC on $FILE"
    fastqc "$FILE" --outdir /gpfs/gibbs/pi/guo/mg2684/ERP131847/fastqs/trimmed_fastq/trimmed_fastqc_results

sbatch trimmed_fastqc.sh

#compile all fastqc files into one
multiqc .
#transfer fastqc file to local downloads folder to check the html
#go to local terminal
#scp -i ~/YaleSSHkey \
mg2684@transfer-mccleary.ycrc.yale.edu:/gpfs/gibbs/pi/guo/mg2684/ERP131847/fastqs/trimmed_fastq/trimmed_fastqc_results/multiqc_report.html \
~/Downloads/multiqc_report_trimmed_$(basename $(dirname $(dirname $(dirname $(dirname /gpfs/gibbs/pi/guo/mg2684/ERP131847/fastqs/trimmed_fastq/trimmed_fastqc_results/multiqc_report.html))))).html
```
#----ALIGNING FASTQ TO REFERENCE GENOME USING STAR----
# Make a directory to store your references
mkdir -p /gpfs/gibbs/pi/guo/mg2684/reference/gencode
cd /gpfs/gibbs/pi/guo/mg2684/reference/gencode

# Genome FASTA (GRCh38 primary assembly)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz

# Annotation GTF
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz
gunzip gencode.v43.annotation.gtf.gz

In bash script directory:
```
#!/bin/bash
#SBATCH --job-name=STAR_index
#SBATCH --output=STAR_index_%j.out
#SBATCH --error=STAR_index_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=60G

module load STAR/2.7.11a-GCC-12.2.0

GENOME_FASTA=/gpfs/gibbs/pi/guo/mg2684/reference/gencode/GRCh38.primary_assembly.genome.fa
ANNOTATION_GTF=/gpfs/gibbs/pi/guo/mg2684/reference/gencode/gencode.v43.annotation.gtf
STAR_INDEX_DIR=/gpfs/gibbs/pi/guo/mg2684/reference/STAR_index_GRCh38

mkdir -p $STAR_INDEX_DIR

STAR --runThreadN 12 \
     --runMode genomeGenerate \
     --genomeDir $STAR_INDEX_DIR \
     --genomeFastaFiles $GENOME_FASTA \
     --sjdbGTFfile $ANNOTATION_GTF \
     --sjdbOverhang 148
```
`sbatch star_index.sh`

```
nano STAR_align_array.sh

#SBATCH -o STAR_align_%A_%a.out
#SBATCH -e STAR_align_%A_%a.err
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mail-type=ALL
#SBATCH --mem=60G
#SBATCH -t 24:00:00
#SBATCH --array=0-23   # <-- adjust this number to match how many samples you have
# ======== Paths ========
FASTQ_DIR=/gpfs/gibbs/pi/guo/mg2684/GSE201407/fastqs/trimmed_fastq
OUT_DIR=/gpfs/gibbs/pi/guo/mg2684/GSE201407/star_output
GENOME_DIR=/gpfs/gibbs/pi/guo/mg2684/reference/STAR_index_GRCh38
# ======== Load environment ========
module load STAR/2.7.11a
# or, if using your conda env:
# source activate rnaseq_tools
# ======== List all SRR samples ========
SAMPLES=($(ls ${FASTQ_DIR}/*_1.fastq.gz | xargs -n 1 basename | sed 's/_1.fastq.gz//'))
# ======== Select current sample based on SLURM_ARRAY_TASK_ID ========
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
# ======== Create output directory if not exists ========
mkdir -p ${OUT_DIR}
# ======== Run STAR alignment ========
STAR \
--runThreadN 12 \
--genomeDir ${GENOME_DIR} \
--readFilesIn ${FASTQ_DIR}/${SAMPLE}_1.fastq.gz ${FASTQ_DIR}/${SAMPLE}_2.fastq.gz \
--readFilesCommand zcat \
--outFileNamePrefix ${OUT_DIR}/${SAMPLE}_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts
```

# --- LOOKING AT SPLICING VIA JUNCTION COUNTS ---
## 1) Extract STAR junction counts (fast, per-sample)
`cd /path/to/star_output   # change this to your actual path`

### make a folder to store the extracted junction count files
`mkdir -p junction_counts`

### loop over all SJ.out.tab files
```
for f in *_SJ.out.tab; do
sample_name=$(basename "$f" _SJ.out.tab)
awk 'BEGIN{OFS="\t"}{print $1":"$2"-"$3":"$4, $6}' "$f" \
> "junction_counts/${sample_name}.junction_counts.tsv"
sed -i '1iJUNC\tCOUNT' "junction_counts/${sample_name}.junction_counts.tsv"
echo "Processed ${sample_name}"
done
```

## 2) Collect total mapped reads per sample (needed for normalization later)
`cd /path/to/star_output`

### create output file for total mapped reads
`echo -e "SAMPLE\tMAPPED_READS" > sample_mapped_reads.tsv`

### extract "Uniquely mapped reads number" from each Log.final.out
```
for f in *_Log.final.out; do
sample_name=$(basename "$f" _Log.final.out)
mapped=$(grep "Uniquely mapped reads number" "$f" | awk '{print $6}')
echo -e "${sample_name}\t${mapped}" >> sample_mapped_reads.tsv
echo "Extracted mapped reads for ${sample_name}"
done
```

## 3) Run featureCounts to get gene-level counts
Make sure you have:
	•	All your sorted BAMs in /path/to/star_output/ (they end with _Aligned.sortedByCoord.out.bam)
	•	Your GTF annotation file (for example: /path/to/annotation/gencode.vM33.annotation.gtf)
	•	featureCounts installed (part of the Subread package)

via sbatch:

```
nano featureCounts

#!/bin/bash
#SBATCH --job-name=featureCounts
#SBATCH --output=featureCounts_%j.out
#SBATCH --error=featureCounts_%j.err
#SBATCH --time=4:00:00
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

module load subread  # if your cluster uses module system
featureCounts -T 8 -p -t exon -g gene_id \
  -a /gpfs/gibbs/pi/guo/mg2684/reference/gencode/gencode.v43.annotation.gtf \
  -o /gpfs/gibbs/pi/guo/mg2684/GSE201407/star_output/gene_counts.txt \
  /gpfs/gibbs/pi/guo/mg2684/GSE201407/star_output/SRR*_Aligned.sortedByCoord.out.bam

sbatch featureCounts
```
```
column -t gene_counts.txt.summary | less -S
ex:
Status                         /gpfs/gibbs/pi/guo/mg2684/GSE201407/star_output/SRR18907>
Assigned                       50791756                                                >
Unassigned_Unmapped            0                                                       >
Unassigned_Read_Type           0                                                       >
Unassigned_Singleton           0                                                       >
Unassigned_MappingQuality      0                                                       >
Unassigned_Chimera             0                                                       >
Unassigned_FragmentLength      0                                                       >
Unassigned_Duplicate           0                                                       >
Unassigned_MultiMapping        5516900                                                 >
Unassigned_Secondary           0                                                       >
Unassigned_NonSplit            0                                                       >
Unassigned_NoFeatures          3482920                                                 >
Unassigned_Overlapping_Length  0                                                       >
Unassigned_Ambiguity           5729377                                                 
```
Overall there are: ~50.8M Assigned / (50.8 + 5.5 + 3.5 + 5.7 = 65.5M total)
→ ≈ 77–78% of reads assigned
Good range for high-quality RNA-seq.

but with head gene_counts.txt you'll see there are no headers, lets fix that for easy mapping later: 
```
awk 'NR==1 {
  printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6; 
  for (i=7; i<=NF; i++) {
    n=split($i,a,"/"); 
    gsub("_Aligned.sortedByCoord.out.bam","",a[n]); 
    printf "\t"a[n]
  }
  print ""
  next
}1' gene_counts.txt > gene_counts_clean.txt

awk 'NR==2 || NR>2 {print}' gene_counts_clean.txt | \
sed 's/\/.*\///g; s/_Aligned.sortedByCoord.out.bam//g' > gene_counts_final.txt
```

# 4) merege junction counts
```
nano merge_junctions_and_compute_JPM.R

#!/usr/bin/env Rscript
# merge_junctions_and_compute_JPM.R
library(data.table)

# paths - change if needed
junc_dir <- "junction_counts"
mapped_file <- "sample_mapped_reads.tsv"   # from previous step
out_merged <- "merged_junction_counts.tsv"

# read per-sample junction files
files <- list.files(junc_dir, pattern="*.junction_counts.tsv", full.names=TRUE)
if(length(files)==0) stop("No junction files found in junction_counts/")

lst <- lapply(files, function(f){
  dt <- fread(f, header=TRUE)
  samp <- sub("\\.junction_counts\\.tsv$","", basename(f))
  setnames(dt, c("JUNC","COUNT"))
  setnames(dt, "COUNT", samp)
  return(dt)
})

# merge all by JUNC
merged <- Reduce(function(x,y) merge(x,y,by="JUNC", all=TRUE), lst)
merged[is.na(merged)] <- 0
fwrite(merged, out_merged, sep="\t")
cat("Wrote", out_merged, "\n")

# compute JPM and write long table for downstream merging
mapped <- fread(mapped_file) # columns: SAMPLE \t MAPPED_READS
mapped_vec <- setNames(mapped$MAPPED_READS, mapped$SAMPLE)

# make JPM table
jpm <- copy(merged)
sample_cols <- colnames(jpm)[-1]
for(s in sample_cols){
  if(! s %in% names(mapped_vec)) stop(paste("Mapped reads missing for", s))
  jpm[[s]] <- (jpm[[s]] * 1e6) / mapped_vec[[s]]
}
fwrite(jpm, "junctions_JPM.tsv", sep="\t")
cat("Wrote junctions_JPM.tsv\n")

nano run_merge_junctions.sh

#!/bin/bash
#SBATCH --job-name=merge_junctions
#SBATCH --output=merge_junctions_%j.out
#SBATCH --error=merge_junctions_%j.err
#SBATCH --time=02:00:00        # adjust as needed
#SBATCH --mem=32G              # adjust if you have many samples
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL

# Load R if it's a module
module load R/4.3.1   # or whatever version your cluster uses

# (Optional) activate conda environment if you’re using one
# source activate rnaseq_tools

# Move to your working directory
cd /gpfs/gibbs/pi/guo/mg2684/GSE201407/star_output

# Run the R script
Rscript /gpfs/gibbs/pi/guo/mg2684/GSE201407/star_output/merge_junctions_and_compute_JPM.R

```
# 5) Normalize junctions to gene RPKM
Junctions per million (JPM) = Junction read count / Gene RPKM
RPKM = gene counts / gene length in kb x total mapped reads in millions

This scales each junction by how expressed its parent gene is, letting you compare junction usage across samples or between known targets.

## make junction to gene map
make sure that base R environment is loaded contains GenomicFeatures:

# Activate your R environment
module load R-bundle-Bioconductor/3.19-foss-2022b-R-4.4.1

make sure you are in cluster not login node! `salloc`

```
nano make_junction_to_gene_map.R

#!/usr/bin/env Rscript

# === make_junction_to_gene_map.R ===
# Create a junction-to-gene mapping table from a GTF annotation.
# Compatible with GenomicFeatures or txdbmaker, including Bioconductor ≥3.18.

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(dplyr)
  library(readr)
  library(BiocParallel)

  args <- commandArgs(trailingOnly = TRUE)

  # Optional: parse a --threads argument
  threads <- 1
  if (any(grepl("--threads", args))) {
    parsed <- suppressWarnings(
      as.numeric(sub("--threads=", "", args[grepl("--threads", args)]))
    )
    if (!is.na(parsed) && parsed > 0) {
      threads <- parsed
    }
  }

  cat(paste0("🧵 Using ", threads, " threads for parallel processing.\n"))
  register(MulticoreParam(threads))

  suppressWarnings({
    if (requireNamespace("txdbmaker", quietly = TRUE)) {
      library(txdbmaker)
    }
  })
})

# === Input files ===
gtf_file <- "/gpfs/gibbs/pi/guo/mg2684/reference/gencode/gencode.v43.annotation.gtf"
junction_file <- "merged_junction_counts.tsv"

message("📘 Loading GTF annotation: ", gtf_file)

# === 1. Load GTF as a TxDb ===
if ("txdbmaker" %in% .packages()) {
  txdb <- txdbmaker::makeTxDbFromGFF(gtf_file, format = "gtf")
} else {
  txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
}
message("✔ TxDb loaded successfully")

# === 2. Extract introns per transcript ===
introns_by_tx <- intronsByTranscript(txdb, use.names = TRUE)

# Convert to data frame with transcript IDs included, skipping empty ones
introns_list <- bplapply(names(introns_by_tx), function(tx) {
  gr <- introns_by_tx[[tx]]
  if (length(gr) == 0) return(NULL)  # skip single-exon transcripts
  df <- as.data.frame(gr)
  df$tx_id <- tx
  df
})
introns_list <- introns_list[!sapply(introns_list, is.null)]

introns_df <- bind_rows(introns_list)
message("✔ Extracted and flattened introns (", nrow(introns_df), " entries)")

# === 3. Clean up and format ===
introns_df <- introns_df %>%
  mutate(
    intron_start = start,
    intron_end = end,
    chrom = as.character(seqnames),
    strand = as.character(strand)
  ) %>%
  select(chrom, intron_start, intron_end, strand, tx_id)

# === 4. Map transcripts to genes (robust version) ===
tx <- transcripts(txdb, columns = c("tx_id", "gene_id"))
tx2gene <- as.data.frame(mcols(tx)) %>%
  select(tx_id, gene_id) %>%
  distinct() %>%
  mutate(tx_id = as.character(tx_id))  # ensure same type as introns_df

introns_df <- left_join(introns_df, tx2gene, by = "tx_id")

# === 5. Collapse to unique intron–gene pairs ===
junction_to_gene <- introns_df %>%
  distinct(chrom, intron_start, intron_end, strand, gene_id) %>%
  mutate(junction_id = paste(chrom, intron_start, intron_end, strand, sep = "_")) %>%
  select(junction_id, gene_id)

# === 6. Write output ===
write_tsv(junction_to_gene, "junction_to_gene.tsv")

# === 7. Log summary ===
message("✅ Wrote junction_to_gene.tsv with ", nrow(junction_to_gene), " junctions.")
message("🧬 Mapped ", length(unique(junction_to_gene$gene_id)), " unique genes.")


nano run_make_junction_to_gene_map.sh

#!/bin/bash
#SBATCH --job-name=jxn_map
#SBATCH --output=jxn_map_%j.out
#SBATCH --error=jxn_map_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-type=ALL

module --force purge
module load R-bundle-Bioconductor/3.19-foss-2022b-R-4.4.1

Rscript make_junction_to_gene_map.R --threads=8
```

that didn't really work only 1 gene found, so now:

Make a gene-span BED from the GTF
```
cd /gpfs/gibbs/pi/guo/mg2684/GSE201407/star_output

# produce genes.bed (0-based start)
awk '$3=="gene" {
  # remove quotes
  gsub(/"/,"",$0);
  split($0,f,"\t");
  # find gene_id in attributes
  match($0, /gene_id [^;]+/);
  gid=substr($0, RSTART+8, RLENGTH-8);
  # print chrom, start-1, end, gene_id, score(0), strand
  print f[1]"\t"($4-1)"\t"$5"\t"gid"\t0\t"f[7]
}' /gpfs/gibbs/pi/guo/mg2684/reference/gencode/gencode.v43.annotation.gtf \
  > genes.bed

head genes.bed
wc -l genes.bed
```
Make a junctions BED from your merged junction table
```
awk -F'\t' 'NR>1 {
  split($1, p, "_");
  chrom=p[1];
  start=p[2];
  end=p[3];
  strand=p[4];
  # BED format: chrom, start-1, end, junction_id, 0, strand
  print chrom"\t"(start-1)"\t"end"\t"$1"\t0\t"strand
}' merged_junction_counts.tsv > junctions.bed

```

Intersect junctions with gene spans (bedtools)
`bedtools intersect -a junctions.bed -b genes.bed > junctions_annotated.bed`
`bedtools intersect -a junctions.bed -b genes.bed -wa -wb > junctions_intersect_genes.bed`

if you get error: ***** WARNING: File *.bed has inconsistent naming convention for record:
it means that there are non conventional reads.

➡ Step 1: Keep all contigs for the initial discovery phase.
You want to capture every potential novel junction.

➡ Step 2: After detection, filter or flag contigs based on:
	•	Read count or junction support (e.g., ≥3 uniquely mapped reads).
	•	Mapping quality (STAR’s NH or MAPQ tags).
	•	Whether the contig appears in your reference GTF or canonical chromosomes.

➡ Step 3: Focus your biological interpretation on canonical chromosomes, but check non-canonical hits for possible artifacts or interesting exceptions.


So we will flag the canonical chromosomes in the junctions file and re-run bed tools

Add a flag column to junctions.bed:
```
awk 'BEGIN{
    OFS="\t";
}
{
    # canonical chr: chr1-22, chrX, chrY, chrM, or bare 1-22,X,Y,M
    if ($1 ~ /^(chr)?([0-9]{1,2}|X|Y|M)$/)
        flag="canonical";
    else
        flag="noncanonical";
    print $0, flag;
}' junctions.bed > junctions_flagged.bed

awk '$7=="canonical"' junctions_flagged.bed > junctions_canonical.bed
awk '$7=="noncanonical"' junctions_flagged.bed > junctions_noncanonical.bed

bedtools intersect \
    -a junctions_flagged.bed \
    -b genes.bed \
    -wa -wb > junctions_annotated_flagged.bed
```
`awk -F'\t' '{print $4"\t"$10}' junctions_intersect_genes.bed | sort -u > junction_to_gene.tsv`
or using annotated file (DO THIS): 
`awk -F'\t' '{
  junc_id = $4;
  flag = $7;
  gene = $11;         # adjust if your genes.bed stores the gene id/name in a different column
  if (gene == "") gene = ".";
  print junc_id"\t"gene"\t"flag;
}' junctions_annotated_flagged.bed | sort -u > junction_to_gene.tsv`

summarize how many genes each junction maps to:
```
awk '{count[$1]++} END {
  multi=0; single=0;
  for (junc in count) {
    if (count[junc]>1) multi++;
    else single++;
  }
  print "Unique junctions mapping to one gene:", single;
  print "Junctions mapping to multiple genes:", multi;
  print "Total junctions:", single+multi;
}' junction_to_gene.tsv
```
i got: Unique junctions mapping to one gene: 695287
Junctions mapping to multiple genes: 295754
Total junctions: 991041
list mutimappers:
`awk '{print $1}' junction_to_gene.tsv | sort | uniq -d > multi_gene_junctions.txt`
`grep -Ff multi_gene_junctions.txt junction_to_gene.tsv | head`

OR USING THE ANNOTATED FILE:
```
# Count unique vs multi-gene junctions, ignoring canonical flag
cut -f1,2 junction_to_gene.tsv | sort -u | \
awk -F'\t' '{
  count[$1]++
}
END {
  single=multi=0
  for (j in count) {
    if (count[j]>1) multi++; else single++;
  }
  total=single+multi
  print "=== Overall (ignoring flag) ==="
  print "Unique junctions mapping to 1 gene:", single
  print "Junctions mapping to multiple genes:", multi
  print "Total junctions:", total
}'

=== Overall (ignoring flag) ===
Unique junctions mapping to 1 gene: 695287
Junctions mapping to multiple genes: 295754
Total junctions: 991041

# Count unique vs multi-gene junctions separately for canonical/noncanonical
awk -F'\t' '{
  junction=$1; gene=$2; flag=$3;
  key=flag"\t"junction;
  count[key]++
}
END {
  uniq["canonical"]=uniq["noncanonical"]=0
  multi["canonical"]=multi["noncanonical"]=0
  total["canonical"]=total["noncanonical"]=0
  for (k in count) {
    split(k,a,"\t"); flag=a[1];
    total[flag]++;
    if (count[k]>1) multi[flag]++; else uniq[flag]++;
  }
  print "=== Per-flag summary ==="
  for (f in total)
    print f,": total=",total[f],", unique=",uniq[f],", multi-gene=",multi[f]
}' junction_to_gene.tsv

=== Per-flag summary ===
canonical : total= 991041 , unique= 695287 , multi-gene= 295754
noncanonical : total= 0 , unique= 0 , multi-gene= 0


# Create helper file grouping gene lists by junction and flag
awk -F'\t' '{print $1"\t"$3"\t"$2}' junction_to_gene.tsv | sort | \
  awk -F'\t' '{
    j=$1; flag=$2; gene=$3;
    key=j"\t"flag;
    if (seen[key,gene] != 1) { seen[key,gene]=1; genes[key]=genes[key]?genes[key]","gene:gene; cnt[key]++ }
  }
  END {
    PROCINFO["sorted_in"]="@ind_str_asc";
    for (k in genes) {
      split(k,a,"\t");
      junc=a[1]; flag=a[2];
      c=cnt[k];
      if (c==1) uniq[flag]++; else multi[flag]++;
      total[flag]++;
    }
    for (f in total) {
      print f": total junctions="total[f]"; unique="uniq[f]"+0"; multi="multi[f]"+0
    }
  }'


awk -F'\t' '$3=="canonical"{print $1"\t"$2}' junction_to_gene.tsv | sort -u > junction_to_gene.canonical.tsv
awk -F'\t' '$3=="noncanonical"{print $1"\t"$2}' junction_to_gene.tsv | sort -u > junction_to_gene.noncanonical.tsv
```

the annoated file actually didn't keep any of the noncanonical junctions!! we need to fix that!:
```
# create LOJ annotated file
bedtools intersect -a junctions_flagged.bed -b genes.bed -wa -wb -loj > junctions_annotated_flagged_loj.bed

# make junction->gene->flag table (gene = '.' if no overlap)
awk -F'\t' '{
  junc=$4; flag=$7; gene=$11;
  if (gene=="" || gene==".") gene=".";
  print junc"\t"gene"\t"flag;
}' junctions_annotated_flagged_loj.bed | sort -u > junction_to_gene.tsv

# summary: counts overall, canonical vs noncanonical, and non-overlaps
echo "Overall counts (ignoring flag):"
cut -f1,2 junction_to_gene.tsv | sort -u | \
  awk -F'\t' '{cnt[$1]++} END{single=multi=0; for(i in cnt){if(cnt[i]==1) single++; else multi++} print "unique:",single; print "multi:",multi; print "total:",single+multi}'

echo; echo "Per-flag non-overlap counts:"
awk -F'\t' '$2=="." {counts[$3]++} END{for(f in counts) print f,counts[f]}' junction_to_gene.tsv
```
Category
Count
Meaning
unique
743,393
Junctions mapping to exactly one gene
multi
295,754
Junctions overlapping multiple genes
total
1,039,147
Total distinct junctions detected (unique + multi)
canonical (no overlap)
44,948
Canonical junctions not overlapping any annotated gene (likely due to incomplete annotation or alternate contigs)
noncanonical (no overlap)
3,158
Noncanonical junctions not overlapping any annotated gene (true novel/cryptic events)

1.	Canonical no-overlaps (44,948)
These are junctions flagged “canonical” by STAR (GT–AG etc.) but that don’t overlap any gene in your genes.bed.
→ Likely belong to unplaced contigs, alternative haplotypes, or unannotated gene models.
→ They’re not necessarily biologically strange — more likely an annotation gap.
2.	Noncanonical no-overlaps (3,158)
These are the interesting ones for your cryptic splicing analysis!
→ They have noncanonical motifs and don’t overlap any known gene.
→ They could represent true novel junctions — potential cryptic or mis-spliced events.
→ Worth investigating which samples they appear in, whether they’re recurrent, and what the local sequence looks like.
3.	Everything else
The rest (roughly 991k) are junctions that successfully overlap annotated genes — those will cover all normal exon–exon junctions and canonical splicing events.


CREATE SUBSET FOR CRYPTIC ANALYSIS
`awk -F'\t' '$3=="noncanonical" && $2=="."' junction_to_gene.tsv > cryptic_junctions.tsv`
look at genomic coordinates
```
head -n 5 cryptic_junctions.tsv
grep -F $(cut -f1 -d$'\t' cryptic_junctions.tsv | head -n 1) junctions_flagged.bed
```
Later, you can feed those junction coordinates into bedtools getfasta to extract sequences for motif or splice site scoring.

now we can normalize to RPKM!

## normalize
these files need to be in the same folder: 
merged_junction_counts.tsv
gene_counts.txt
junction_to_gene.tsv
```
nano normalize_junctions_to_RPKM.R
# === normalize_junctions_to_RPKM.R ===
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

# === 1. Load data ===
message("📂 Loading input files...")
junction_counts <- read_tsv("merged_junction_counts.tsv")
gene_counts <- read_tsv("gene_counts.txt", comment = "#")
junction_to_gene <- read_tsv("junction_to_gene.tsv", col_names = c("junction_id", "gene_id"))

# === 2. Clean up gene_counts ===
message("🧹 Cleaning gene count column names...")
gene_counts <- gene_counts %>%
  rename(gene_id = Geneid) %>%
  select(gene_id, Length, matches("Aligned.sortedByCoord.out.bam"))

# 🪄 Clean up column names to match SRR IDs in junction file
colnames(gene_counts) <- gsub(".*/|_Aligned.sortedByCoord.out.bam", "", colnames(gene_counts))

# === 3. Compute total mapped reads per sample ===
message("📊 Computing total mapped reads per sample...")

# Identify all SRR-like sample columns
srr_cols <- grep("^SRR\\d+$", colnames(gene_counts), value = TRUE)

if (length(srr_cols) == 0) {
  stop("❌ No SRR sample columns found in gene_counts! Check your column names: ",
       paste(colnames(gene_counts), collapse = ", "))
}

# Ensure the sample columns are numeric
for (c in srr_cols) {
  if (!is.numeric(gene_counts[[c]])) {
    gene_counts[[c]] <- as.numeric(gene_counts[[c]])
  }
}

# Compute total mapped reads per sample
gene_counts_numeric <- gene_counts %>% select(all_of(srr_cols))
total_counts <- colSums(gene_counts_numeric, na.rm = TRUE)

# === 4. Compute RPKM per gene per sample (create *_gene columns) ===
message("⚙️ Computing gene RPKM columns...")
rpkm <- gene_counts %>% select(gene_id, Length)

for (s in names(total_counts)) {
  new_col <- paste0(s, "_gene")
  rpkm[[new_col]] <- (gene_counts[[s]] * 1e9) / (gene_counts$Length * total_counts[[s]])
}

# === 5. Load canonical/noncanonical flag data ===
message("🏷️ Adding canonical/noncanonical flags (if available)...")

if (file.exists("junctions_flagged.bed")) {
  junction_flags <- read_tsv(
    "junctions_flagged.bed",
    col_names = c("chrom", "start", "end", "junction_id", "score", "strand", "flag"),
    col_types = "ciicccc"
  ) %>% select(junction_id, flag)
} else {
  message("⚠️ No 'junctions_flagged.bed' found — proceeding without flags.")
  junction_flags <- tibble(junction_id = junction_counts$junction_id, flag = "unknown")
}

# === 6. Clean and standardize IDs + merge ===
message("🧹 Cleaning and standardizing IDs before merging...")

# Strip version numbers from Ensembl gene IDs
junction_to_gene <- junction_to_gene %>%
  mutate(gene_id = str_replace(gene_id, "\\.\\d+$", ""))

gene_counts <- gene_counts %>%
  mutate(gene_id = str_replace(gene_id, "\\.\\d+$", ""))

# Merge everything together
message("🔗 Merging junctions, gene IDs, RPKMs, and flags...")

junction_norm <- junction_counts %>%
  left_join(junction_to_gene, by = "junction_id") %>%
  left_join(rpkm, by = "gene_id", suffix = c("_junction", "_gene")) %>%
  left_join(junction_flags, by = "junction_id")

# Check merge success
mapped <- sum(!is.na(junction_norm$gene_id))
total <- nrow(junction_norm)
message(sprintf("\n🔍 Junctions mapped to genes: %d / %d (%.2f%%)",
                mapped, total, 100 * mapped / total))

# === 7. Compute normalized junction ratio (Junctions Per Million RPKM) ===
message("⚙️ Computing normalized junction ratios...")

# Detect sample IDs from junction_counts (excluding junction_id)
samples <- setdiff(colnames(junction_counts), "junction_id")

# Figure out which columns in junction_norm correspond to gene RPKMs
rpkm_cols <- grep("_gene$", colnames(junction_norm), value = TRUE)

# ✅ Diagnostic check
message("🧪 Checking column name alignment...")
message("Junction samples: ", paste(samples, collapse = ", "))
message("RPKM columns: ", paste(rpkm_cols, collapse = ", "))

for (s in samples) {
  # Find matching junction and gene columns
  junction_col <- s
  gene_col <- rpkm_cols[grepl(s, rpkm_cols)]

  if (length(gene_col) == 1 && junction_col %in% colnames(junction_norm)) {
    norm_col <- paste0(s, "_JPM")
    junction_norm[[norm_col]] <- junction_norm[[junction_col]] / (junction_norm[[gene_col]] + 1e-6)
  } else {
    message(sprintf("⚠️ Skipping sample %s (missing columns)", s))
  }
}

# === 8. Quick summary ===
message("\n📈 Summary of merged data:")
summary_df <- junction_norm %>%
  count(flag) %>%
  rename(Junctions = n)
print(summary_df)

if ("gene_id" %in% colnames(junction_norm)) {
  mapped <- sum(!is.na(junction_norm$gene_id))
  total <- nrow(junction_norm)
  message(sprintf("\n🔍 Junctions mapped to genes: %d / %d (%.2f%%)",
                  mapped, total, 100 * mapped / total))
}

# === 9. Save output ===
write_tsv(junction_norm, "junctions_normalized_to_RPKM.tsv")
message("\n✅ Normalized junction file written to 'junctions_normalized_to_RPKM.tsv'")

nano run_normalize_junctions_to_RPKM.sh

#!/bin/bash
#SBATCH --job-name=jxn_norm
#SBATCH --output=jxn_norm_%j.out
#SBATCH --error=jxn_norm_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-type=ALL

# === Environment setup ===
module --force purge
module load R-bundle-Bioconductor/3.19-foss-2022b-R-4.4.1

echo "======================================="
echo " Job: Junction normalization to RPKM"
echo " Started: $(date)"
echo " Node: $(hostname)"
echo "======================================="

# === Check input files ===
echo "Checking input files..."
for f in merged_junction_counts.tsv gene_counts.txt junction_to_gene.tsv junctions_flagged.bed; do
    if [ ! -f "$f" ]; then
        echo "❌ ERROR: Missing required file: $f"
        echo "Aborting job."
        exit 1
    fi
done
echo "✅ All input files present."

# === Run the R script ===
echo "Running normalization script..."
Rscript normalize_junctions_to_RPKM.R

# === Completion message ===
status=$?
if [ $status -eq 0 ]; then
    echo "✅ Normalization completed successfully!"
else
    echo "❌ R script exited with an error (exit code $status)"
fi

echo "Finished: $(date)"
echo "======================================="
```
junction_to_gene <- read.delim("junction_to_gene.tsv", header=TRUE)
counts <- read.delim("gene_counts.txt", header=TRUE, comment.char="#")

# Check name overlap
cat("junction_to_gene IDs:", nrow(junction_to_gene), "\n")
cat("counts IDs:", nrow(counts), "\n")

# Check format
head(junction_to_gene$junction_id)
head(counts$Geneid)

# Check actual overlap
sum(counts$Geneid %in% junction_to_gene$junction_id)
