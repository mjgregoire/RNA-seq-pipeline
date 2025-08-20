# RNA-seq Pipeline Master Script 
**Purpose: Reproducible workflow for RNA-seq analysis.**

## Log in to HPC and load conda
Conda environment: rnaseq_tools (Packages: sra-tools, entrez-direct, fastqc, multiqc, samtools, etc...)

` module load miniconda `
` conda activate rnaseq_tools `

This was installed via: 

` conda install -c bioconda {package} ` ` conda create -n{name} -c bioconda {packages separated by spaces} `

You can check what packages are you the environment with `conda list`.

With the environment loaded you can install things needed in it later with the ` conda install -c bioconda {package} ` code

## Set up/go to working directory
mkdir {path to your directory here, e.g.: /gpfs/gibbs/pi/guo/mg2684/RNAseq/ProjectName
cd {path to your directory}

### For downloading SRA projects or ENA projects with SRA tools
Get the SRR accession list from SRA Project (SRP)
You don't want to use the login node on HPC b/c will take up too much compute power so others can't do stuff, so request some time in a cluster to do manual work: ` salloc --mem=1G --time=1:00:00 `
Use the following code to get SRR accession list: ` esearch -db sra -query {SRA name here} | efetch -format runinfo | cut -d ',' -f 1 | grep SRR > SRR_Acc_List.txt `
Or for European files: ` -qO - {link to ENA e.g.: "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=ERP131847&result=read_run&fields=run_accession"} | tail -n + 2 > ENA_Run_List.txt


### Run FASTQ download script
# Script downloads all SRR files and compresses them and converts .sra to .fastq files
echo "Starting FASTQ download with fasterq-dump..."
sbatch fastq_array.sh

    #!/bin/bash
    ###Run on day partition###
    #SBATCH -p day
    ###name job###
    #SBATCH --job-name=fastq_array
    ###request 1 GB memory###
    #SBATCH --mem-per-cpu=1G
    ###request 6 hours worth of time### #or might need more time depending on size and number of files
    #SBATCH -t 6:00:00
    ###create email trail###
    #SBATCH --mail-type=ALL
    #SBATCH --output=logs/fastq_%A_%a.out
    #SBATCH --error=logs/fastq_%A_%a.err
    #SBATCH --array=0-24     # Adjust if you have more or fewer SRRs
    module load miniconda
    conda activate rnaseq_tools
    cd ~/palmer_scratch/rnaseq_download
    SRR=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" SRR_Acc_List.txt)
    echo "[$(date)] Downloading $SRR"
    fasterq-dump --split-files --threads 4 "$SRR"
    gzip ${SRR}_1.fastq ${SRR}_2.fastq

# -------------------------------
# Step 4: Verify output
# -------------------------------
# Optional: run ls or check file counts
ls -lh fastq_files | grep '.fastq.gz'
echo "Download complete: $(ls fastq_files/*.fastq.gz | wc -l) fastq files found."

#or?
squeue --me
tail -f ~/palmer_scratch/RNAseq_download/logs/*.out

wc -l SRR_Acc_List.txt
ls *_1.fastq.gz | wc -l

ls *_1.fastq.gz | sed 's/_1.fastq.gz//' | sort > downloaded.txt
comm -23 <(sort SRR_Acc_List.txt) downloaded.txt > not_downloaded.txt

# -------------------------------
# Step 5: FASTQC
# -------------------------------
# run fastqc as as bash script
echo "Starting FASTQC..."
mkdir fastqc_results
sbatch fastqc.sh

    #!/bin/bash
    #SBATCH -p day
    #SBATCH --job-name=fastqc
    #SBATCH --mail-type=ALL
    #SBATCH --mem=4G
    #SBATCH -t 2:00:00
    #SBATCH --array=0-49
    #SBATCH --output=logs/fastqc_%A_%a.out
    #SBATCH --error=logs/fastqc_%A_%a.err
    module load miniconda
    conda activate rnaseq_tools
    cd ~/palmer_scratch/RNAseq_download
    # Get the nth fastq.gz file
    FILE=$(ls *.fastq.gz | sort | sed -n "$((SLURM_ARRAY_TASK_ID+1))p")
    echo "Running FastQC on $FILE"
    fastqc "$FILE" --outdir ~/palmer_scratch/RNAseq_download/fastqc_results

#compile all fastqc files into one
multiqc .
#transfer fastqc file to local downloads folder to check the html
#go to local terminal
#if ssh not set up:
#ssh-keygen
#cat ~/YaleSSHkey.pub #upload to yale: https://sshkeys.ycrc.yale.edu/ 
#ssh -i ~/YaleSSHkey mg2684@mccleary.ycrc.yale.edu
#nano ~/.ssh/config: Host mccleary
    HostName mccleary.ycrc.yale.edu
    User mg2684
    IdentityFile ~/YaleSSHkey
#chmod 600 ~/.ssh/config
#chmod 600 ~/YaleSSHkey
#chmod 700 ~/.ssh
scp -i ~/YaleSSHkey mg2684@transfer-mccleary.ycrc.yale.edu:~/palmer_scratch/RNAseq_download/fastqc_results/multiqc_report.html ~/Downloads/

# -------------------------------
# Step 6: Trim and QC
# -------------------------------
echo "Staring trimming..."
sbatch trim_fastp.sbatch.
    #!/bin/bash
    #SBATCH --job-name=fastp_trim
    #SBATCH --output=fastp_trim.out
    #SBATCH --error=fastp_trim.err
    #SBATCH --time=06:00:00
    #SBATCH --mem=16G
    #SBATCH --cpus-per-task=4
    #SBATCH --mail-user=michelle.gregoire@yale.edu  
    #SBATCH --mail-type=END,FAIL

    # Load modules or activate conda
    module purge
    module load miniconda
    source activate your_fastp_env  # Replace with your conda env name that has fastp installed

    # Define directories
    IN_DIR=~/palmer_scratch/RNAseq_download
    OUT_DIR=${IN_DIR}/trimmed_fastq
    mkdir -p "$OUT_DIR"

    # Trimming loop
    for R1 in "$IN_DIR"/*_1.fastq.gz; do
        BASE=$(basename "$R1" _1.fastq.gz)
        R2="${IN_DIR}/${BASE}_2.fastq.gz"

        echo "Trimming $BASE..."

        fastp \
          -i "$R1" \
          -I "$R2" \
          -o "$OUT_DIR/${BASE}_R1.trimmed.fastq.gz" \
          -O "$OUT_DIR/${BASE}_R2.trimmed.fastq.gz" \
          --detect_adapter_for_pe \
          --cut_front \
          --cut_tail \
          --cut_window_size 4 \
          --cut_mean_quality 20 \
          --length_required 20 \
          --thread 4 \
          --html "$OUT_DIR/${BASE}_fastp.html" \
          --json "$OUT_DIR/${BASE}_fastp.json"
    done

    echo "All trimming done."

#check it is working while running with:
squeue --me
tail -f fastp_trim.out
wc -l fastp_trim.out

# -------------------------------
# Step 7: Re-Run FASTQC
# -------------------------------
# run fastqc as as bash script
echo "Starting FASTQC..."
mkdir fastqc_results
sbatch fastqc.sh

    #!/bin/bash
    #SBATCH -p day
    #SBATCH --job-name=fastqc
    #SBATCH --mail-type=ALL
    #SBATCH --mem=4G
    #SBATCH -t 2:00:00
    #SBATCH --array=0-49
    #SBATCH --output=logs/fastqc_%A_%a.out
    #SBATCH --error=logs/fastqc_%A_%a.err
    module load miniconda
    conda activate rnaseq_tools
    cd ~/palmer_scratch/RNAseq_download
    # Get the nth fastq.gz file
    FILE=$(ls *.fastq.gz | sort | sed -n "$((SLURM_ARRAY_TASK_ID+1))p")
    echo "Running FastQC on $FILE"
    fastqc "$FILE" --outdir ~/palmer_scratch/RNAseq_download/fastqc_results

#compile all fastqc files into one
multiqc .
#transfer fastqc file to local downloads folder to check the html
#go to local terminal
#ssh -i ~/YaleSSHkey mg2684@mccleary.ycrc.yale.edu
scp -i ~/YaleSSHkey mg2684@transfer-mccleary.ycrc.yale.edu:~/palmer_scratch/RNAseq_download/trimmed_fastq/fastqc_results/multiqc_report.html ~/Downloads/

# -------------------------------
# Next steps:
# -------------------------------
# - align with STAR or HISAT2
# - quantify with featureCounts or Salmon
# - downstream DE analysis with DESeq2 / edgeR

echo "[$(date)] RNA-seq data download pipeline complete."
