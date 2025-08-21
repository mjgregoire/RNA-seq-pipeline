# RNA-seq Pipeline Master Script 
**Purpose: Reproducible workflow for RNA-seq analysis. This project is currently in progress and will be updated regularly from Agust 2025 through and to November 2025**

## Log in to HPC and load conda
Conda environment: rnaseq_tools (Packages: sra-tools, entrez-direct, fastqc, multiqc, samtools, etc...)

`module load miniconda`
`conda activate rnaseq_tools`

This was installed via: 

`conda install -c bioconda {package}` `conda create -n{name} -c bioconda {packages separated by spaces}`

You can check what packages are you the environment with `conda list`.

With the environment loaded you can install things needed in it later with the `conda install -c bioconda {package}` code

## Set up/go to working directory
`mkdir {path to your directory here, make an RNA seq folder and then sub folders for each project}`
`cd {path to your directory}`

### For downloading FASTQ from SRA projects or ENA projects with SRA tools, make accession list first
Get the SRR accession list from SRA Project (SRP). But you *don't want to use the login node* on HPC b/c will take up too much compute power and others can't do stuff, so request some time in a cluster to do manual work: `salloc --mem=1G --time=1:00:00`

Use the following code to get SRR accession list: `esearch -db sra -query {SRA name here} | efetch -format runinfo | cut -d ',' -f 1 | grep SRR > Run_Acc_List.txt`

Or for European files: `-qO - {link to ENA e.g.: "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=ERP131847&result=read_run&fields=run_accession"} | tail -n + 2 > Run_Acc_List.txt`

### Run FASTQ download script
This bash script will convert .sra to .fastq files if applicable, and downloads and compresses fastq files.
Create the script in your folder with: `nano fastq_array.sh`

```
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
cd ~{path to your directory here}
SRR=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" SRR_Acc_List.txt) 
echo "[$(date)] Downloading $SRR"
fasterq-dump --split-files --threads 4 "$SRR" 
gzip ${SRR}_1.fastq ${SRR}_2.fastq
```

Run the script with the following commmand: `sbatch fastq_array.sh`
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

To open the multiqc .html report you need to transfer the file to your local downloads folder. To do this you need to go to your local terminal. If you have not already set up an ssh key, do the following: 

`ssh-keygen`
`cat ~/YaleSSHkey.pub` and upload the key to Yale: https://sshkeys.ycrc.yale.edu/ 

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

On your terminal type the following to SSH into the HPC and download the file locally (this downloads to your Downloads folder): 
`scp -i ~/YaleSSHkey {your ID}@transfer-mccleary.ycrc.yale.edu:{path to your accesion here}/fastqs/fastqc_results/multiqc_report.html \
~/Downloads/multiqc_report_$(basename $(dirname $(dirname $(dirname /{path to your accession here}/fastqs/fastqc_results/multiqc_report.html)))).html`

## Trimming FASTQ files
To trim poor sequence quality or adapters on FASTQ files make and run the following bash script. 
`nano trim_fastp.sh`
`sbatch trim_fastp.sh`

```
#!/bin/bash
#SBATCH --job-name=fastp_trim
#SBATCH --output=fastp_trim.out
#SBATCH --error=fastp_trim.err
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL 
# Load modules or activate conda
module purge
module load miniconda
source activate rnaseq_tools  
# Define directories
IN_DIR=~{path to your directory here}
OUT_DIR=${IN_DIR}/trimmed_fastq
mkdir -p "$OUT_DIR"
# Trimming loop
for R1 in "$IN_DIR"/*_1.fastq.gz; do
    BASE=$(basename "$R1" _1.fastq.gz)
    R2="${IN_DIR}/${BASE}_2.fastq.gz"

    echo "Trimming $BASE..."

    # Estimate read length from the first read in the R1 file
    # zcat prints the gzipped file, head -n 2 gets the first sequence and its header
    READ_LEN=$(zcat "$R1" | head -n 2 | tail -n 1 | wc -c)
    READ_LEN=$((READ_LEN - 1))  # subtract newline

    # Set minimum length based on read length
    if [ "$READ_LEN" -le 50 ]; then
        MIN_LEN=20
    else
        MIN_LEN=50
    fi

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
      --html "$OUT_DIR/${BASE}_fastp.html" \
      --json "$OUT_DIR/${BASE}_fastp.json"

done
```
Check the script is working while running with:
```
squeue --me
tail -f fastp_trim.out
wc -l fastp_trim.out
```
Now you can re-Run FASTQC using the same script as before and check that the sequences have better quality now:

```
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
    cd ~{your path here}
    # Get the nth fastq.gz file
    FILE=$(ls *.fastq.gz | sort | sed -n "$((SLURM_ARRAY_TASK_ID+1))p")
    echo "Running FastQC on $FILE"
    fastqc "$FILE" --outdir ~{your path here/fastqc_results}

#compile all fastqc files into one
multiqc .
#transfer fastqc file to local downloads folder to check the html
#go to local terminal
#scp -i ~/YaleSSHkey {your ID}@transfer-mccleary.ycrc.yale.edu:{path to your accesion here}/fastqs/fastqc_results/multiqc_report.html \
~/Downloads/multiqc_report_$reRun(basename $(dirname $(dirname $(dirname /{path to your accession here}/fastqs/fastqc_results/multiqc_report.html)))).html
```


# Next steps:
# -------------------------------
# - align with STAR or HISAT2
# - quantify with featureCounts or Salmon
# - downstream DE analysis with DESeq2 / edgeR

echo "[$(date)] RNA-seq data download pipeline complete."
