### For downloading FASTQ from SRA projects or ENA projects with SRA tools, make accession list first
Get the SRR accession list from SRA Project (SRP). But you *don't want to use the login node* on HPC b/c will take up too much compute power and others can't do stuff, so request some time in a cluster to do manual work: `salloc --mem=1G --time=1:00:00`

Use the following code to get SRR accession list: `esearch -db sra -query {SRA name here} | efetch -format runinfo | cut -d ',' -f 1 | grep SRR > Run_Acc_List.txt`

Or for European files: `wget -qO - {link to ENA e.g.: "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=ERP131847&result=read_run&fields=run_accession"} | tail -n + 2 > Run_Acc_List.txt`

### Run FASTQ download script
This bash script will convert .sra to .fastq files if applicable, and downloads and compresses fastq files.
Create the script in your folder with: `nano download_fastq.sh`

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

Run the script with the following commmand: `sbatch --array=0-$(( $(wc -l < Run_Acc_List.txt) - 1 )) download_fastq.sh`
Using sbatch on bash scripts will send them to an available cluster and they will run in the background until completed, the time runs out, or they run into an error. You can at this point run other scripts or leave the command line and your script should still run. 

**Verify output from download**
Check that the bash script is working with: `squeue --me`
Or with the output files using: `tail -f ~{your file path here}/logs/*.out`

Run ls or check file counts: `ls -lh fastq_files | grep '.fastq.gz'`
`wc -l SRR_Acc_List.txt`
`ls *_1.fastq.gz | wc -l`

`ls *_1.fastq.gz | sed 's/_1.fastq.gz//' | sort > downloaded.txt`
`comm -23 <(sort SRR_Acc_List.txt) downloaded.txt > not_downloaded.txt`
