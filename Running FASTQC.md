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
