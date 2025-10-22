# ALIGNING FASTQ TO REFERENCE GENOME USING STAR
## Make a directory to store your references
```
mkdir -p /gpfs/gibbs/pi/guo/mg2684/reference/gencode
cd /gpfs/gibbs/pi/guo/mg2684/reference/gencode
```

## Genome FASTA (GRCh38 primary assembly)
```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
```

## Annotation GTF
```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz
gunzip gencode.v43.annotation.gtf.gz
```

## Make reference index
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

## Align
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
`sbatch STAR_align_array.sh`
