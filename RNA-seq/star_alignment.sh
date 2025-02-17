#!/bin/bash -l
#SBATCH --job-name=star_alignment
#SBATCH --cluster=genius
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8       
#SBATCH --time=48:00:00        
#SBATCH --mem=32G              
#SBATCH --output=/scratch/leuven/362/vsc36201/RNA-seq/logs/star_%j.out
#SBATCH --error=/scratch/leuven/362/vsc36201/RNA-seq/logs/star_%j.err
#SBATCH -A lp_svbelleghem       

# Load STAR module 
module load STAR/2.7.3a-GCCcore-6.4.0

# Create output directory
mkdir -p /scratch/leuven/362/vsc36201/RNA-seq/alignment

# Loop through all paired-end FASTQ files
for sample in /scratch/leuven/362/vsc36201/RNA-seq/fastq-files/*_1.fq.gz; do
    sample_name=$(basename "$sample" "_1.fq.gz")
    fastq1="/scratch/leuven/362/vsc36201/RNA-seq/fastq-files/${sample_name}_1.fq.gz"
    fastq2="/scratch/leuven/362/vsc36201/RNA-seq/fastq-files/${sample_name}_2.fq.gz"
    
    echo "Processing sample: $sample_name"
    
    STAR \
        --genomeDir /scratch/leuven/362/vsc36201/RNA-seq/genome_index/ \
        --runThreadN 8 \
        --readFilesIn "$fastq1" "$fastq2" \
        --readFilesCommand zcat \
        --outFileNamePrefix "/scratch/leuven/362/vsc36201/RNA-seq/alignment/${sample_name}." \
        --outSAMtype BAM SortedByCoordinate

    echo "Finished processing: $sample_name"
done

echo "STAR alignment completed."
