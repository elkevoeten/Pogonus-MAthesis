#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name=star_genome_index
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8       
#SBATCH --time=12:00:00         
#SBATCH --mem=32G               
#SBATCH --output=/scratch/leuven/362/vsc36201/RNA-seq/logs/genome_index_%j.out
#SBATCH --error=/scratch/leuven/362/vsc36201/RNA-seq/logs/genome_index_%j.err
#SBATCH -A lp_svbelleghem      

# Load STAR module (if needed)
module load STAR/2.7.3a-GCCcore-6.4.0

# Create output directory for the genome index
mkdir -p /scratch/leuven/362/vsc36201/RNA-seq/genome_index

# Run STAR genome index generation
STAR --runMode genomeGenerate \
     --genomeDir /scratch/leuven/362/vsc36201/RNA-seq/genome_index \
     --genomeFastaFiles /scratch/leuven/362/vsc36201/RNA-seq/sorted_prim_dud.fasta \
     --runThreadN 8 \
     --sjdbGTFfile /scratch/leuven/362/vsc36201/RNA-seq/annotation/braker.gtf \  # Optional: GTF file for known gene annotations
     --sjdbOverhang 100  # This is generally set to (ReadLength - 1), e.g., 100 for 101 bp reads

echo "STAR genome index generation completed."
