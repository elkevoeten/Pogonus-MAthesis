#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name=fastqc
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8 
#SBATCH --time=24:00:00 
#SBATCH --output=/scratch/leuven/362/vsc36201/RNA-seq/logs/fastqc.%j.out
#SBATCH --error=/scratch/leuven/362/vsc36201/RNA-seq/logs/fastqc.%j.err
#SBATCH -A lp_svbelleghem

# Load FastQC module
module load FastQC/0.11.9-Java-11  

# Create output directory 
mkdir -p /scratch/leuven/362/vsc36201/RNA-seq/fastqc

# Run FastQC on all .fq.gz files 
for sample in /scratch/leuven/362/vsc36201/RNA-seq/fastq-files/*.fq.gz; do
    echo "Processing $sample..."
    fastqc -o /scratch/leuven/362/vsc36201/RNA-seq/fastqc "$sample" &
done

wait

echo "FastQC analysis completed!"
