#!/bin/bash -l
#SBATCH --job-name=featurecounts
#SBATCH --output=/scratch/leuven/362/vsc36201/RNA-seq/logs/counts_%j.log
#SBATCH --error=/scratch/leuven/362/vsc36201/RNA-seq/logs/counts_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=16GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH -A lp_svbelleghem
#SBATCH --cluster=genius

echo "Starting featureCounts job at $(date)"

conda activate featurecounts

ANNOTATION="/scratch/leuven/362/vsc36201/RNA-seq/annotation/braker.gtf"
OUTPUT_DIR="/scratch/leuven/362/vsc36201/RNA-seq/counts"
ALIGNMENT_DIR="/scratch/leuven/362/vsc36201/RNA-seq/alignment"

mkdir -p "$OUTPUT_DIR"

featureCounts -T 4 -p --countReadPairs -Q 10 -g gene_id -a "$ANNOTATION" -o "$OUTPUT_DIR/pogonus.counts" "$ALIGNMENT_DIR"/*.filtered.sorted.bam

echo "featureCounts job finished at $(date)"
