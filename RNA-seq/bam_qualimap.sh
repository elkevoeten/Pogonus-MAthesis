#!/bin/bash
#SBATCH --job-name=bam_qualimap
#SBATCH --cluster=genius
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8       
#SBATCH --time=48:00:00         
#SBATCH --mem=32G               
#SBATCH --output=/scratch/leuven/362/vsc36201/RNA-seq/logs/qualimap_%j.out
#SBATCH --error=/scratch/leuven/362/vsc36201/RNA-seq/logs/qualimap_%j.err
#SBATCH -A lp_svbelleghem       

# Load necessary modules
module load cluster/genius/amd
module load Qualimap/2.3-foss-2022b-R-4.2.2

# Define paths
BAM_DIR="/scratch/leuven/362/vsc36201/RNA-seq/alignment"   # Change to your BAM file directory
OUT_DIR="/scratch/leuven/362/vsc36201/RNA-seq/qualimap"    # Change to your output directory

# Create output directory
mkdir -p "$OUT_DIR"

# Loop through BAM files and run Qualimap
for bam_file in "$BAM_DIR"/*.sorted.filtered.bam; do
    sample_name=$(basename "$bam_file" .sorted.filtered.bam)

    echo "Running Qualimap on $sample_name..."
    
    qualimap bamqc -bam "$bam_file" -outdir "$OUT_DIR" -outformat PDF:HTML
    
    echo "Finished processing $sample_name"
done

echo "Qualimap analysis complete!"
