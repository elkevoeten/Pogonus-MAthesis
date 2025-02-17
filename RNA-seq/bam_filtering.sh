#!/bin/bash
#SBATCH --job-name=bam_filtering
#SBATCH --cluster=genius
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8       
#SBATCH --time=48:00:00         
#SBATCH --mem=32G               
#SBATCH --output=/scratch/leuven/362/vsc36201/RNA-seq/logs/filtering_%j.out
#SBATCH --error=/scratch/leuven/362/vsc36201/RNA-seq/logs/filtering_%j.err
#SBATCH -A lp_svbelleghem      

# Load the SAMtools module
module load cluster/genius/batch
module load SAMtools/1.18-GCC-12.3.0
module load libdeflate/1.19-GCCcore-13.2.0

# Define the directory containing the BAM files
bam_dir="/scratch/leuven/362/vsc36201/RNA-seq/alignment"

# Loop through all BAM files in the alignment directory
for bamfile in "$bam_dir"/*.Aligned.sortedByCoord.out.bam; do
    # Extract the sample name by removing the file extension
    sample_name=$(basename "$bamfile" .Aligned.sortedByCoord.out.bam)
    
    # Construct the output filtered BAM file path
    output_bam="${bam_dir}/${sample_name}.sorted.filtered.bam"
    
    # Run samtools view command to filter the BAM files
    samtools view -f 0x02 -q 30 -b "$bamfile" > "$output_bam"
    
    # Print a message when done processing each file
    echo "Filtered $bamfile -> $output_bam"
done

echo "Filtering complete for all BAM files."
