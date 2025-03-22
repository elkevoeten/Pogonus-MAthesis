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

# Load modules
module load cluster/genius/batch
module load SAMtools/1.18-GCC-12.3.0
module load libdeflate/1.19-GCCcore-13.2.0

# Define BAM-files directory
bam_dir="/scratch/leuven/362/vsc36201/RNA-seq/alignment"

# Loop through all BAM-files for filtering and sorting
for bamfile in "$bam_dir"/*.Aligned.sortedByCoord.out.bam; do
    sample_name=$(basename "$bamfile" .Aligned.sortedByCoord.out.bam)
    filtered_bam="${bam_dir}/${sample_name}.filtered.bam"
    sorted_bam="${bam_dir}/${sample_name}.filtered.sorted.bam"
    samtools view -f 0x02 -q 30 -b "$bamfile" > "$filtered_bam"
    samtools sort -o "$sorted_bam" "$filtered_bam"
    echo "Filtered and sorted $bamfile -> $sorted_bam"
done
