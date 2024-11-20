#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name 2_msmc2
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=12
#SBATCH --time=72:00:00
#SBATCH --output=2_msmc2.%j.out
#SBATCH --error=2_msmc2.%j.err
#SBATCH --array=1-10
#SBATCH -A lp_svbelleghem

source /data/leuven/362/vsc36201/miniconda3/etc/profile.d/conda.sh

module load BCFtools/1.15.1-GCC-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0
module load SciPy-bundle/2023.07-gfbf-2023a

#You first need to create an environment called msmc2 and install msmc2 in that env:
#conda create -n msmc2 msmc2
conda activate msmc2

# Set PYTHONPATH to include msmc-tools directory
export PYTHONPATH=/data/leuven/362/vsc36201/msmc-tools:$PYTHONPATH

#you need to submit with sbatch -a 1-10 (or however many samples in the array) so that all individuals run at the same time 
indID=$((SLURM_ARRAY_TASK_ID -1))

REF=sorted_prim_dud.fasta
IN=/scratch/leuven/362/vsc36201/MSMC_analysis/msmc-in
OUT=/scratch/leuven/362/vsc36201/MSMC_analysis/msmc-out
OUT2=/scratch/leuven/362/vsc36201/MSMC_analysis/msmc-final

bamCaller="/data/leuven/362/vsc36201/msmc-tools/bamCaller.py"
generate_multihetsep="/data/leuven/362/vsc36201/msmc-tools/generate_multihetsep.py"

SAMPLE=(GC129388 GC129395 GC129400 GC129406 GC129413 GC129417 GC136116 GC136110 GC136117 GC136123)
samtools index -M $IN/$(echo "${SAMPLE[indID]}").dudPrim.filtered.sorted.nd.bam $IN/$(echo "${SAMPLE[indID]}").dudPrim.filtered.sorted.nd.bam.bai

# Loop through each scaffold (SCAF)
for SCAF in CM008230.1_RagTag CM008231.1_RagTag CM008233.1_RagTag CM008234.1_RagTag CM008235.1_RagTag CM008236.1_RagTag CM008237.1_RagTag CM008238.1_RagTag CM008239.1_RagTag CM008240.1_RagTag
do
    # Call variants for the current scaffold
    bcftools mpileup -q 20 -Q 20 -C 50 -r $SCAF -f $REF $IN/"${SAMPLE[indID]}".dudPrim.filtered.sorted.nd.bam | \
    bcftools call -c - | \
    python "$bamCaller" 30 $OUT/"${SAMPLE[indID]}_$SCAF.mask.bed.gz" | \
    gzip -c > $OUT/"${SAMPLE[indID]}_$SCAF.vcf.gz"

    # Generate multi-heterozygosity separation
    python "$generate_multihetsep" --mask=$OUT/"${SAMPLE[indID]}_$SCAF.mask.bed.gz" $OUT/"${SAMPLE[indID]}_$SCAF.vcf.gz" > $OUT/"${SAMPLE[indID]}_$SCAF.txt"
done

# Collect all generated text file paths for further processing
COMMAND=""
for SCAF in CM008230.1_RagTag CM008231.1_RagTag CM008233.1_RagTag CM008234.1_RagTag CM008235.1_RagTag CM008236.1_RagTag CM008237.1_RagTag CM008238.1_RagTag CM008239.1_RagTag CM008240.1_RagTag
do
    COMMAND="$COMMAND $OUT/${SAMPLE[indID]}_$SCAF.txt"
done

#this script was written for msmc v2.1.4 
msmc2_Linux -t 12 -o $OUT2/$(echo "${SAMPLE[indID]}") $COMMAND

#Sample names
for FINAL in GC129388 GC129395 GC129400 GC129406 GC129413 GC129417 GC136116 GC136110 GC136117 GC136123
do
python MSMC_plotInput.py -I $OUT2/$FINAL.final.txt -u 2.1e-09 -g 1 > $OUT2/$FINAL.final.Rin.txt
done
