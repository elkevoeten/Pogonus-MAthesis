#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name msmc2_103
#SBATCH --nodes=5 
#SBATCH --ntasks-per-node=21
#SBATCH --time=72:00:00
#SBATCH --output=/scratch/leuven/362/vsc36201/MSMC_analysis/logs/msmc2.%j.out
#SBATCH --error=/scratch/leuven/362/vsc36201/MSMC_analysis/logs/msmc2.%j.err
#SBATCH --array=1-103
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

SAMPLE=(GC129388 GC129389 GC129390 GC129391 GC129392 GC129393 \
	GC129394 GC129395 GC129396 GC129397 GC129398 GC129399 \
	GC129400 GC129401 GC129402 GC129403 GC129404 GC129405 \
	GC129406 GC129407 GC129408 GC129409 GC129410 GC129411 \
	GC129412 GC129413 GC129414 GC129415 GC129416 GC129417 \
	GC129418 GC129419 GC129420 GC129421 GC129422 GC129423 \
	GC129424 GC129425 GC129426 GC129427 GC129428 GC129429 \
	GC129430 GC129431 GC129432 GC129433 GC129434 GC129435 \
	GC129437 GC129438 GC129439 GC129440 GC136078 GC136079 \
	GC136080 GC136081 GC136082 GC136083 GC136084 GC136085 \
	GC136086 GC136087 GC136088 GC136089 GC136090 GC136091 \
	GC136092 GC136093 GC136094 GC136095 GC136096 GC136097 \
	GC136098 GC136099 GC136100 GC136101 GC136102 GC136103 \
	GC136104 GC136105 GC136106 GC136107 GC136108 GC136109 \
	GC136110 GC136111 GC136112 GC136113 GC136114 GC136115 \
	GC136116 GC136117 GC136118 GC136119 GC136120 GC136121 \
	GC136122 GC136123 GC136124 GC136125 GC136126 GC136127 \
	GC136128)

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
for FINAL in "${SAMPLE[@]}"
do
python MSMC_plotInput.py -I $OUT2/$FINAL.final.txt -u 2.1e-09 -g 1 > $OUT2/$FINAL.final.Rin.txt
done
