# RNA-Seq Analysis Pipeline

This repository contains scripts for processing RNA-Seq data. The scripts were executed in the following order:

#### 1. Quality check of `.fq.qz`-files
Run FastQC to assess the quality of raw sequencing reads: `fastqc.sh`

#### 2. Indexing reference genome `sorted_prim_dud.fasta`
Prepare the reference genome for alignment using STAR: `star_reference_index.sh`

#### 3. Alignment to the reference genome using STAR
Map the reads to the reference genome: `star_alignment.sh`

#### 4. Filter and sort BAM-files
Process the aligned reads: `bam_filtering.sh`

#### 5. Check the alignment stats
```sh
for bamfile in alignment/*.filtered.sorted.bam; do
    echo "Processing $bamfile..."
    samtools stats "$bamfile" > "${bamfile%.filtered.sorted.bam}.stats"
done
echo "All files have been processed!"
```

#### 6. Index BAM-files
Index the filtered and sorted BAM files:
```sh
for bam_file in *.filtered.sorted.bam; do
    echo "Indexing $bam_file"
    samtools index "$bam_file"
    echo "Finished $bam_file"
done
```

#### 7. FeatureCounts Analysis
Run featureCounts to quantify gene expression: `featurecounts.sh`

Format the counts file:
```sh
grep -v "^#" pogonus.counts | sed '/^Status/d' \
    | sed 's|/scratch/leuven/362/vsc36201/RNA-seq/alignment/||g' \
    | sed 's/.filtered.sorted.bam//g' \
    > pogonus.counts.formatted
```

#### 8. Differential Expression Analysis with DESeq2
Run DESeq2 analysis in R: `PCBAR_Deseq2.R`

##### Extract Up-Regulated Genes
Write up-regulated genes (log2FoldChange > 1, p-value < 0.05) to a file:
```sh
awk '$3 != "NA" && $3 > 1 && $7 < 0.05 {print}' deseq.results.tsv | sort -k3,3nr > PCBAR-up-logFC1-p05.txt
```
Extract gene names:
```sh
awk '{print $1}' PCBAR-up-logFC1-p05.txt > PCBAR-up-logFC1-p05-geneid.txt
```

##### Extract Down-Regulated Genes
Write down-regulated genes (log2FoldChange < -1, p-value < 0.05) to a file:
```sh
awk '$3 != "NA" && $3 < -1 && $7 < 0.05 {print}' deseq.results.tsv | sort -k3,3n > PCBAR-down-logFC1-p05.txt
```
Extract gene names:
```sh
awk '{print $1}' PCBAR-down-logFC1-p05.txt > PCBAR-down-logFC1-p05-geneid.txt
```

##### Count Number of Differentially Expressed Genes
```sh
wc -l PCBAR-up-logFC1-p05.txt
wc -l PCBAR-down-logFC1-p05.txt
```

---

