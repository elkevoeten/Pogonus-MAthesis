cat("\014")   #clear console
rm(list=ls()) #clear environment

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("GenomicRanges")
# BiocManager::install("rtracklayer")

library(rstudioapi)
main_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(main_path))
getwd()

############################################################################################################
# script to extract repeat content of the dudzele reference genome(based on Fst calculations Maria's MA 
# thesis)
############################################################################################################
library(GenomicRanges)
library(rtracklayer)
library(dplyr)

gff_file <- "../data/sorted_prim_dud.fasta.out.gff.gz"

gff_data <- read.table(gff_file, sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(gff_data) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

hist((gff_data$end-gff_data$start), breaks = 100, xlim = c(0,340000),ylim = c(0,200),main=NA,xlab="Repeat length (bp)")

# Calculate the histogram data without plotting it
hist_data <- hist(gff_data$end - gff_data$start, breaks = 100, plot = FALSE)

# Generate a color gradient from light pink to dark red
color_gradient <- colorRampPalette(c("lightpink1", "lightpink4"))(length(hist_data$counts))

# Plot the histogram with gradient colors
plot(hist_data, col = color_gradient, border = "black", xlim = c(0, 340000), ylim = c(0, 200),
     main = NA, xlab = "Repeat length (bp)", ylab = "Frequency")

# Create GRanges object
gr <- GRanges(seqnames = gff_data$seqname,
              ranges = IRanges(start = gff_data$start, end = gff_data$end),
              strand = gff_data$strand)

# Get unique scaffolds
scaffolds <- unique(seqnames(gr))
scaffolds <- unique(scaffolds[grep("RagTag", scaffolds)])[c(1:10)]

result_list <- lapply(scaffolds, function(scaffold) {
  # Subset the GRanges object for the current scaffold
  gr_scaffold <- gr[seqnames(gr) == scaffold]
  
  # Determine the length of the scaffold
  scaffold_length <- max(end(gr_scaffold))
  
  # Create sliding windows for the scaffold
  windows <- slidingWindows(GRanges(seqnames = scaffold, ranges = IRanges(start = 1, end = scaffold_length)), 
                            width = 100000, step = 100000)[[1]]
  
  # Calculate the amount covered by features in each window
  coverage <- sapply(seq_along(windows), function(i) {
    window <- windows[i]
    overlaps <- pintersect(rep(window, length(gr_scaffold)), gr_scaffold)
    sum(width(overlaps))
  })
  
  # Create a data frame for the results
  data.frame(
    seqnames = scaffold,
    start = start(windows),
    end = end(windows),
    coverage = coverage
  )
})

# Combine the results for all scaffolds
result <- do.call(rbind, result_list)

# Print the result
print(result)

result$mid <- (result$start+result$end)/2
result$perc <- result$coverage/100000

resultS <- subset(result, result$seqnames == 'CM008230.1_RagTag')
head(resultS)

plot(resultS$start, resultS$perc, pch=19, cex = 0.5)

scaffolds <- c("CM008230.1_RagTag", "CM008231.1_RagTag", "CM008233.1_RagTag", "CM008234.1_RagTag", "CM008235.1_RagTag", "CM008236.1_RagTag", "CM008237.1_RagTag", "CM008238.1_RagTag","CM008239.1_RagTag","CM008240.1_RagTag")

for (scaffold in scaffolds) {
  result_scaffolds <- subset(result, result$seqnames == scaffold)
  
  plot(result_scaffolds$start, result_scaffolds$perc, pch=19, cex = 0.6, col="orange4",main = scaffold, xlab = "Chromosome Position",ylab="Repeat Content")
}

result$chromosome <- dplyr::case_when(
  result$seqnames == "CM008230.1_RagTag" ~ "Chr1",
  result$seqnames == "CM008233.1_RagTag" ~ "Chr2",
  result$seqnames == "CM008234.1_RagTag" ~ "Chr3",
  result$seqnames == "CM008235.1_RagTag" ~ "Chr4",
  result$seqnames == "CM008236.1_RagTag" ~ "Chr5",
  result$seqnames == "CM008237.1_RagTag" ~ "Chr6",
  result$seqnames == "CM008238.1_RagTag" ~ "Chr7",
  result$seqnames == "CM008239.1_RagTag" ~ "Chr8",
  result$seqnames == "CM008240.1_RagTag" ~ "Chr9",
  result$seqnames == "CM008231.1_RagTag" ~ "Chr10",
  TRUE ~ NA_character_
)

# write.table(result, file = "intermediate-output/repeat-content.tsv", sep = "\t", quote = FALSE, row.names = T)
