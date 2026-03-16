cat("\014")   #clear console
rm(list=ls()) #clear environment

library(rstudioapi)
main_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(main_path))

############################################################################################################
# script to annotate the DESeq2.CG.resLFC.sig.DEG.tsv dataset to be able to plot the differentially 
# expressed genes along their genome position.
############################################################################################################

# These paths should exist as we're already at step 05, but you can run this for the sake of validation:
# output_path <- "output-figures"
# if (!dir.exists(output_path)) {
#   dir.create(output_path)
#   cat("Output path created:", output_path, "\n")
# } else {
#   cat("Output path already exists:", output_path, "\n")
# }
# 
# intermediate_output_path <- "intermediate-output"
# if (!dir.exists(intermediate_output_path)) {
#   dir.create(intermediate_output_path)
#   cat("Output path created:", intermediate_output_path, "\n")
# } else {
#   cat("Output path already exists:", intermediate_output_path, "\n")
# }

# LOAD DATA ---------------------------------------------------------------------------------------------------
CG.resLFC.sig.DEG <- read.delim("intermediate-output/DESeq2.CG.resLFC.sig.DEG.tsv", sep="\t")

CG.resLFC.sig.DEG.plot <- CG.resLFC.sig.DEG %>% 
  mutate(diff_expr = ifelse(log2FoldChange > 0, "upregulation", "downregulation"))

inversions <- read.delim("../data/crossings_inversions.tsv", header = TRUE, sep = "\t")
regions <- read.delim("../data/pogonus_regions.tsv", header = TRUE, sep = "\t")
fst_intervals <- read.delim("intermediate-output/Fst-intervals.tsv", sep = "\t")

# ANNOTATE DATAFRAME --------------------------------------------------------------------------------------------
# function to assign each gene to a region:
assign_location <- function(gene_chr, gene_start, inv_coords) {
  inv_chr <- inv_coords[inv_coords$Chr == gene_chr, ]
  in_inversion <- any(gene_start >= inv_chr$start & gene_start <= inv_chr$end)
  if (in_inversion) {
    return("inversion")
  } else {
    return("background")
  }
}

CG.resLFC.sig.DEG.plot$location <- mapply(assign_location,
                                       gene_chr = CG.resLFC.sig.DEG.plot$Chr,
                                       gene_start = CG.resLFC.sig.DEG.plot$Start,
                                       MoreArgs = list(inv_coords = inversions))
head(CG.resLFC.sig.DEG.plot)

# do the same for the regions (could be handy for boxplots)
assign_region <- function(gene_chr, gene_start, region_data) {
  matching_regions <- region_data[region_data$Chr == gene_chr, ]
  in_region <- matching_regions$start <= gene_start & gene_start <= matching_regions$end
  if (any(in_region)) {
    return(matching_regions$region[which(in_region)[1]])
  } else {
    return(NA)
  }
}

CG.resLFC.sig.DEG.plot$region <- mapply(assign_region,
                                     gene_chr = CG.resLFC.sig.DEG.plot$Chr,
                                     gene_start = CG.resLFC.sig.DEG.plot$Start,
                                     MoreArgs = list(region_data = regions))
str(CG.resLFC.sig.DEG.plot)
CG.resLFC.sig.DEG.plot <- CG.resLFC.sig.DEG.plot %>% mutate_if(is.character, as.factor) 
CG.resLFC.sig.DEG.plot$GeneID <- as.character(CG.resLFC.sig.DEG.plot$GeneID)
CG.resLFC.sig.DEG.plot$region <- factor(CG.resLFC.sig.DEG.plot$region, levels = regions$region)

CG.resLFC.sig.DEG.plot$chromosome <- dplyr::case_when(
  CG.resLFC.sig.DEG.plot$Chr == "CM008230.1_RagTag" ~ "Chr1",
  CG.resLFC.sig.DEG.plot$Chr == "CM008233.1_RagTag" ~ "Chr2",
  CG.resLFC.sig.DEG.plot$Chr == "CM008234.1_RagTag" ~ "Chr3",
  CG.resLFC.sig.DEG.plot$Chr == "CM008235.1_RagTag" ~ "Chr4",
  CG.resLFC.sig.DEG.plot$Chr == "CM008236.1_RagTag" ~ "Chr5",
  CG.resLFC.sig.DEG.plot$Chr == "CM008237.1_RagTag" ~ "Chr6",
  CG.resLFC.sig.DEG.plot$Chr == "CM008238.1_RagTag" ~ "Chr7",
  CG.resLFC.sig.DEG.plot$Chr == "CM008239.1_RagTag" ~ "Chr8",
  CG.resLFC.sig.DEG.plot$Chr == "CM008240.1_RagTag" ~ "Chr9",
  CG.resLFC.sig.DEG.plot$Chr == "CM008231.1_RagTag" ~ "Chr10",
  TRUE ~ NA_character_
)

CG.resLFC.sig.DEG.plot$chromosome <- factor(CG.resLFC.sig.DEG.plot$chromosome,
                                    levels = paste0("Chr", 1:10))

# write.table(CG.resLFC.sig.DEG.plot, file = "intermediate-output/DESeq2.CG.resLFC.sig.DEG.annotated.tsv", sep = "\t", quote = FALSE, row.names = T)