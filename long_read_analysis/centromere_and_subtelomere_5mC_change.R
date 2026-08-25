# Make plots showing 5mC change at centromeres and subtelomeres

# Load required packages and functions
library(methodical)
source("../helper_scripts/ggplot_functions.R")

# Load 5mC RSE and mask methylation of CpGs covered by < 10 reads or >= 350 reads
chlamy_5mc_rse = HDF5Array::loadHDF5SummarizedExperiment("CC2937_T2T_5mC_rse")
assay(chlamy_5mc_rse, 1)[assay(chlamy_5mc_rse, 2) < 10 | assay(chlamy_5mc_rse, 2) >= 350] = NA

### Plot 5mC change at centromeres

# Import locations of centromeres
centromeres_gr = rtracklayer::import.bed("../genome_files/CC2937_T2T.centromeres.bed")

# Get 5mC for all centromeres
centromere_methylation = summarizeRegionMethylation(chlamy_5mc_rse, genomic_regions = centromeres_gr)
centromere_methylation$chromosome = as.numeric(gsub("chromosome_", "", seqnames(GRanges(row.names(centromere_methylation)))))
centromere_methylation$chromosome = factor(centromere_methylation$chromosome, 
  levels = rev(as.character(as.numeric(gsub("chromosome_", "", seqnames(GRanges(row.names(centromere_methylation))))))))

# Add a column with the change in the deletion samples
centromere_methylation$deletion_change = (centromere_methylation$chlamy_7A + centromere_methylation$chlamy_18A)/2 - 
  (centromere_methylation$chlamy_8A + centromere_methylation$chlamy_15A)/2

# Make a barplot for centromere 5mC change
centromere_barplot = ggplot(centromere_methylation, aes(x = chromosome, y = deletion_change)) +
  geom_col(fill = "#68abb8", colour = "black")
centromere_barplot = customize_ggplot_theme(centromere_barplot, 
  title = NULL, xlab = "Chromosome", ylab = "Centromere Mean mCG Change") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) + 
  geom_hline(yintercept = 0) +
  coord_flip()
centromere_barplot
ggsave(plot = centromere_barplot, "plots/centromere_meth_change_barplot.pdf", width = 6, height = 9, device = cairo_pdf)

### Plot 5mC change at subtelomeres

# Import locations of subtelomers and indicate if they are on the p or q arm
subtelomeres_gr = rtracklayer::import.bed("../genome_files/CC2937_T2T.subtelomeres.bed")
subtelomeres_gr$arm = "q"
subtelomeres_gr[nearest(centromeres_gr, subtelomeres_gr)]$arm = "p"

# Get 5mC for all subtelomeres
subtelomere_methylation = summarizeRegionMethylation(chlamy_5mc_rse, genomic_regions = subtelomeres_gr)
subtelomere_methylation$chromosome = as.numeric(gsub("chromosome_", "", seqnames(GRanges(row.names(subtelomere_methylation)))))
subtelomere_methylation$arm = factor(subtelomeres_gr$arm, levels = c("q", "p"))
subtelomere_methylation$chromosome = factor(subtelomere_methylation$chromosome, 
  levels = rev(unique(as.character(as.numeric(gsub("chromosome_", "", seqnames(GRanges(row.names(subtelomere_methylation)))))))))

# Add a column with the change in the deletion samples
subtelomere_methylation$deletion_change = (subtelomere_methylation$chlamy_7A + subtelomere_methylation$chlamy_18A)/2 - 
  (subtelomere_methylation$chlamy_8A + subtelomere_methylation$chlamy_15A)/2

# Make a barplot for subtelomere 5mC change
subtelomere_barplot = ggplot(subtelomere_methylation, aes(x = chromosome, y = deletion_change, fill = arm)) +
  geom_col(colour = "black", position = "dodge")
subtelomere_barplot = customize_ggplot_theme(subtelomere_barplot, 
  title = NULL, xlab = "Chromosome", ylab = "Subtelomere Mean mCG Change", fill_title = "Arm") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) + 
  scale_fill_manual(breaks = c("p", "q"), values = c(p = "#d1eeea", q = "#2a5674")) +
  geom_hline(yintercept = 0) + 
  coord_flip()
subtelomere_barplot
ggsave(plot = subtelomere_barplot, "plots/subtelomere_meth_change_barplot.pdf", width = 6, height = 9, device = cairo_pdf)