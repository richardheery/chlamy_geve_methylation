# Make plot of 6mA across GEVE in GEVE+ and GEVE- samples

# Load required packages and functions
library(methodical)
library(dplyr)
source("../helper_scripts/ggplot_functions.R")

# Create a GRanges for whole chlamy genome and get location of GEVE
chlamy_genome_CC2937 = Biostrings::readDNAStringSet("../genome_files/CC2937_T2T.fa")[1:17]
chlamy_genome_gr = GRanges(seqinfo(chlamy_genome_CC2937))[1:17]
geve_gr = rtracklayer::import.bed("../genome_files/GEVE_T2T.bed")

# Load 6mA RSE and mask methylation of sites covered by < 10 reads or >= 350 reads
chlamy_6ma_rse = HDF5Array::loadHDF5SummarizedExperiment("CC2937_T2T_6mA_rse")
assay(chlamy_6ma_rse, 1)[assay(chlamy_6ma_rse, 2) < 10 | assay(chlamy_6ma_rse, 2) >= 350] = NA

# Expand GEVE and get 100bp windows across the region
geve_gr_expanded = expand_granges(geve_gr, 50000, 50000)
geve_windows = tile(geve_gr_expanded, width = 100)[[1]]

# Get mean mat change in windows and add column with position
geve_window_mat = summarizeRegionMethylation(meth_rse = chlamy_6ma_rse, genomic_regions = geve_windows)
geve_window_mat$position = start(GRanges(row.names(geve_window_mat)))

# Convert geve_window_mat to long format and put samples in order
geve_window_mat_long = tidyr::pivot_longer(geve_window_mat, -position, values_to = "mat")
geve_window_mat_long$name = factor(geve_window_mat_long$name, 
  levels = c("chlamy_8A", "chlamy_15A", "chlamy_7A", "chlamy_18A"))

# Create a list with GRanges for BVP, DVP, the GEVE
chlamy_genes_gr = readRDS("../genome_files/chlamy_genes_gr.rds")
gene_class_list = readRDS("../genome_files/gene_class_list.rds")
geve_genes_annotation_grl = lapply(gene_class_list[c("BVP", "DVP")], function(x) chlamy_genes_gr[x])
geve_genes_annotation_grl = c(geve_genes_annotation_grl, list(GEVE = geve_gr))

# Create a plot annotating the GEVE
geve_annotation_plot = annotatePlot(plotRegionValues(geve_window_mat, sample_name = "chlamy_15A", geom_smooth_params = list(color = NA)), 
  annotation_grl = geve_genes_annotation_grl, annotation_plot_proportion = 0.25, ylab = NULL, annotation_plot_only = T, 
  grl_colours = c(GEVE = "#8A4F7D", DVP = "#d1eeea", BVP = "#2a5674")) +
  labs(x = "Chromosome 15 (bp)")
geve_annotation_plot

# Create a barplot with methylation across 4 samples
labels = c("GEVE+ (G1)", "GEVE+ (G2)", "GEVE- (D1)", "GEVE- (D2)")
geve_at_meth_plot = ggplot(geve_window_mat_long, aes(x = position, y = mat*100)) +
  geom_col(fill = "#710627", color = "#710627")
geve_at_meth_plot = customize_ggplot_theme(geve_at_meth_plot, xlab = "Chromosome 15 (bp)", 
  ylab = "Mean mAT per 100 bp", show_legend = F, 
  facet = "name", facet_nrow = 4, facet_labels = labels, facet_scales = "fixed") +
  scale_x_continuous(labels = scales::comma, expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 50, 100)) +
  theme(plot.margin = margin(5, 15, 5, 5), strip.background = element_blank())
geve_at_meth_plot_annotated = cowplot::plot_grid(
  geve_at_meth_plot  + 
    theme(legend.position = "none", axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), axis.title.x = element_blank()), 
  geve_annotation_plot + 
    theme(axis.text = element_text(size = 14)), 
  nrow = 2, align = "v", rel_heights = c(0.85, 0.15))
geve_at_meth_plot_annotated 
ggsave(plot = geve_at_meth_plot_annotated, "plots/geve_at_methylation_4_samples.pdf", width = 9, height = 9, device = cairo_pdf)