# Plot 5mC across relicts in GEVE- and GEVE+ samples

# Load required packages and functions
library(methodical)
library(dplyr)
source("../auxillary_scripts/ggplot_functions.R")

# Load 5mC RSE and mask methylation of sites covered by < 10 reads or >= 350 reads
chlamy_5mc_rse = HDF5Array::loadHDF5SummarizedExperiment("CC2937_T2T_5mC_rse")
assay(chlamy_5mc_rse, 1)[assay(chlamy_5mc_rse, 2) < 10 | assay(chlamy_5mc_rse, 2) >= 350] = NA

# Get coordinates for the chr16 relict
relict_gr = rtracklayer::import.bed("../genome_files/GEVE_and_relicts_T2T.bed")[3]

# Expand relict and get 100 bp windows across the region
relict_gr_expanded = expand_granges(relict_gr, 1000, 1000)
relict_windows = unlist(tile(relict_gr_expanded, width = 100))

# Get mean mCG in relict windows and add a column with start position of each window
relict_windows_mCG = summarizeRegionMethylation(meth_rse = chlamy_5mc_rse, genomic_regions = relict_windows)
relict_windows_mCG$position = start(GRanges(row.names(relict_windows_mCG)))

# Convert chlamy_relicts_cg_meth to long format and put samples in oder
relicts_window_mCG_long = tidyr::pivot_longer(relict_windows_mCG, cols = -position, values_to = "mCG")
relicts_window_mCG_long = filter(relicts_window_mCG_long, name != "Punui2")
relicts_window_mCG_long$name = factor(relicts_window_mCG_long$name, 
  levels = c("chlamy_8A", "chlamy_15A", "chlamy_7A", "chlamy_18A"))

# Create a list for annotating relicts with genes and mAT sites
chlamy_genes_gr = readRDS("../genome_files/chlamy_genes_gr.rds")
mat_sites = readRDS("mAT_sites_whole_genome.rds")
relict_gene_class_list = readRDS("../genome_files/relict_gene_class_list.rds")
relict_genes_list_grl = lapply(relict_gene_class_list[c("BVP", "DVP")], function(x) chlamy_genes_gr[x])
relict_genes_list_grl = c(relict_genes_list_grl, list(Relict = relict_gr, mAT = expand_granges(mat_sites, 50, 50)))
relict_genes_list_grl = relict_genes_list_grl[c("mAT", "BVP", "DVP", "Relict")]

# Create a data.frame for relict genes for use with gggenes
relict_gene_class_vec = unlist(with(reshape2::melt(relict_gene_class_list), split(L1, value)))
relict_genes_df = data.frame(unlist(GRangesList(relict_genes_list_grl[c("BVP", "DVP")])))
relict_genes_df$region_type = relict_gene_class_vec[relict_genes_df$gene_id]

# Create a plot annotating the relict 
levels(relict_genes_df$region_type) = c("mAT", "BVP", "DVP", "Relict")
relict_annotation_plot = annotatePlot(plotRegionValues(relict_windows_mCG, sample_name = "chlamy_15A", geom_smooth_params = list(color = NA)), 
  annotation_grl = relict_genes_list_grl, annotation_plot_proportion = 0.25, ylab = NULL, annotation_plot_only = T,
  grl_colours = c(Relict = "#8A4F7D", DVP = NA, BVP = NA, mAT = "#C49E85")) +
  gggenes::geom_gene_arrow(data = relict_genes_df , mapping = aes(xmin = start, xmax = end, y = region_type, fill = region_type, forward = strand == "+"), 
    arrow_body_height = grid::unit(3.75, "mm"), arrowhead_height = grid::unit(3.75, "mm")) +
  labs(x = "Chromosome 16 (bp)") + scale_fill_manual(values = c("#C49E85", "#2a5674", "#d1eeea", "#8A4F7D"))
relict_annotation_plot

# Create a barplot with 5mC across 4 samples
labels = c("GEVE+ (G1)", "GEVE+ (G2)", "GEVE- (D1)", "GEVE- (D2)")
relict_cg_meth_plot = ggplot(relicts_window_mCG_long, aes(x = position, y = mCG*100)) +
  geom_col(fill = "#104C91", color = "#104C91")
relict_cg_meth_plot = customize_ggplot_theme(relict_cg_meth_plot, xlab = "Chromosome 16 (bp)", 
  ylab = "Mean mCG per 100 bp", show_legend = F, 
  facet = "name", facet_nrow = 5, facet_labels = labels, facet_scales = "fixed") +
  scale_x_continuous(labels = scales::comma, expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 50, 100), limits = c(0, 100)) +
  theme(plot.margin = margin(5, 15, 5, 5), strip.background = element_blank())
relict_cg_meth_plot_annotated = cowplot::plot_grid(
  relict_cg_meth_plot  + 
    theme(legend.position = "none", axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), axis.title.x = element_blank()), 
  relict_annotation_plot + 
    theme(axis.text = element_text(size = 14)), 
  nrow = 2, align = "v", rel_heights = c(0.825, 0.175))
relict_cg_meth_plot_annotated 
ggsave(plot = relict_cg_meth_plot_annotated, "plots/relict_5mC_4_samples.pdf", width = 9, height = 9, device = cairo_pdf)