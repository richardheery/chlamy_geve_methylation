### Make plots showing CG and AT methylation in GEVE

# Load required packages and functions
library(methodical)
library(HDF5Array)
library(dplyr)
source("../auxillary_scripts/ggplot_functions.R")

# Load GRanges for annotating chlamy genome
geve_gr = import.bed("../genome_files/GEVE_T2T.bed")
chlamy_genes_gr = readRDS("../genome_files/chlamy_genes_gr.rds")

# Load RSEs for 5mC and 6mA and mask methylation of sites covered by < 10 reads or >= 350 reads
chlamy_cg_meth_rse = loadHDF5SummarizedExperiment("CC2937_T2T_5mC_rse")
chlamy_at_meth_rse = loadHDF5SummarizedExperiment("CC2937_T2T_6mA_rse")
assay(chlamy_cg_meth_rse, 1)[assay(chlamy_cg_meth_rse, 2) < 10 | assay(chlamy_cg_meth_rse, 2) >= 350] = NA
assay(chlamy_at_meth_rse, 1)[assay(chlamy_at_meth_rse, 2) < 10 | assay(chlamy_at_meth_rse, 2) >= 350] = NA

### Make plot showing methylation of BVP and DVP genes in GEVE

# Create a GRangesList with BVP and DVP genes and the GEVE
gene_class_list = readRDS("../genome_files/gene_class_list.rds")
geve_genes_annotation_grl = lapply(gene_class_list[c("BVP", "DVP")], function(x) chlamy_genes_gr[x])
geve_genes_annotation_grl = c(geve_genes_annotation_grl, list(GEVE = geve_gr))

# Expand GEVE region and get 1kb windows across the region
geve_gr_expanded = expand_granges(geve_gr, 50000, 50000)
geve_windows = tile(geve_gr_expanded, width = 1000)[[1]]

# Get mean mCG and mAT in 1kb windows in GEVE in chlamy_15A
geve_window_mCG = summarizeRegionMethylation(meth_rse = chlamy_cg_meth_rse[, "chlamy_15A"], genomic_regions = geve_windows)
geve_window_mCG$position = start(GRanges(row.names(geve_window_mCG)))
geve_window_mAT = summarizeRegionMethylation(meth_rse = chlamy_at_meth_rse[, "chlamy_15A"], genomic_regions = geve_windows)
geve_window_mAT$position = start(GRanges(row.names(geve_window_mAT)))

# Combine geve_window_mAT and geve_window_mCG
geve_windows_combined = bind_rows(
  list(CG = geve_window_mCG, AT = geve_window_mAT),
  .id = "site"
)

# Create a plot annotating the GEVE
geve_annotation_plot = annotatePlot(plotRegionValues(geve_window_mCG, sample_name = "chlamy_15A", geom_smooth_params = list(color = NA)), 
  annotation_grl = geve_genes_annotation_grl, annotation_plot_proportion = 0.25, ylab = NULL, annotation_plot_only = T, 
  grl_colours = c("#8A4F7D", "#d1eeea", "#2a5674")) +
  labs(x = "Chromosome 15 (bp)") 
geve_annotation_plot

# Create a scatter plot with ApT and CpG methylation in the GEVE in chlamy_15A 
geve_gene_meth_plot = ggplot(geve_windows_combined, aes(x = position, y = chlamy_15A*100, fill = site)) +
  geom_point(shape = 21, size = 3, alpha = 1)
geve_gene_meth_plot = customize_ggplot_theme(geve_gene_meth_plot, xlab = "Chromosome 15 (bp)", 
  ylab = "Methylation Percentage per 1 kb", show_legend = T, fill_colors = c("#C49E85", "#710627")) +
  theme(plot.margin = margin(5, 15, 5, 5), strip.background = element_blank())
geve_gene_meth_plot_annotated = cowplot::plot_grid(geve_gene_meth_plot  + theme(axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), axis.title.x = element_blank()), 
  geve_annotation_plot, nrow = 2, align = "v", rel_heights = c(0.825, 0.175)) 
geve_gene_meth_plot_annotated 
ggsave(plot = geve_gene_meth_plot_annotated, "plots/geve_gene_mCG_vs_mAT_plot.pdf", width = 16, height = 9, device = cairo_pdf)

### Zoom in on region between chromosome_15:1150000-1170000 and show CG and AT methylation for individual sites

# Define GRanges for regions to magnify
zoom_region_gr = GRanges("chromosome_15:1150000-1170000")

# Get site-level mCG and mAT values in zoom in chlamy_8A and chlamy_18A
zoom_region_mCG = extractGRangesMethSiteValues(meth_rse = chlamy_cg_meth_rse[, c("chlamy_8A", "chlamy_18A")], genomic_regions = zoom_region_gr)
zoom_region_mAT = extractGRangesMethSiteValues(meth_rse = chlamy_at_meth_rse[, c("chlamy_8A", "chlamy_18A")], genomic_regions = zoom_region_gr)

# Combine mCG and mAT values for zoom_region_gr
zoom_region_gr_combined = bind_rows(
  list(CG = zoom_region_mCG, AT = zoom_region_mAT),
  .id = "site"
)

# Add a column with position of site
zoom_region_gr_combined$position = start(GRanges(row.names(zoom_region_gr_combined)))

# Convert zoom_region_gr_combined to long format 
zoom_region_gr_combined_long = data.frame(tidyr::pivot_longer(zoom_region_gr_combined, -c(site, position), 
  names_to = "sample", values_to = "value"))
zoom_region_gr_combined_long = filter(zoom_region_gr_combined_long, value > 0 | site == "CG")
zoom_region_gr_combined_long$sample = factor(zoom_region_gr_combined_long$sample,  levels = c("chlamy_8A", "chlamy_18A"))

# Create a GRangesList with the GEVE and all chlamy genes and the GEVE and a data.frame with all genes on chr15 for use with gggenes
zoom_annotation_grl = list(Genes = chlamy_genes_gr, GEVE = geve_gr)
geve_genes_df = data.frame(chlamy_genes_gr)
geve_genes_df = filter(geve_genes_df, seqnames == "chromosome_15")
geve_genes_df$region_type = "Genes"

# Create annotation plot for zoom region
zoom_annotation_plot = annotatePlot(plotRegionValues(zoom_region_gr_combined, sample_name = "chlamy_8A", geom_smooth_params = list(color = NA)), 
  annotation_grl = zoom_annotation_grl, annotation_plot_proportion = 0.25, ylab = NULL, annotation_plot_only = T, 
  grl_colours = c(GEVE = "#8A4F7D", Genes = NA)) +
  labs(x = "Chromosome 15 (bp)") +
  gggenes::geom_gene_arrow(data = geve_genes_df, mapping = aes(xmin = start, xmax = end, y = region_type, fill = region_type, forward = strand == "+"), 
    arrow_body_height = grid::unit(3.75, "mm"), arrowhead_height = grid::unit(3.75, "mm")) +
  scale_fill_manual(values = c(Genes = "#2a5674"))
zoom_annotation_plot

# Create a scatter plot with ApT and CpG methylation in zoom region in chlamy_8A and chlamy_18A
zoom_region_meth_plot = ggplot(zoom_region_gr_combined_long, aes(x = position, y = value*100, fill = site)) +
  geom_point(shape = 21, size = 3, alpha = 1) 
zoom_region_meth_plot = customize_ggplot_theme(zoom_region_meth_plot, xlab = "Chromosome 15 (bp)", 
  ylab = "Methylation % per Site", show_legend = T, fill_colors = c("#C49E85", "#710627"),
  facet = "sample", facet_nrow = 2, facet_labels = c("GEVE+ (8A)", "GEVE- (18A)"), facet_scales = "fixed") +
  scale_x_continuous(labels = scales::comma, expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 50, 100)) +
  theme(plot.margin = margin(5, 15, 5, 5), strip.background = element_blank())
zoom_region_meth_plot_annotated = cowplot::plot_grid(
  zoom_region_meth_plot  + 
    theme(axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), axis.title.x = element_blank()), 
  zoom_annotation_plot + 
    theme(axis.text = element_text(size = 14)),
  nrow = 2, align = "v", rel_heights = c(1.35, 0.35))
zoom_region_meth_plot_annotated
ggsave(plot = zoom_region_meth_plot_annotated, "plots/zoom_region_meth_plot_5mC_vs_6mA.pdf", width = 16, height = 6.8, device = cairo_pdf)
