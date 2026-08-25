# Plot coverage and 5mC across GEVE for Erazo-Garcia et al. 2025 samples

# Load required packages and functions
library(rtracklayer)
library(GenomicAlignments)
library(methodical)
source("../auxillary_scripts/ggplot_functions.R")

# Load GRanges with GEVE, expand by 50 kb at both ends and break GEVE into 1 kb windows
geve_gr = rtracklayer::import.bed("../genome_files/GEVE_T2T.bed")
geve_gr_expanded = expand_granges(geve_gr, 50000, 50000)
geve_windows = tile(geve_gr_expanded, width = 100)[[1]]

# Load GRanges with TET gene and create a GRangesList with GEVE and gene
tet_gene = rtracklayer::import.bed("../genome_files/TET_GEVEv2_173.bed")
chlamy_annotation_grl = GRangesList("TET" = tet_gene, GEVE = geve_gr)

# Get GRanges with deletion
deletion_gr = rtracklayer::import.bed("../genome_files/GEVE_deletion_T2T.bed")

# Load RSE for 5mC for original runs and liftover to T2T
Erazo_2025_5mC_rse = HDF5Array::loadHDF5SummarizedExperiment("Erazo_2025_5mC_rse")
chain = import.chain("../genome_files/Erazo_2025_to_T2T.chain", exclude = NA)
Erazo_2025_5mC_rse = liftoverMethRSE(Erazo_2025_5mC_rse, chain, )

# Mask methylation of CpGs covered by < 10 reads 
assay(Erazo_2025_5mC_rse, 1)[assay(Erazo_2025_5mC_rse, 2) < 10] = NA

# Get mean values per 100 bp
chlamy_geve_coverage = summarizeRegionMethylation(meth_rse = Erazo_2025_5mC_rse, genomic_regions = geve_windows, assay = "Cov")

# Add position as a column
chlamy_geve_coverage$position = start(GRanges(row.names(chlamy_geve_coverage)))

# Convert chlamy_geve_coverage to long format
chlamy_geve_coverage_long = tidyr::pivot_longer(chlamy_geve_coverage, -position, values_to = "coverage")

# Create a plot annotating the GEVE
geve_annotation_plot = annotatePlot(plotRegionValues(chlamy_geve_coverage, sample_name = "FirstRun", geom_smooth_params = list(color = NA)), 
  annotation_grl = chlamy_annotation_grl, annotation_plot_proportion = 0.25, ylab = NULL, annotation_plot_only = T, 
  grl_colours = c(TET = "#2a5674", GEVE = "#8A4F7D")) +
  labs(x = "Chromosome 15 (bp)") 

# Create a barplot with CpG coverage across the 3 runs and add annotation plot
geve_coverage_plot = ggplot(chlamy_geve_coverage_long, aes(x = position, y = coverage)) +
  geom_col(fill = "#66717E") + 
  geom_vline(xintercept = c(start(deletion_gr), end(deletion_gr)), linetype = "dashed", color = "#444444")
geve_coverage_plot = customize_ggplot_theme(geve_coverage_plot, xlab = "Chromosome 15 (bp)", 
  ylab = "Mean Coverage per 100 bp", show_legend = F, 
  facet = "name", facet_nrow = 3, facet_labels = c("First Run", "Second Run", "Third Run"), facet_scales = "free_y") +
  scale_x_continuous(labels = scales::comma, expand = c(0, 0)) + 
  theme(plot.margin = margin(5, 15, 5, 5), strip.background = element_blank())
geve_coverage_plot_annotated = cowplot::plot_grid(geve_coverage_plot  + theme(legend.position = "none", axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), axis.title.x = element_blank()), 
  geve_annotation_plot, nrow = 2, align = "v", rel_heights = c(0.85, 0.15))
geve_coverage_plot_annotated 
ggsave(plot = geve_coverage_plot_annotated, "plots/geve_coverage_plot_annotated.pdf", width = 9, height = 9, device = cairo_pdf)

# Get mean mCG values in windows
chlamy_geve_cg_meth = summarizeRegionMethylation(meth_rse = Erazo_2025_5mC_rse, genomic_regions = geve_windows, assay = "beta")

# Add position as a column
chlamy_geve_cg_meth$position = start(GRanges(row.names(chlamy_geve_cg_meth)))

# Convert chlamy_geve_cg_meth to long format
chlamy_geve_cg_meth_long = tidyr::pivot_longer(chlamy_geve_cg_meth, -position, values_to = "mCG")

# Create a dot plot with CpG methylation across the 3 runs and add annotation plot
geve_cg_meth_plot = ggplot(chlamy_geve_cg_meth_long, aes(x = position, y = mCG*100)) +
  geom_col(fill = "#710627", color = "#710627") +
  geom_vline(xintercept = c(start(deletion_gr), end(deletion_gr)), linetype = "dashed", color = "#444444")
geve_cg_meth_plot = customize_ggplot_theme(geve_cg_meth_plot, xlab = "Chromosome 15 (bp)", 
  ylab = "Mean mCG % per 100 bp", show_legend = F, 
  facet = "name", facet_nrow = 3, facet_labels = c("First Run", "Second Run", "Third Run"), facet_scales = "fixed") +
  scale_x_continuous(labels = scales::comma, expand = c(0, 0)) +
  theme(plot.margin = margin(5, 15, 5, 5), strip.background = element_blank())
geve_cg_meth_plot_annotated = cowplot::plot_grid(geve_cg_meth_plot  + theme(legend.position = "none", axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), axis.title.x = element_blank()), 
  geve_annotation_plot, nrow = 2, align = "v", rel_heights = c(0.85, 0.15))
geve_cg_meth_plot_annotated 
ggsave(plot = geve_cg_meth_plot_annotated, "plots/geve_cg_meth_plot_annotated.pdf", width = 9, height = 9, device = cairo_pdf)