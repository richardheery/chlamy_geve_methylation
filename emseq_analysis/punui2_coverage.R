### Create coverage plot for Punui2

# Load required packages and functions
library(methodical)
library(SummarizedExperiment)
source("../auxillary_scripts/ggplot_functions.R")

# Read in CC2937 T2T as a DNAStringSet
chlamy_genome_CC2937 = Biostrings::readDNAStringSet("../genome_files/CC2937_T2T.fa")[1:17]

# Get meth sites for CC2937 genome
CC2937_cpgs_gr = extractMethSitesFromGenome(chlamy_genome_CC2937, standard_seqs_only = F)

# Get path to CGmap file for viral Punui2 sample
punui2_cgmap_file = "cgmap_files/Punui2.CGmap.gz"

# Create colData for Punui2
coldata_viral = data.frame(
  status = "Viral",
  row.names = "Punui2"
)

# Create an RSE for punui2_cgmap_file
punui2_5mc_rse = makeMethRSEFromInputFiles(meth_files = punui2_cgmap_file, 
  seqnames_col = 1, start_col = 3, total_reads_col = 8, meth_reads_col = 7, 
  zero_based = F, collapse_strands = T, sequence_context = "CG", meth_sites = CC2937_cpgs_gr, 
  sample_metadata = coldata_viral, hdf5_dir = "punui2_CC2937_T2T_5mC_rse")

# Mask methylation of CpGs covered by < 10 reads or >= 350 reads
assay(punui2_5mc_rse, 1)[assay(punui2_5mc_rse, 2) < 10 | assay(punui2_5mc_rse, 2) >= 350] = NA

# Load GRanges with GEVE, expand by 50 kb at both ends and break GEVE into 1 kb windows
geve_gr = rtracklayer::import.bed("../genome_files/GEVE_T2T.bed")
chlamy_annotation_grl = GRangesList(GEVE = geve_gr)
geve_gr_expanded = expand_granges(geve_gr, 200000, 200000)
geve_windows = tile(geve_gr_expanded, width = 100)[[1]]

# Calculate the mean CpG coverage and methylation within the GEVE
punui2_5mc_rse_geve = subsetByOverlaps(punui2_5mc_rse, geve_gr, ignore.strand = T)
colMeans(assay(punui2_5mc_rse_geve, "Cov"), na.rm = T)
(sum(assay(punui2_5mc_rse_geve, "Cov") * assay(punui2_5mc_rse_geve, "beta"), na.rm = T))/
  (sum(assay(punui2_5mc_rse_geve, "Cov"), na.rm = T))

# Get mean values per 100 bp
punui_geve_coverage = summarizeRegionMethylation(meth_rse = punui2_5mc_rse, genomic_regions = geve_windows, assay = "Cov")

# Add position as a column
punui_geve_coverage$position = start(GRanges(row.names(punui_geve_coverage)))

# Convert punui_geve_coverage to long format
punui_geve_coverage_long = tidyr::pivot_longer(punui_geve_coverage, -position, values_to = "coverage")

# Create a plot annotating the GEVE
geve_annotation_plot = annotatePlot(plotRegionValues(punui_geve_coverage, sample_name = "Punui2", geom_smooth_params = list(color = NA)), 
  annotation_grl = chlamy_annotation_grl, annotation_plot_proportion = 0.25, ylab = NULL, annotation_plot_only = T, 
  grl_colours = c("#8A4F7D")) +
  labs(x = "Chromosome 15 (bp)")

# Create a barplot with CpG coverage across the 3 runs and add annotation plot
geve_coverage_plot = ggplot(punui_geve_coverage_long, aes(x = position, y = coverage)) +
  geom_col(fill = "#66717E") 
geve_coverage_plot = customize_ggplot_theme(geve_coverage_plot, xlab = "Chromosome 15 (bp)", 
  ylab = "Mean Coverage\nper 100 bp", show_legend = F) +
  scale_x_continuous(labels = scales::comma, expand = c(0, 0)) + 
  theme(plot.margin = margin(5, 15, 5, 5), strip.background = element_blank())
geve_annotation_plot = geve_annotation_plot + coord_cartesian(xlim = ggplot_build(geve_coverage_plot)$layout$panel_params[[1]]$x.range)
geve_coverage_plot_annotated = cowplot::plot_grid(geve_coverage_plot  + theme(legend.position = "none", axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), axis.title.x = element_blank()), 
  geve_annotation_plot, nrow = 2, align = "v", rel_heights = c(0.8, 0.25))
geve_coverage_plot_annotated = geve_coverage_plot_annotated + coord_cartesian(clip = "on")
ggsave(plot = geve_coverage_plot_annotated, "plots/punui_coverage_plot_annotated.pdf", width = 16, height = 4.5, device = cairo_pdf)