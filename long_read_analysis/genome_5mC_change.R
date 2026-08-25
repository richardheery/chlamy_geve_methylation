# Find windows gaining and losing 5mC change across genome in GEVE- samples and make plots showing 5mC change 

# Load required packages and functions
library(methodical)
library(dplyr)
library(RIdeogram)
library(karyoploteR)
source("../helper_scripts/ggplot_functions.R")

# Create a GRanges for whole chlamy genome and get location of GEVE
chlamy_genome_CC2937 = Biostrings::readDNAStringSet("../genome_files/CC2937_T2T.fa")[1:17]
chlamy_genome_gr = GRanges(seqinfo(chlamy_genome_CC2937))[1:17]
geve_gr = rtracklayer::import.bed("../genome_files/GEVE_T2T.bed")

# Load 5mC RSE and mask methylation of sites covered by < 10 reads or >= 350 reads
chlamy_5mc_rse = HDF5Array::loadHDF5SummarizedExperiment("CC2937_T2T_5mC_rse")
assay(chlamy_5mc_rse, 1)[assay(chlamy_5mc_rse, 2) < 10 | assay(chlamy_5mc_rse, 2) >= 350] = NA

# Create a GRanges with 1 kb windows across whole genome
genome_windows = unlist(tile(chlamy_genome_gr, width = 1000))
names(genome_windows) = NULL

# Extract CG methylation values for genome windows
genome_window_5mC = summarizeRegionMethylation(meth_rse = chlamy_5mc_rse, genomic_regions = genome_windows, 
  BPPARAM = BiocParallel::SerialParam())

# Add difference in 5mCG between GEVE− and GEVE+ samples
genome_window_5mC$deletion_change = 
  rowMeans(genome_window_5mC[c("chlamy_7A", "chlamy_18A")], na.rm = T) - 
  rowMeans(genome_window_5mC[c("chlamy_15A", "chlamy_8A")], na.rm = T)

# Define hypermethylated windows as those with a methylation change > 0.2, hypomethylated windows with a change < -0.2
hypermethylated_windows = GRanges(row.names(genome_window_5mC)[which(genome_window_5mC$deletion_change > 0.15)])
hypomethylated_windows = GRanges(row.names(genome_window_5mC)[which(genome_window_5mC$deletion_change < -0.15)])
unchanged_windows = GRanges(row.names(genome_window_5mC)[which(abs(genome_window_5mC$deletion_change) < 0.05)])
rtracklayer::export.bed(hypermethylated_windows, "deletion_hypermethylated_windows.bed")
rtracklayer::export.bed(hypomethylated_windows, "deletion_hypomethylated_windows.bed")
rtracklayer::export.bed(unchanged_windows, "deletion_unchanged_windows.bed")

### Make karyotype plot of 5mC change across whole genome

# Create a GRanges with 5mC change across genome
genome_window_5mC_gr = GRanges(row.names(genome_window_5mC), value = genome_window_5mC$deletion_change)

# Load GRanges with centromeres and GEVE
centromeres_gr = rtracklayer::import.bed("../genome_files/CC2937_T2T.centromeres.bed")
geve_gr = rtracklayer::import.bed("../genome_files/GEVE_T2T.bed")

# Create a GRanges to use as cytobands to display the centromeres and GEVE
centromeres_gr$gieStain = "centromere"
non_centromeres_gr = setdiff(chlamy_genome_gr, centromeres_gr)
non_centromeres_gr$gieStain = "none"
geve_gr$gieStain = "GEVE"
cytobands_gr = c(centromeres_gr, non_centromeres_gr, geve_gr)

# Update seqlevels of GRanges to shorten sequence names
seqlevels(chlamy_genome_gr) = paste0("chr", 1:17)
seqlevels(cytobands_gr) = paste0("chr", 1:17)
seqlevels(genome_window_5mC_gr) = paste0("chr", 1:17)

# Make karyoplot showing 5mC change across genome
cairo_pdf("plots/karyoplot_5mC_change.pdf", 16, 9)
kp = plotKaryotype(genome = chlamy_genome_gr, cytobands = cytobands_gr)
kpAddBaseNumbers(kp, tick.dist = 1000000, units = "Mb", add.units = T)
kpAxis(kp, r0 = 0, r1 = 1, ymin = -0.7, ymax = 0.7, tick.pos = 0)
kpLines(kp, genome_window_5mC_gr, col = "#CD2626", r0 = 0, r1 = 1, ymin = -0.7, ymax = 0.7)
kpAddCytobands(kp, color.table = c(centromere = "#D2C465", GEVE = "#8A4F7D", none = "grey95"))
dev.off()

### Make ideogram plot of 5mC change across whole genome

# Make a data.frame for the karyotype
karyotype_df = as.data.frame(chlamy_genome_gr, row.names = NULL)[1:3]
names(karyotype_df) = c("Chr", "Start", "End")
karyotype_df$Chr = as.character(karyotype_df$Chr)

# Merge adjacent regions with same score
genome_window_5mC_gr$rounded_value = round(genome_window_5mC_gr$value, digits = 2)
genome_window_5mC_merged_gr = sort(unlist(GRangesList(lapply(
  split(genome_window_5mC_gr, genome_window_5mC_gr$rounded_value), function(x) reduce(x)))))
genome_window_5mC_merged_gr$Values = as.numeric(names(genome_window_5mC_merged_gr))

# Convet genome_window_5mC_merged_gr to a data.frame
genome_values_df = transmute(data.frame(genome_window_5mC_merged_gr), 
  Chr = seqnames, Start = start, End = end, Values = Values)

# Split geve_gr into 100 kb regions and make into a data.frame including display parameters
geve_gr_tiled = unlist(tile(geve_gr, width = 100000))
geve_df = transmute(data.frame(geve_gr_tiled), 
  Type = "GEVE", Shape = "box", Chr = "chr15", Start = start, End = end, color = "8A4F7D")

# Make a data.frame from centromeres including display parameters and combine with GEVE to make annotation data.frame for ideogram
centromeres_df = transmute(data.frame(centromeres_gr), 
  Type = "Centromere", Shape = "box", Chr = paste0("chr", gsub("chromosome_", "", as.numeric(seqnames))), 
  Start = start, End = end, color = "D2C465")
annotation_df = rbind(geve_df, centromeres_df)

# Make ideogram plot showing 5mC change across genome 
ideogram(karyotype = karyotype_df, overlaid = genome_values_df, label = annotation_df, label_type = "marker",
  colorset1 = c(low = "#53868B", mid = "white", high = "#CD2626"), output = "plots/deletion_5mc_change_ideogram.svg")

### Make line plots of 5mC change for individual chromosomes

# Import coordinates for satellites and filter for those over 50 kb
satellites_gr = rtracklayer::import.bed("../genome_files/CC2937_T2T.satellites.bed")
satellites_gr = satellites_gr[width(satellites_gr) > 50000]

# Split genome_window_5mC by chromosome
genome_window_5mC$position = start(GRanges(row.names(genome_window_5mC)))
genome_window_5mC_list = split(genome_window_5mC, seqnames(GRanges(row.names(genome_window_5mC))))
names(genome_window_5mC_list) = paste("Chromosome", seq_along(genome_window_5mC_list))

# Make a GRangesList with centromeres, satellites and the GEVE
chlamy_annotation_grl = GRangesList(
  Centromere = centromeres_gr,
  GEVE = geve_gr,
  Satellites = satellites_gr
)

# Make line plots of 5mC change for each chromosome 
chromosome_plots = lapply(names(genome_window_5mC_list), function(x) {
  
  # Create a line plot of 5mC change
  chr_change_plot = ggplot(genome_window_5mC_list[[x]], aes(x = position, y = deletion_change*100)) +
    geom_line(color = "#710627") +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed")
  chr_change_plot = customize_ggplot_theme(chr_change_plot, xlab = paste(x, "(bp)"), 
    ylab = "mCG Change in 1 kb Windows", show_legend = F) +
    scale_x_continuous(labels = scales::comma, expand = expansion(add = c(100000, 100000))) +
    scale_y_continuous(limits = c(-70, 70), expand = expansion(0, 0)) +
    theme(plot.margin = margin(5, 15, 5, 5), strip.background = element_blank()) 
  
  # Annotate the plot with chlamy_annotation_grl
  chr_change_plot_annotated = annotatePlot(chr_change_plot, annotation_grl = chlamy_annotation_grl, ylab = NULL, annotation_plot_proportion = 0.2,
    grl_colours = c(Centromere = "#D2C465", GEVE = "#8A4F7D", Satellites = "#DBCBD8"))
  chr_change_plot_annotated
  
})

# Save plots for chromosomes 5, 9 13 and 15
ggsave(plot = chromosome_plots[[5]], "plots/deletion_5mc_change_chr5_lineplot.pdf.pdf", width = 16, height = 9)
ggsave(plot = chromosome_plots[[9]], "plots/deletion_5mc_change_chr9_lineplot.pdf.pdf", width = 16, height = 9)
ggsave(plot = chromosome_plots[[13]], "plots/deletion_5mc_change_chr13_lineplot.pdf.pdf", width = 16, height = 9)
ggsave(plot = chromosome_plots[[15]], "plots/deletion_5mc_change_chr15_lineplot.pdf", width = 16, height = 9)

### Make plot of 5mC in GEVE+ and GEVE- samples along with Punui virus

# Load meth RSE for Punui2 and add to chlamy_5mc_rse
punui2_5mc_rse = HDF5Array::loadHDF5SummarizedExperiment("../emseq_analysis/punui2_CC2937_T2T_5mC_rse/")
assay(punui2_5mc_rse, 1)[assay(punui2_5mc_rse, 2) < 10 | assay(punui2_5mc_rse, 2) >= 350] = NA
chlamy_5mc_rse = cbind(chlamy_5mc_rse, punui2_5mc_rse)

# Expand GEVE and get 100bp windows across the region
geve_gr_expanded = expand_granges(geve_gr, 50000, 50000)
geve_windows = tile(geve_gr_expanded, width = 100)[[1]]

# Get mean mCG in GEVE windows and add column with position
geve_window_mCG = summarizeRegionMethylation(meth_rse = chlamy_5mc_rse, genomic_regions = geve_windows)
geve_window_mCG$position = start(GRanges(row.names(geve_window_mCG)))

# Convert geve_window_mCG_long to long format and put samples in order
geve_window_mCG_long = tidyr::pivot_longer(geve_window_mCG, -position, values_to = "mCG")
geve_window_mCG_long$name = factor(geve_window_mCG_long$name, 
  levels = c("chlamy_8A", "chlamy_15A", "chlamy_7A", "chlamy_18A", "Punui2"))

# Create a list with GRanges for BVP, DVP, the GEVE
chlamy_genes_gr = readRDS("../genome_files/chlamy_genes_gr.rds")
gene_class_list = readRDS("../genome_files/gene_class_list.rds")
geve_genes_annotation_grl = lapply(gene_class_list[c("BVP", "DVP")], function(x) chlamy_genes_gr[x])
geve_genes_annotation_grl = c(geve_genes_annotation_grl, list(GEVE = geve_gr))

# Create a plot annotating the GEVE
geve_annotation_plot = annotatePlot(plotRegionValues(geve_window_mCG, sample_name = "chlamy_15A", geom_smooth_params = list(color = NA)), 
  annotation_grl = geve_genes_annotation_grl, annotation_plot_proportion = 0.25, ylab = NULL, annotation_plot_only = T, 
  grl_colours = c(GEVE = "#8A4F7D", DVP = "#d1eeea", BVP = "#2a5674")) +
  labs(x = "Chromosome 15 (bp)")
geve_annotation_plot

# Create a barplot with methylation across 4 samples and Punuivirus
labels = c("GEVE+ (G1)", "GEVE+ (G2)", "GEVE- (D1)", "GEVE- (D2)", "Punuivirus Virion DNA")
geve_cg_meth_plot = ggplot(geve_window_mCG_long, aes(x = position, y = mCG*100)) +
  geom_col(fill = "#710627", color = "#710627")
geve_cg_meth_plot = customize_ggplot_theme(geve_cg_meth_plot, xlab = "Chromosome 15 (bp)", 
  ylab = "Mean mCG per 100 bp", show_legend = F, 
  facet = "name", facet_nrow = 5, facet_labels = labels, facet_scales = "fixed") +
  scale_x_continuous(labels = scales::comma, expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 50, 100)) +
  theme(plot.margin = margin(5, 15, 5, 5), strip.background = element_blank())
geve_cg_meth_plot_annotated = cowplot::plot_grid(
  geve_cg_meth_plot  + 
    theme(legend.position = "none", axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), axis.title.x = element_blank()), 
  geve_annotation_plot + 
    theme(axis.text = element_text(size = 14)), 
  nrow = 2, align = "v", rel_heights = c(0.85, 0.15))
geve_cg_meth_plot_annotated 
ggsave(plot = geve_cg_meth_plot_annotated , "plots/geve_methylation_4_samples_with_virus.pdf", width = 9, height = 9, device = cairo_pdf)

### Make boxplots of mCG change in 100 bp GEVE windows with and without mAT

# Import location of mAT sites 
mat_sites = readRDS("mAT_sites_whole_genome.rds")

# Create a data.frame with 5mC change for each GEVE window and indicate if they overlap mAT sites 
geve_window_mCG_change_vs_mAT_df = data.frame(
  deletion_change = 
    rowMeans(geve_window_mCG[c("chlamy_7A", "chlamy_18A")], na.rm = T) -
    rowMeans(geve_window_mCG[c("chlamy_15A", "chlamy_8A")], na.rm = T),
  mAT = ifelse(GRanges(row.names(geve_window_mCG)) %over% mat_sites, "mAT+", "mAT-")
)

# Create boxplots for mCG change in windows based on the presence of mAT sites
set.seed(123)
geve_window_5mC_change_boxplots = ggplot(geve_window_mCG_change_vs_mAT_df, aes(x = mAT, y = deletion_change, fill = mAT)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.25) +
  geom_jitter(shape = 21)
geve_window_5mC_change_boxplots = customize_ggplot_theme(geve_window_5mC_change_boxplots, 
  ylab = "5mC Change in \n100 bp Windows", show_legend = F, fill_colors = c("#4C5760", "#C49E85")) +
  theme(aspect.ratio = 1)
geve_window_5mC_change_boxplots
ggsave(plot = geve_window_5mC_change_boxplots, filename = "plots/geve_window_5mC_change_vs_6ma_boxplots.pdf", width = 9, height = 9)