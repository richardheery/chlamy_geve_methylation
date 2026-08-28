# Find differentially expressed TEs between GEVE- and GEVE+ samples

# Load required packages and functions
library(GenomicFeatures)
library(rtracklayer)
library(DESeq2)
library(dplyr)
source("../helper_scripts/ggplot_functions.R")

# Get coordinates of GEVE and the deletion
geve_gr = import.bed("../genome_files/GEVE_T2T.bed")
deletion_gr = rtracklayer::import.bed("../genome_files/GEVE_deletion_T2T.bed")

# Get GTF with TEs and construct names from TE locus ID (transcript ID), subfamily (gene_id), family and class
te_gtf = rtracklayer::import.gff2("../genome_files/CC2937_T2T.TEs_v3_8.bed.gtf")
names(te_gtf) = paste(te_gtf$transcript_id, te_gtf$gene_id, te_gtf$family_id, te_gtf$class_id, sep = ":")
geve_tes_gr = subsetByOverlaps(te_gtf, geve_gr, ignore.strand = T)
non_geve_tes = names(subsetByOverlaps(te_gtf, geve_gr, ignore.strand = T, invert = T))

# Filter for TEs at least 500 bp long and remove TEs overlapping the deletion
te_gtf = te_gtf[width(te_gtf) >= 500]
te_gtf = subsetByOverlaps(te_gtf, deletion_gr, invert = T)

# Get TEs on chromosome 15 (outside GEVE)
chr15_tes = names(subsetByOverlaps(te_gtf[seqnames(te_gtf) == "chromosome_15"], geve_gr, invert = T))

# Make a table with TElocal counts for all samples
telocal_files = list.files("TElocal_results", full.names = T)
names(telocal_files) = basename(tools::file_path_sans_ext(telocal_files))
telocal_counts = data.frame(lapply(telocal_files, function(x) data.frame(data.table::fread(x), row.names = 1)))
names(telocal_counts) = names(telocal_files)

# Subset telocal_counts only for TEs
telocal_counts = telocal_counts[names(te_gtf), ]

# Get list of TEs based on methylation status
te_methylation_list = readRDS("../long_read_analysis/te_methylation_list.rds")
te_methylation_vec = unlist(with(reshape2::melt(te_methylation_list), split(L1, value)))

# Create data.frame with metadata for samples for use with DESeq2
deseq_metadata = data.frame(
  sample = gsub(".*_", "", names(telocal_counts)),
  strain = gsub(".*_", "", stringr::str_sub(names(telocal_counts), end = -2)),
  day = factor(ifelse(grepl("d3", names(telocal_counts)), "D3", "D5"), levels = c("D3", "D5")),
  replicate = gsub(".*_", "", stringr::str_sub(names(telocal_counts), start = -1, end = -1)),
  row.names = names(telocal_counts)
)

# Add column for GEVE status
deseq_metadata$geve_status = factor(ifelse(deseq_metadata$strain %in% c("15", "8"), "Full", "Deletion"), levels = c("Full", "Deletion"))

# Create a DESeqDataSet with TElocal counts and fit additive models
telocal_dds = DESeq(DESeqDataSetFromMatrix(countData = telocal_counts, colData = deseq_metadata, design = ~ day + geve_status))

# Get results and indicate if changes are significant (p-value < 0.05 and absolute fold change > 1)
telocal_deletion_vs_full_results = data.frame(results(telocal_dds, contrast = c("geve_status", "Deletion", "Full")))
telocal_deletion_vs_full_results$significant = telocal_deletion_vs_full_results$padj < 0.05 & abs(telocal_deletion_vs_full_results$log2FoldChange) > 1

# Add methylation status of TEs to results
telocal_deletion_vs_full_results$methylation_status = factor(te_methylation_vec[row.names(telocal_deletion_vs_full_results)], c("mCG", "mAT", "Dual", "Unmarked"))

# Add columns with -log10 adjusted p-values and chromosome name
telocal_deletion_vs_full_results$minus_log10_padj = -log10(telocal_deletion_vs_full_results$padj)
telocal_deletion_vs_full_results$chromosome = factor(as.character(seqnames(te_gtf[row.names(telocal_deletion_vs_full_results)])), levels = seqlevels(te_gtf))

# Add name of TE family 
telocal_deletion_vs_full_results$family = ifelse(telocal_deletion_vs_full_results$significant, te_gtf$family_id, "")

# Create a volcano plot for TEs on chromosome 15 (excluding those in the GEVE)
chr15_te_volcano_plot = ggplot(telocal_deletion_vs_full_results[chr15_tes, ], aes(x = log2FoldChange, y = minus_log10_padj, fill = methylation_status, label = family)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = 1.30103, linetype = "dashed") +
  geom_point(shape = 21, size = 5, alpha = ifelse(telocal_deletion_vs_full_results[chr15_tes, ]$significant, 1, 0.25)) + 
  ggrepel::geom_text_repel(size = 3)
chr15_te_volcano_plot = customize_ggplot_theme(chr15_te_volcano_plot, 
  xlab = "Log2 Fold Change", ylab = "-Log10 Adjusted p-value",
  fill_colors = c("#C49E85", "#710627", "#011936", "#708090")) +
  theme(aspect.ratio = 1)
chr15_te_volcano_plot
ggsave(plot = chr15_te_volcano_plot, "plots/chr15_de_tes_volcano_labelled.pdf", width = 9, height = 9, device = cairo_pdf)

# Filter for significant results and add direction of change
telocal_deletion_vs_full_results_significant = filter(telocal_deletion_vs_full_results[non_geve_tes, ], significant)
telocal_deletion_vs_full_results_significant$direction = factor(ifelse(telocal_deletion_vs_full_results_significant$log2FoldChange > 0, "Upregulated", "Downregulated"))

# Summarize the number of significant TEs by chromosome, methylation status and direction
telocal_deletion_vs_full_results_significant_summary = summarise(group_by(telocal_deletion_vs_full_results_significant, chromosome, methylation_status, direction), count = n())
telocal_deletion_vs_full_results_significant_summary$direction = factor(telocal_deletion_vs_full_results_significant_summary$direction)
telocal_deletion_vs_full_results_significant_summary = tidyr::complete(ungroup(telocal_deletion_vs_full_results_significant_summary), 
  chromosome, methylation_status, direction, fill = list(count = 0))
telocal_deletion_vs_full_results_significant_summary = filter(telocal_deletion_vs_full_results_significant_summary, methylation_status %in% c("mCG", "Unmarked"))

# Create a barplot with the number of DE TEs per chromosome
de_te_barplot = ggplot(telocal_deletion_vs_full_results_significant_summary, aes(x = chromosome, y = count, fill = direction)) +
  geom_col(position = position_dodge2(width = 1, preserve = "single", padding = 0.1), color = "black") 
de_te_barplot = customize_ggplot_theme(de_te_barplot, facet = "methylation_status", x_labels = 1:17, xlab = "Chromosome",
  ylab = "Number of TEs", fill_colors = c("#4B878BFF", "#D01C1FFF"),facet_nrow = 2, facet_scales = "fixed") + 
  theme(strip.background = element_rect(fill = "grey95"), axis.text.x = element_text(size = 12))
de_te_barplot
ggsave(plot = de_te_barplot, "plots/de_te_barplot.pdf", width = 9, height = 9, device = cairo_pdf)

### Find differentially expressed TE subfamilies using TEcount output

# Make a table with class, family and subfamily for each TE and use subfamily:family:class as row names
te_metadata = dplyr::select(data.frame(te_gtf), class = class_id, family = family_id, subfamily = gene_id)
te_metadata$id = paste(te_metadata$subfamily, te_metadata$family, te_metadata$class, sep = ":")
te_metadata = filter(te_metadata, !duplicated(te_metadata$id))
te_metadata = tibble::column_to_rownames(te_metadata, "id")

# Make a table with TEcount output for all samples
tecount_files = list.files("TEcount_results", full.names = T, recursive = T)
names(tecount_files) = basename(gsub("/TEcount_out.cntTable", "", tecount_files))
tecount_counts = data.frame(lapply(tecount_files, function(x) data.frame(data.table::fread(x), row.names = 1)))
names(tecount_counts) = names(tecount_files)

# Subset tecount_counts only for TE subfamilies and put in the same order as te_metadata
tecount_counts = tecount_counts[row.names(te_metadata), ]

# Create a DESeqDataSet with TEcount counts and fit additive models
tecount_dds = DESeq(DESeqDataSetFromMatrix(countData = tecount_counts, colData = deseq_metadata, design = ~ day + geve_status))

# Get results and indicate if changes are significant (p-value < 0.05 and absolute fold change > 1)
tecount_deletion_vs_full_results = data.frame(results(tecount_dds, contrast = c("geve_status", "Deletion", "Full")))
tecount_deletion_vs_full_results$significant = tecount_deletion_vs_full_results$padj < 0.05 & abs(tecount_deletion_vs_full_results$log2FoldChange) > 1

# Add columns with -log10 adjusted p-values
tecount_deletion_vs_full_results$minus_log10_padj = -log10(tecount_deletion_vs_full_results$padj)

# Add family and subfamily as columns
tecount_deletion_vs_full_results$family = te_metadata[row.names(tecount_deletion_vs_full_results), ]$family
tecount_deletion_vs_full_results$subfamily = row.names(tecount_deletion_vs_full_results)

# Create a volcano plot for diffential TE subfamily expression
te_subfamily_volcano_plot = ggplot(tecount_deletion_vs_full_results, aes(x = log2FoldChange, y = minus_log10_padj, fill = family, label = subfamily)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = 1.30103, linetype = "dashed") +
  geom_point(shape = 21, size = 5, alpha = ifelse(tecount_deletion_vs_full_results$significant, 1, 0.25)) +
  ggrepel::geom_text_repel(data = filter(tecount_deletion_vs_full_results, significant), size = 2.5, max.overlaps = 15)
te_subfamily_volcano_plot = customize_ggplot_theme(te_subfamily_volcano_plot, 
  xlab = "Log2 Fold Change", ylab = "-Log10 Adjusted p-value", fill_title = "TE Family", title = "Differentially Expressed TE Subfamilites")
te_subfamily_volcano_plot
ggsave(plot = te_subfamily_volcano_plot, "plots/te_subfamily_volcano_plot.pdf", width = 16, height = 9, device = cairo_pdf)

### Plot subfamilies by mean expression in GEVE+ and GEVE-

# Get normalized counts for TE subfamilies
tecount_counts = data.frame(counts(te_count_dds, normalized = T))

# Separate values for GEVE- and GEVE+ samples
tecount_counts_deletion = dplyr::select(tecount_counts, matches("_7|_18"))
tecount_counts_wt = dplyr::select(tecount_counts, matches("_8|_15"))

# Get the mean counts for each subfamily in GEVE- and GEVE+ samples
tecount_counts_summary = data.frame(
  deletion = rowMeans(tecount_counts_deletion),
  wt = rowMeans(tecount_counts_wt)
)
tecount_counts_summary = tibble::rownames_to_column(tecount_counts_summary, "te_subfamily")

# Convert to long format 
tecount_counts_summary_long = tidyr::pivot_longer(tecount_counts_summary, cols = -te_subfamily, names_to = "status", values_to = "counts")

# Rank subfamilies from most to least expressed
tecount_counts_summary_long = arrange(mutate(group_by(tecount_counts_summary_long, status), rank = rank(-counts)), rank)
tecount_counts_summary_long = filter(tecount_counts_summary_long, rank <= 20)

# Update GEVE status
tecount_counts_summary_long$status = factor(ifelse(tecount_counts_summary_long$status == "deletion", "GEVE-", "GEVE+"), levels = c("GEVE+", "GEVE-"))

# Group results by status and subfamily and make a vector linking the groups to their associated rank
tecount_counts_summary_long$grouping = paste(tecount_counts_summary_long$status, tecount_counts_summary_long$te_subfamily, sep = ":")
grouping_to_rank = setNames(tecount_counts_summary_long$rank, tecount_counts_summary_long$grouping)

# Convert tecount_cpm to long format and add status and day for samples
tecount_counts_long = tidyr::pivot_longer(tibble::rownames_to_column(tecount_counts, "te_subfamily"), 
  cols = -te_subfamily, names_to = "sample", values_to = "counts")
tecount_counts_long$status = factor(ifelse(grepl("_7|_18", tecount_counts_long$sample), "GEVE-", "GEVE+"), c("GEVE+", "GEVE-"))
tecount_counts_long$day = ifelse(grepl("d3", tecount_counts_long$sample), "D3", "D5")

# Filter for groups in tecount_counts_summary_long and add associated rank
tecount_counts_long$grouping = paste(tecount_counts_long$status, tecount_counts_long$te_subfamily, sep = ":")
tecount_counts_long = filter(tecount_counts_long,  grouping %in% tecount_counts_summary_long$grouping)
tecount_counts_long$rank = grouping_to_rank[tecount_counts_long$grouping]

# Create barplot of TE subfamilies by expression in GEVE- and GEVE+ samples
te_rank_plot = ggplot(tecount_counts_summary_long, aes(x = tidytext::reorder_within(te_subfamily, -rank, status), y = counts, fill = status)) +
  geom_col(colour = "black", alpha = 0.75) + 
  geom_point(data = tecount_counts_long, aes(color = day), size = 3) + guides(fill = "none")
te_rank_plot = customize_ggplot_theme(te_rank_plot, 
  xlab = "TE Subfamily", x_labels_angle = 45, ylab = "Mean Normalized Counts", 
  fill_colors = c("GEVE-" = "#b97978", "GEVE+" = "#e1d7AA"), colors = c("D3" = "#BCBDDC", "D5" = "#54278F"),
  show_legend = T,
  facet = "status", facet_nrow = 2, facet_scales = "free_x") +
  theme(strip.background = element_rect(fill = "grey95")) +
  tidytext::scale_x_reordered() +
  theme(axis.text.x = element_text(size = 12)) 
te_rank_plot 
ggsave(plot = te_rank_plot, "plots/te_rank_plot.pdf", width = 9, height = 9, device = cairo_pdf)