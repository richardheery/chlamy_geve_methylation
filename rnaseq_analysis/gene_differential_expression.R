# Find differentially expressed genes between GEVE- and GEVE+ samples

# Load required packages
library(dplyr)
library(GenomicFeatures)
library(rtracklayer)
library(DESeq2)
source("../helper_scripts/ggplot_functions.R")
source("../helper_scripts/enrichment_tests.R")

# Get locations of genes 
chlamy_genes_gr = readRDS("../genome_files/chlamy_genes_gr.rds")
chlamy_genes = names(chlamy_genes_gr)

# Get names of GEVE genes and genes in the deletion
geve_genes_gr = rtracklayer::import.bed("../genome_files/GEVE_all_genes.bed")
geve_genes = geve_genes_gr$name
deletion_gr = rtracklayer::import.bed("../genome_files/GEVE_deletion_T2T.bed")
deletion_genes = subsetByOverlaps(geve_genes_gr, deletion_gr)$name
non_geve_genes = setdiff(chlamy_genes, geve_genes)

# Subset geve_genes for those not in the deletion
geve_genes = setdiff(geve_genes, deletion_genes)

# Get genes on chromosome 15 (outside GEVE)
chr15_genes = setdiff(chlamy_genes_gr[seqnames(chlamy_genes_gr) == "chromosome_15"]$gene_id, c(geve_genes, deletion_genes))

# Get paths to all TEcount files
tecount_files = list.files("TEcount_results", full.names = T, recursive = T)

# Read in files as a table
tecount_table = data.frame(lapply(tecount_files, function(x) data.table::fread(x)[[2]]), 
  row.names = data.table::fread(tecount_files[1])[[1]])
names(tecount_table) = basename(dirname(tecount_files))

# Subset tecount_table for genes
all_gene_counts = tecount_table[chlamy_genes, ]

# Create data.frame with metadata for samples for use with DESeq2
deseq_metadata = data.frame(
  sample = gsub(".*_", "", names(tecount_table)),
  strain = gsub(".*_", "", stringr::str_sub(names(tecount_table), end = -2)),
  day = factor(ifelse(grepl("d3", names(tecount_table)), "D3", "D5"), levels = c("D3", "D5")),
  replicate = gsub(".*_", "", stringr::str_sub(names(tecount_table), start = -1, end = -1)),
  row.names = names(tecount_table)
)

# Add column for GEVE status
deseq_metadata$geve_status = ifelse(deseq_metadata$strain %in% c("15", "8"), "Full", "Deletion")

# Create a DESeqDataSet for all genes and get VST normalized counts
all_gene_counts_dds = DESeq(DESeqDataSetFromMatrix(countData = all_gene_counts, colData = deseq_metadata, design = ~ day + geve_status))
all_gene_counts_vst_normalized = assay(vst(all_gene_counts_dds, blind = FALSE))

### Make PCA plot of samples 

# Remove deletion genes and filter for genes expressed in at least 1 sample
non_deletion_gene_counts_normalized = all_gene_counts_vst_normalized[setdiff(chlamy_genes, deletion_genes), ]
non_deletion_gene_counts_normalized = non_deletion_gene_counts_normalized[rowSums(non_deletion_gene_counts_normalized) > 0, ]

# Perform PCA
gene_counts_pca = prcomp(t(non_deletion_gene_counts_normalized), center = T, scale. = F)

# Calculate percentage of variance explained by each PC
pca_var_per <- round(gene_counts_pca$sdev ^ 2 / sum(gene_counts_pca$sdev^2) * 100, 1)

# Create a data.frame with PC1 and PC2 values added to deseq_metadata
pca_df = mutate(deseq_metadata, PC1 = gene_counts_pca$x[, 1], PC2 = gene_counts_pca$x[, 2])

# Create PCA plot
gene_counts_pca_plot = ggplot(pca_df, aes(x = PC1, y = PC2, colour = geve_status, shape = day, label = sample)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 5, colour = "black")
gene_counts_pca_plot = customize_ggplot_theme(gene_counts_pca_plot, 
  colors =  c("#b97978", "#e1d7AA"), 
  xlab = paste("PC1", ": ", pca_var_per[1], "%", sep=""), 
  ylab = paste("PC2", ": ", pca_var_per[2], "%", sep="")
) + labs(shape = NULL) + 
  theme(legend.position = c(0.2, 0.1), legend.background = element_rect(colour = "black", fill = "white", linewidth = 0.5), legend.box = "horizontal")
gene_counts_pca_plot
ggsave(plot = gene_counts_pca_plot, "plots/gene_counts_pca_plot.pdf",
  width = 9, height = 9, device = cairo_pdf)

###  Plot expression rank of GEVEv2_173 (TET) in GEVE+ samples

# Get normalized counts for all GEVE genes
all_geve_gene_counts_normalized = data.frame(counts(all_gene_counts_dds, normalized = T))[c(geve_genes, deletion_genes), ]

# Get names of GEVE+ samples and subset for these samples
full_samples = row.names(filter(deseq_metadata, geve_status == "Full"))
all_geve_gene_counts_normalized_full_samples = all_geve_gene_counts_normalized[, full_samples]

# Create a data.frame with mean gene expression of genes across GEVE+ samples
all_geve_gene_counts_normalized_full_samples_mean = data.frame(gene_id = row.names(all_geve_gene_counts_normalized_full_samples), 
  mean_expression = log10(rowMeans(all_geve_gene_counts_normalized_full_samples)))
all_geve_gene_counts_normalized_full_samples_mean = filter(all_geve_gene_counts_normalized_full_samples_mean, is.finite(mean_expression))

# Rank genes in order of descending mean expression
all_geve_gene_counts_normalized_full_samples_mean = arrange(all_geve_gene_counts_normalized_full_samples_mean, desc(mean_expression))
all_geve_gene_counts_normalized_full_samples_mean$rank = 1:nrow(all_geve_gene_counts_normalized_full_samples_mean)

# Create plot of the rank and mean expression of genes in GEVE+ samples highlighting the position of GEVEv2_173
geve_gene_rank_plot = ggplot(all_geve_gene_counts_normalized_full_samples_mean, aes(x = rank, y = mean_expression)) + 
  geom_point(alpha = 0.25, shape = 21, fill = "black", size = 2) +
  geom_point(data = all_geve_gene_counts_normalized_full_samples_mean["GEVEv2_173", ], mapping = aes(x = rank, y = mean_expression), shape = 21, fill = "red", size = 3) +
  geom_text(data = all_geve_gene_counts_normalized_full_samples_mean["GEVEv2_173", ], mapping = aes(x = rank + 10, y = mean_expression + 0.25, label = "TET"), size = 8)
geve_gene_rank_plot = customize_ggplot_theme(geve_gene_rank_plot, xlab = "Rank", ylab = "Mean Log10\nNormalized Counts", title = NULL) +
  scale_x_continuous(expand = c(0, 0), labels = scales::comma)
geve_gene_rank_plot
ggsave(plot = geve_gene_rank_plot, "plots/tet_rank_plot.pdf", width = 9, height = 9, device = cairo_pdf)

### Create boxplots showing the mean normalized counts in GEVE+ and GEVE- samples at day 3 and day5

# Get normalized counts for GEVE genes outside the deletion
non_deletion_geve_gene_counts_normalized = data.frame(counts(all_gene_counts_dds, normalized = T))[geve_genes, ]

# Get the mean normalized counts for GEVE genes in GEVE+ and GEVE- samples at day 3 and day5
non_deletion_geve_gene_counts_normalized_summary = data.frame(
  deletion_d3 = rowMeans(dplyr::select(non_deletion_geve_gene_counts_normalized, matches("d3_7|d3_18"))),
  wt_d3 = rowMeans(dplyr::select(non_deletion_geve_gene_counts_normalized, matches("d3_8|d3_15"))),
  deletion_d5 = rowMeans(dplyr::select(non_deletion_geve_gene_counts_normalized, matches("d5_7|d5_18"))),
  wt_d5 = rowMeans(dplyr::select(non_deletion_geve_gene_counts_normalized, matches("d5_8|d5_15")))
)

# Add a column with methylation status of genes
gene_methylation_list = readRDS("../long_read_analysis/gene_methylation_list.rds")
gene_methylation_vec = unlist(with(reshape2::melt(gene_methylation_list), split(L1, value)))
non_deletion_geve_gene_counts_normalized_summary$methylation_status = 
  factor(gene_methylation_vec[row.names(non_deletion_geve_gene_counts_normalized_summary)], 
    levels = c("mAT", "mCG", "Dual", "Unmarked"))

# Convert to long format and add columns with day and deletion status
non_deletion_geve_gene_counts_normalized_summary_long = tidyr::pivot_longer(non_deletion_geve_gene_counts_normalized_summary, 
  -methylation_status, names_to = "Status", values_to = "normalized_counts")
non_deletion_geve_gene_counts_normalized_summary_long$day = ifelse(grepl("d3", non_deletion_geve_gene_counts_normalized_summary_long$Status), "D3", "D5")
non_deletion_geve_gene_counts_normalized_summary_long$sample = gsub("_.*", "", non_deletion_geve_gene_counts_normalized_summary_long$Status)
non_deletion_geve_gene_counts_normalized_summary_long$sample = factor(ifelse(non_deletion_geve_gene_counts_normalized_summary_long$sample == "deletion", "GEVE-", "GEVE+"), c("GEVE+", "GEVE-"))

# Create boxlplots
gene_boxplots = ggplot(non_deletion_geve_gene_counts_normalized_summary_long, aes(x = sample, y = log10(normalized_counts+1), fill = methylation_status)) +
  geom_boxplot()
gene_boxplots = customize_ggplot_theme(gene_boxplots, 
  xlab = "Status", ylab = "Mean Log10\nNormalized Counts", fill_title = "Methylation\nStatus",
  fill_colors = c("#C49E85", "#710627", "#011936", "#708090"),
  facet = "day", facet_nrow = 1, facet_scales = "free_x") +
  theme(strip.background = element_rect(fill = "grey95"), legend.position = "top")
gene_boxplots 
ggsave(plot = gene_boxplots, "plots/geve_gene_boxplots.pdf", width = 9, height = 9, device = cairo_pdf)

### Plot DEG results

# Remove GEVE deletion genes from all_gene_counts and create a DESeqDataSet
non_deletion_gene_counts = all_gene_counts[setdiff(chlamy_genes, deletion_genes), ]
non_deletion_gene_counts_dds = DESeq(DESeqDataSetFromMatrix(countData = non_deletion_gene_counts, colData = deseq_metadata, design = ~ day + geve_status))

# Get results and indicate if changes are significant (p-value < 0.05 and absolute fold change > 1)
deletion_vs_full_geve_deseq_results = data.frame(results(non_deletion_gene_counts_dds, contrast = c("geve_status", "Deletion", "Full")))
deletion_vs_full_geve_deseq_results$significant = deletion_vs_full_geve_deseq_results$padj < 0.05 & abs(deletion_vs_full_geve_deseq_results$log2FoldChange) > 1
data.table::fwrite(tibble::rownames_to_column(deletion_vs_full_geve_deseq_results, "gene_id"), 
  "deletion_vs_full_geve_deseq_results.tsv.gz", sep = "\t")

# Add methylation status of genes to results
deletion_vs_full_geve_deseq_results$methylation_status = factor(gene_methylation_vec[row.names(deletion_vs_full_geve_deseq_results)], c("mCG", "mAT", "Dual", "Unmarked"))

# Add columns with -log10 adjusted p-values and chromosome name
deletion_vs_full_geve_deseq_results$minus_log10_padj = -log10(deletion_vs_full_geve_deseq_results$padj)
deletion_vs_full_geve_deseq_results$chromosome = factor(as.character(seqnames(chlamy_genes_gr[row.names(deletion_vs_full_geve_deseq_results)])), levels = seqlevels(chlamy_genes_gr))

# Filter for significant results and add direction of change
deletion_vs_full_geve_deseq_results_significant = filter(deletion_vs_full_geve_deseq_results[non_geve_genes, ], significant)
deletion_vs_full_geve_deseq_results_significant$direction = factor(ifelse(deletion_vs_full_geve_deseq_results_significant$log2FoldChange > 0, "Upregulated", "Downregulated"))

# Summarize the number of significant genes by chromosome, methylation status and direction 
deletion_vs_full_geve_deseq_results_significant_summary = summarise(group_by(deletion_vs_full_geve_deseq_results_significant, 
  chromosome, methylation_status, direction), count = n())
deletion_vs_full_geve_deseq_results_significant_summary$direction = factor(deletion_vs_full_geve_deseq_results_significant_summary$direction)
deletion_vs_full_geve_deseq_results_significant_summary = tidyr::complete(ungroup(deletion_vs_full_geve_deseq_results_significant_summary), 
  chromosome, methylation_status, direction, fill = list(count = 0))
deletion_vs_full_geve_deseq_results_significant_summary = filter(deletion_vs_full_geve_deseq_results_significant_summary, chromosome != "cpDNA")

# Create a barplot with the number of DEGs per chromosome
de_gene_barplot = ggplot(deletion_vs_full_geve_deseq_results_significant_summary, aes(x = chromosome, y = count, fill = direction)) +
  geom_col(position = position_dodge2(width = 1, preserve = "single", padding = 0.1), color = "black") 
de_gene_barplot = customize_ggplot_theme(de_gene_barplot, facet = "methylation_status", x_labels = 1:17, xlab = "Chromosome",
  ylab = "Number of Genes", fill_colors = c("#4B878BFF", "#D01C1FFF"),facet_nrow = 4, facet_scales = "free") + 
  theme(strip.background = element_rect(fill = "grey95"), axis.text.x = element_text(size = 12)) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) +
  theme(legend.position = "top")
de_gene_barplot
ggsave(plot = de_gene_barplot, "plots/de_gene_barplot.pdf", width = 9, height = 9, device = cairo_pdf)

# Create a volcano plot for DE GEVE genes
geve_volcano_plot = ggplot(deletion_vs_full_geve_deseq_results[geve_genes, ], aes(x = log2FoldChange, y = minus_log10_padj, fill = methylation_status)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = 1.30103, linetype = "dashed") +
  geom_point(shape = 21, size = 5, alpha = ifelse(deletion_vs_full_geve_deseq_results[geve_genes, ]$significant, 1, 0.25))
geve_volcano_plot = customize_ggplot_theme(geve_volcano_plot, 
  xlab = "Log2 Fold Change", ylab = "-Log10 Adjusted p-value", title = "Differentially Expressed Genes in GEVE",
  fill_colors = c("#C49E85", "#710627", "#011936", "#708090"))
geve_volcano_plot
ggsave(plot = geve_volcano_plot, "plots/geve_de_genes_volcano.pdf", width = 9, height = 9, device = cairo_pdf)

# Create a volcano plot for DEGs on chromosome 15 (excluding those in the GEVE)
chr15_volcano_plot = ggplot(deletion_vs_full_geve_deseq_results[chr15_genes, ], aes(x = log2FoldChange, y = minus_log10_padj, fill = methylation_status)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = 1.30103, linetype = "dashed") +
  geom_point(shape = 21, size = 5, alpha = ifelse(deletion_vs_full_geve_deseq_results[chr15_genes, ]$significant, 1, 0.25))
chr15_volcano_plot = customize_ggplot_theme(chr15_volcano_plot, 
  xlab = "Log2 Fold Change", ylab = "-Log10 Adjusted p-value", , title = "Differentially Expressed Genes on Chromosome 15",
  fill_colors = c("#C49E85", "#710627", "#011936", "#708090"))
chr15_volcano_plot
ggsave(plot = chr15_volcano_plot, "plots/chr15_de_genes_volcano.pdf", width = 9, height = 9, device = cairo_pdf)