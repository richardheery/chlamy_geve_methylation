# Test enrichment of KEGG pathways among differentially expressed genes

# Load required packages
library(dplyr)
source("../helper_scripts/ggplot_functions.R")
source("../helper_scripts/enrichment_tests.R")

# Import DESeq2 results
deletion_vs_full_geve_deseq_results = data.frame(data.table::fread("deletion_vs_full_geve_deseq_results.tsv.gz"), row.names = 1)

# Get names of significant DEGs and genes that were actually tested
significant_genes = row.names(filter(deletion_vs_full_geve_deseq_results, significant))
significant_upregulated_genes = row.names(filter(deletion_vs_full_geve_deseq_results, significant, log2FoldChange > 0))
significant_downregulated_genes = row.names(filter(deletion_vs_full_geve_deseq_results, significant, log2FoldChange < 0))

# Get names of genes that were actually tested
tested_genes = row.names(filter(deletion_vs_full_geve_deseq_results, !is.na(pvalue)))

# Import list with genes in KEGG pathways and all genes with at least one KEGG pathway
genes_per_kegg_pathway = readRDS("../genome_files/kegg_pathway_genes_list.rds")
genes_per_kegg_pathway = genes_per_kegg_pathway[lengths(genes_per_kegg_pathway) >= 10]
kegg_genes = unique(unlist(genes_per_kegg_pathway))

# Define the universe as all tested genes which were associated with at least one KEGG pathway
universe_genes = intersect(kegg_genes, tested_genes)

# Get names of significant upregulated and downregulated genes associated with KEGG pathways
significant_upregulated_genes_kegg = intersect(significant_upregulated_genes, kegg_genes)
significant_downregulated_genes_kegg = intersect(significant_downregulated_genes, kegg_genes)

# Test enrichment of KEGG pathways among upregulated genes and rank the results by relative_enrichment
upregulated_genes_kegg_pathway_enrichment = fisher_test_apply(test = significant_upregulated_genes_kegg, query_list = genes_per_kegg_pathway, universe = universe_genes)
upregulated_genes_kegg_pathway_enrichment = arrange(upregulated_genes_kegg_pathway_enrichment, relative_enrichment)
upregulated_genes_kegg_pathway_enrichment$query_name = factor(upregulated_genes_kegg_pathway_enrichment$query_name, levels = upregulated_genes_kegg_pathway_enrichment$query_name)
upregulated_genes_kegg_pathway_enrichment$ranking = 1:nrow(upregulated_genes_kegg_pathway_enrichment)

# Test enrichment of KEGG pathways among downregulated genes and rank the results by relative_enrichment
downregulated_genes_kegg_pathway_enrichment = fisher_test_apply(test = significant_downregulated_genes_kegg, query_list = genes_per_kegg_pathway, universe = universe_genes)
downregulated_genes_kegg_pathway_enrichment = arrange(downregulated_genes_kegg_pathway_enrichment, relative_enrichment)
downregulated_genes_kegg_pathway_enrichment$query_name = factor(downregulated_genes_kegg_pathway_enrichment$query_name, levels = downregulated_genes_kegg_pathway_enrichment$query_name)
downregulated_genes_kegg_pathway_enrichment$ranking = 1:nrow(downregulated_genes_kegg_pathway_enrichment)

# Combine pathway enrichment results for upregulated and downregulated genes
combined_kegg_pathway_enrichment = bind_rows(list(
  Upregulated = upregulated_genes_kegg_pathway_enrichment, Downregulated = downregulated_genes_kegg_pathway_enrichment), .id = "Direction")
combined_kegg_pathway_enrichment$Direction = factor(combined_kegg_pathway_enrichment$Direction, levels = c("Upregulated", "Downregulated"))

# Recorrect p-values for combined results and filter for significant pathways
combined_kegg_pathway_enrichment$q_value = p.adjust(combined_kegg_pathway_enrichment$p_value, method = "fdr")
combined_kegg_pathway_enrichment = filter(combined_kegg_pathway_enrichment, q_value < 0.05)

# Create plot of significantly enriched pathways for upregulated and downregulated genes
degs_kegg_pathway_enrichment_plot = ggplot(combined_kegg_pathway_enrichment, 
  aes(x = relative_enrichment, y = tidytext::reorder_within(query_name, ranking, Direction), color = q_value, size = test_overlap_size)) + 
  geom_point() + 
  geom_vline(xintercept = 1, linetype = "dashed")
degs_kegg_pathway_enrichment_plot = customize_ggplot_theme(degs_kegg_pathway_enrichment_plot, 
  xlab = "Relative Enrichment", ylab = "KEGG Pathway", color_title = "Adjusted\np-value",
  facet = "Direction", facet_ncol = 1, facet_scales = "free_y", axis_title_size = 14,
  legend_title_size = 14, legend_text_size = 10, legend_key_size = 0.5) +
  theme(axis.text.y = element_text(size = 10), strip.background = element_blank(), strip.text = element_text(size = 14)) +
  labs(size = "Overlap Size") +
  scale_colour_continuous(guide = guide_colorbar(order = 1, reverse = T)) +
  tidytext::scale_y_reordered() 
degs_kegg_pathway_enrichment_plot
ggsave(plot = degs_kegg_pathway_enrichment_plot, "plots/degs_kegg_pathway_enrichment_plot.pdf",
  width = 7.5, height = 7, device = cairo_pdf())