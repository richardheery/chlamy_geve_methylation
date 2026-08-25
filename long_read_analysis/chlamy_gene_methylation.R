# Define genes and TEs marked with mAT and mCG

# Load required packages and functions
library(methodical)
source("../helper_scripts/ggplot_functions.R")

# Get locations of genes and their TSS and names of GEVE genes
chlamy_genes_gr = readRDS("../genome_files/chlamy_genes_gr.rds")
chlamy_tss = promoters(chlamy_genes_gr, upstream = 0, downstream = 1)
geve_genes = rtracklayer::import.bed("../genome_files/GEVE_all_genes.bed")$name

### Define mAT genes

# Load RSE for 6mA and mask sites covered by < 10 or >= 350 reads
chlamy_6ma_rse = HDF5Array::loadHDF5SummarizedExperiment("CC2937_T2T_6mA_rse")
colnames(chlamy_6ma_rse) = c("G2", "D2", "D1", "G1")
assay(chlamy_6ma_rse, 1)[assay(chlamy_6ma_rse, 2) < 10 | assay(chlamy_6ma_rse, 2) >= 350] = NA

# Get methylation values for all AT sites across the genome and convert to long format
all_6ma_values = data.frame(assay(chlamy_6ma_rse))
all_6ma_values_long = data.frame(tidyr::pivot_longer(all_6ma_values, cols = everything(), names_to = "sample", values_to = "mAT"))

#  Plot distribution of AT methylation values in the 4 samples. 0.5 seems to separate methylated from unmethylated mAT sites
all_6ma_values_hist = ggplot(all_6ma_values_long, aes(x = mAT, fill = sample)) +
  geom_histogram(alpha = 0.5,  position = "identity", color = "black") 
all_6ma_values_hist = customize_ggplot_theme(all_6ma_values_hist, xlab = "Proportion Methylated", ylab = "Number of AT Sites", 
  fill_colors = c("#FEE08B", "#ABDDA4", "#3288BD", "#9E0142")) + 
  scale_x_log10(labels = scales::label_number(drop0trailing=TRUE), breaks = c(0.005, 0.01, 0.025, 0.05, 0.1, 0.2, 0.5, 1)) +
  theme(legend.position = "top")
all_6ma_values_hist
ggsave(plot = all_6ma_values_hist, "plots/6mA_distribution_all_sites.pdf", width = 9, height = 9)

# Define mAT sites as those with > 0.5 methylation both chlamy_8A and chlamy_15A
mAT_sites = rowRanges(chlamy_6ma_rse)[which(rowMins(assay(chlamy_6ma_rse[, c("G1", "G2")])) > 0.5)]
saveRDS(mAT_sites, "mAT_sites_whole_genome.rds")

# Find distribution of distances of mAT sites to their nearest TSS
tss_dist_to_mAT_sites_genes = strandedDistance(mAT_sites, chlamy_tss[nearest(mAT_sites, chlamy_tss)])

# 75% of mAT sites are within -400 to +1000 of TSS
hist(tss_dist_to_mAT_sites_genes[abs(tss_dist_to_mAT_sites_genes) < 2000])
prop.table(table(tss_dist_to_mAT_sites_genes > -400 & tss_dist_to_mAT_sites_genes < 1000))

# Define promoters as -400 to +1000 bp (or the end of the gene for genes shorter than 1000 bp)
chlamy_promoters_gr = expand_granges(chlamy_tss, upstream = 400, 
  downstream = pmin(width(chlamy_genes_gr) - 1, 999))

# Define mAT genes as those with promoters overlapping mAT sites
mAT_genes = names(subsetByOverlaps(chlamy_promoters_gr, mAT_sites))

### Define genes which gain mCG in deletion samples

# Load RSE for 5mC and mask sites covered by < 10 or >= 350 reads
chlamy_5mc_rse = HDF5Array::loadHDF5SummarizedExperiment("CC2937_T2T_5mC_rse/")
colnames(chlamy_5mc_rse) = c("G2", "D2", "D1", "G1")
assay(chlamy_5mc_rse, 1)[assay(chlamy_5mc_rse, 2) < 10 | assay(chlamy_5mc_rse, 2) >= 350] = NA

# Get combined promoter and gene body region for each gene
chlamy_genes_and_promoters_gr = expand_granges(chlamy_genes_gr, upstream = 400)

# Get mean 5mC for combined gene and promoter regions
gene_mean_5mC = summarizeRegionMethylation(chlamy_5mc_rse, genomic_regions = chlamy_genes_and_promoters_gr)

# Get change in 5mC between deletion and WT samples
gene_mean_5mC$deletion_change = ((gene_mean_5mC$D1 + gene_mean_5mC$D2)/2) - 
  ((gene_mean_5mC$G1 + gene_mean_5mC$G2)/2)

# Define mCG genes as those gaining > 0.15 5mC in gene and promoter regions in deletion samples
# 80% are in the GEVE. The majority of the non-GEVE mCG genes are still on chromosome 15
hist(gene_mean_5mC[geve_genes, ]$deletion_change)
mCG_genes_deletion = row.names(gene_mean_5mC)[which(gene_mean_5mC$deletion_change > 0.15)]

### Export sets of genes based on mAT and mCG status

# Define unmarked genes as those without mCG or mAT
unmarked_genes = setdiff(names(chlamy_genes_gr), c(mCG_genes_deletion, mAT_genes))

# Define dual genes as those with mAT and that gain mCG in deletion samples
dual_genes = intersect(mCG_genes_deletion, mAT_genes)

# Remove dual genes from mCG and mAT genes
mCG_genes = setdiff(mCG_genes_deletion, dual_genes)
mAT_genes = setdiff(mAT_genes, dual_genes)

# Export a list with the different categories of genes by methylation status
gene_methylation_list = list(mAT = mAT_genes, mCG = mCG_genes, Dual = dual_genes, Unmarked = unmarked_genes)
saveRDS(gene_methylation_list, "gene_methylation_list.rds")

### Define mAT and mCG TEs

# Load TE coordinates from GTF
chlamy_tes_gr = rtracklayer::import.gff2("../genome_files/CC2937_T2T.TEs_v3_8.bed.gtf")
names(chlamy_tes_gr) = paste(chlamy_tes_gr$transcript_id, chlamy_tes_gr$gene_id, chlamy_tes_gr$family_id, chlamy_tes_gr$class_id, sep = ":")

# Extract mean 5mC for TEs
te_mean_5mC = summarizeRegionMethylation(chlamy_5mc_rse, genomic_regions = chlamy_tes_gr)

# Calculate change between deletion and WT samples
te_mean_5mC$deletion_change = ((te_mean_5mC$D1 + te_mean_5mC$D2)/2) - 
  ((te_mean_5mC$G1 + te_mean_5mC$G2)/2)

# Define TEs gaining 5mC in deletion samples using 0.15 as threshold. 
mCG_tes_deletion = row.names(te_mean_5mC)[which(te_mean_5mC$deletion_change > 0.15)]

# Find TEs overlapping mAT sites
mAT_tes = names(subsetByOverlaps(chlamy_tes_gr, mAT_sites))

# Define unmarked TEs as those without mCG or mAT
unmarked_tes = setdiff(names(chlamy_tes_gr), c(mCG_tes_deletion, mAT_tes))

# Define dual TEs as those with mAT and that gain mCG in deletion samples
dual_tes = intersect(mCG_tes_deletion, mAT_tes)

# Remove dual TEs from mCG and mAT TEs
mCG_tes = setdiff(mCG_tes_deletion, dual_tes)
mAT_tes = setdiff(mAT_tes, dual_tes)

# Export a list with the different categories of TEs by methylation status
te_methylation_list = list(mAT = mAT_tes, mCG = mCG_tes, Dual = dual_tes, Unmarked = unmarked_tes)
saveRDS(te_methylation_list, "te_methylation_list.rds")