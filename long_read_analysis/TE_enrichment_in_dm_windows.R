# Test enrichment of TE subfamilies in hypermethylated and hypomethylated windows

# Load required packages
library(dplyr)
library(regioneR)
source("../auxillary_scripts/ggplot_functions.R")

# Create a GRanges for the Chlamy genome
chlamy_genome_gr = GRanges(seqinfo(Biostrings::readDNAStringSet("../genome_files/CC2937_T2T.fa")))[1:17]

# Load GRanges with GEVE and relicts
geve_and_relicts_gr = rtracklayer::import.bed("../genome_files/GEVE_and_relicts_T2T.bed")

# Load GRanges with TEs and subset for those outside GEVE
chlamy_tes_gr = rtracklayer::import.gff("../genome_files/CC2937_T2T.TEs_v3_8.bed.gtf")
chlamy_tes_gr = subsetByOverlaps(chlamy_tes_gr, geve_and_relicts_gr, invert = T)

# Make a table with class, family and subfamily for each TE
te_metadata = data.frame(class = chlamy_tes_gr$class_id, family = chlamy_tes_gr$family_id, subfamily = chlamy_tes_gr$gene_id, row.names = names(chlamy_tes_gr))
te_metadata = filter(te_metadata, !duplicated(subfamily))

# Make vectors for matching TE subfamily to family and class
te_subfamily_to_family = setNames(te_metadata$family, te_metadata$subfamily)
te_subfamily_to_class = setNames(te_metadata$class, te_metadata$subfamily)

# Set the number of permutations to run
nperm = 10000

# Import hypermethylated and hypomethylated windows in GEVE_ samples and subset for elements outside the GEVE/relicts
hyper_windows = subsetByOverlaps(rtracklayer::import.bed("deletion_hypermethylated_windows.bed"), geve_and_relicts_gr, invert = T)
hypo_windows = subsetByOverlaps(rtracklayer::import.bed("deletion_hypomethylated_windows.bed"), geve_and_relicts_gr, invert = T)

### Test enrichment of repeat subfamilies in hypermethylated windows

# Calculate the number of TEs in each subfamily overlapping hypermethylated windows
hyper_subfamily_overlaps = data.frame(table(subsetByOverlaps(chlamy_tes_gr, hyper_windows)$gene_id))
hyper_subfamily_overlaps = tibble::column_to_rownames(hyper_subfamily_overlaps, "Var1")

# Create a GRangeList shuffling hyper_windows across genome (except GEVE)
set.seed(123)
system.time({random_hyper_windows = GRangesList(parallel::mclapply(seq.int(nperm), function(x)
    randomizeRegions(A = hyper_windows, genome = chlamy_genome_gr, 
      allow.overlaps = F, per.chromosome = F, mask = geve_and_relicts_gr), mc.cores = 12))})
saveRDS(random_hyper_windows, "random_hyper_windows.rds")
random_hyper_windows = readRDS("random_hyper_windows.rds")

# Name GRanges with permuation number and convert to a flat GRanges
names(random_hyper_windows) = paste0("permutation_", seq_along(random_hyper_windows))
random_hyper_windows_flat = unlist(random_hyper_windows)

# Find overlaps between random_hyper_windows and repeat families and subfamilies
random_hyper_windows_repeat_overlaps = data.frame(findOverlaps(random_hyper_windows_flat, chlamy_tes_gr))
random_hyper_windows_repeat_overlaps$permutation = names(random_hyper_windows_flat)[random_hyper_windows_repeat_overlaps$queryHits]
random_hyper_windows_repeat_overlaps$family = chlamy_tes_gr$family_id[random_hyper_windows_repeat_overlaps$subjectHits]
random_hyper_windows_repeat_overlaps$subfamily = chlamy_tes_gr$gene_id[random_hyper_windows_repeat_overlaps$subjectHits]

# Count number of hits for each subfamily of repeat for each permutation
random_hyper_windows_subfamily_summary = 
  summarise(group_by(random_hyper_windows_repeat_overlaps, permutation, subfamily), count = n())

# Convert permutation and subfamily to factors 
random_hyper_windows_subfamily_summary$permutation = 
  factor(random_hyper_windows_subfamily_summary$permutation, levels = names(random_hyper_windows))
random_hyper_windows_subfamily_summary$subfamily = 
  factor(random_hyper_windows_subfamily_summary$subfamily, levels = sort(unique(chlamy_tes_gr$gene_id)))

# Convert random_hyper_windows_family_summary to wide format
random_hyper_windows_subfamily_summary = tidyr::pivot_wider(random_hyper_windows_subfamily_summary, 
  names_from = permutation, values_from = count, values_fill = 0, names_expand = T, id_expand = T)
random_hyper_windows_subfamily_summary = tibble::column_to_rownames(random_hyper_windows_subfamily_summary, "subfamily")

# Ensure repeats are in same order as hypermethylated_repeat_overlaps
random_hyper_windows_subfamily_summary = random_hyper_windows_subfamily_summary[row.names(hyper_subfamily_overlaps), ]

# Indicate if each random permutation has the same or more overlaps with each subfamily of repeat 
expected_exceeds_observed_hyper_window_subfamily = random_hyper_windows_subfamily_summary >= hyper_subfamily_overlaps[, 1]

# Calculate the number of times the randomized regions have an equal or greater number of overlaps than the actual hypermethylated regions
hyper_window_te_subfamily_perm_test_results = data.frame(
  subfamily = row.names(expected_exceeds_observed_hyper_window_subfamily),
  observed_overlaps = hyper_subfamily_overlaps[, 1],
  mean_random_overlaps = rowMeans(random_hyper_windows_subfamily_summary),
  expected_exceeds_observed_counts = rowSums(expected_exceeds_observed_hyper_window_subfamily),
  row.names = NULL)

# Add p-values 
hyper_window_te_subfamily_perm_test_results$p_value = (hyper_window_te_subfamily_perm_test_results$expected_exceeds_observed_counts + 1)/(nperm+1)

### Test enrichment of repeat subfamilies in hypomethylated windows

# Calculate the number of TEs in each subfamily overlapping hypomethylated windows
hypo_subfamily_overlaps = data.frame(table(subsetByOverlaps(chlamy_tes_gr, hypo_windows)$gene_id))
hypo_subfamily_overlaps = tibble::column_to_rownames(hypo_subfamily_overlaps, "Var1")

# Create a GRangeList shuffling hypo_windows across genome (except GEVE)
set.seed(123)
system.time({random_hypo_windows = GRangesList(parallel::mclapply(seq.int(nperm), function(x)
    randomizeRegions(A = hypo_windows, genome = chlamy_genome_gr, 
      allow.overlaps = F, per.chromosome = F, mask = geve_and_relicts_gr), mc.cores = 12))})
saveRDS(random_hypo_windows, "random_hypo_windows.rds")
random_hypo_windows = readRDS("random_hypo_windows.rds")

# Name GRanges with permuation number and convert to a flat GRanges
names(random_hypo_windows) = paste0("permutation_", seq_along(random_hypo_windows))
random_hypo_windows_flat = unlist(random_hypo_windows)

# Find overlaps between random_hypo_windows and repeat families and subfamilies
random_hypo_windows_repeat_overlaps = data.frame(findOverlaps(random_hypo_windows_flat, chlamy_tes_gr))
random_hypo_windows_repeat_overlaps$permutation = names(random_hypo_windows_flat)[random_hypo_windows_repeat_overlaps$queryHits]
random_hypo_windows_repeat_overlaps$family = chlamy_tes_gr$family_id[random_hypo_windows_repeat_overlaps$subjectHits]
random_hypo_windows_repeat_overlaps$subfamily = chlamy_tes_gr$gene_id[random_hypo_windows_repeat_overlaps$subjectHits]

# Count number of hits for each subfamily of repeat for each permutation
random_hypo_windows_subfamily_summary = 
  summarise(group_by(random_hypo_windows_repeat_overlaps, permutation, subfamily), count = n())

# Convert permutation and subfamily to factors 
random_hypo_windows_subfamily_summary$permutation = 
  factor(random_hypo_windows_subfamily_summary$permutation, levels = names(random_hypo_windows))
random_hypo_windows_subfamily_summary$subfamily = 
  factor(random_hypo_windows_subfamily_summary$subfamily, levels = sort(unique(chlamy_tes_gr$gene_id)))

# Convert random_hypo_windows_family_summary to wide format
random_hypo_windows_subfamily_summary = tidyr::pivot_wider(random_hypo_windows_subfamily_summary, 
  names_from = permutation, values_from = count, values_fill = 0, names_expand = T, id_expand = T)
random_hypo_windows_subfamily_summary = tibble::column_to_rownames(random_hypo_windows_subfamily_summary, "subfamily")

# Ensure repeats are in same order as hypomethylated_repeat_overlaps
random_hypo_windows_subfamily_summary = random_hypo_windows_subfamily_summary[row.names(hypo_subfamily_overlaps), ]

# Indicate if each random permutation has the same or more overlaps with each subfamily of repeat 
expected_exceeds_observed_hypo_window_subfamily = random_hypo_windows_subfamily_summary >= hypo_subfamily_overlaps[, 1]

# Calculate the number of times the randomized regions have an equal or greater number of overlaps than the actual hypomethylated regions
hypo_window_te_subfamily_perm_test_results = data.frame(
  subfamily = row.names(expected_exceeds_observed_hypo_window_subfamily),
  observed_overlaps = hypo_subfamily_overlaps[, 1],
  mean_random_overlaps = rowMeans(random_hypo_windows_subfamily_summary),
  expected_exceeds_observed_counts = rowSums(expected_exceeds_observed_hypo_window_subfamily),
  row.names = NULL)

# Add p-values 
hypo_window_te_subfamily_perm_test_results$p_value = (hypo_window_te_subfamily_perm_test_results$expected_exceeds_observed_counts + 1)/(nperm+1)

### Make plots combining hyper and hypo results

# Combine hyper and hypo results, correct p-values and add empty rows for missing subfamily/change combinations
combined_window_te_subfamily_perm_test_results = bind_rows(list(
  hyper = hyper_window_te_subfamily_perm_test_results,
  hypo = hypo_window_te_subfamily_perm_test_results), .id = "Change")
combined_window_te_subfamily_perm_test_results$q_value = p.adjust(combined_window_te_subfamily_perm_test_results$p_value, method = "fdr")
combined_window_te_subfamily_perm_test_results$significance = sig_sym(combined_window_te_subfamily_perm_test_results$q_value)
combined_window_te_subfamily_perm_test_results = tidyr::complete(combined_window_te_subfamily_perm_test_results, subfamily, Change, 
  fill = list(observed_overlaps = 0, mean_random_overlaps = 0, significance = ""))

# Filter combined_window_te_subfamily_perm_test_results just for significant results
combined_subfamily_overlap_enriched = filter(combined_window_te_subfamily_perm_test_results, significance != "")
combined_subfamily_overlap_enriched = tidyr::complete(combined_subfamily_overlap_enriched, subfamily, Change, 
  fill = list(observed_overlaps = 0, mean_random_overlaps = 0, significance = ""))
combined_subfamily_overlap_enriched$family = factor(te_subfamily_to_family[combined_subfamily_overlap_enriched$subfamily])
combined_subfamily_overlap_enriched$class = factor(te_subfamily_to_class[combined_subfamily_overlap_enriched$subfamily])

# Create barplot for subfamilies that are significant in either hyper or hypo windows
combined_subfamily_enriched_barplot = ggplot(combined_subfamily_overlap_enriched, aes(x = subfamily, y = observed_overlaps-mean_random_overlaps, fill = Change)) +
  geom_col(color = "black", position = position_dodge(width = 0.9)) +
  geom_text(mapping = aes(label = significance, y = observed_overlaps-mean_random_overlaps), size = 8, position = position_dodge(width = 0.9))
combined_subfamily_enriched_barplot = customize_ggplot_theme(combined_subfamily_enriched_barplot, 
  xlab = "Repeat Subfamily", ylab = "Observed - Expected Counts", 
  x_labels_angle = 45, axis_text_size = 16, 
  fill_colors = c("#D01C1FFF", "#4B878BFF"), fill_labels = c("Increase", "Decrease"), fill_title = "5mC Change") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  ggh4x::facet_nested(. ~ class + family, 
    scales = "free_x", space = "free_x", 
    strip = ggh4x::strip_nested(text_x = list(class = element_text(size = 12, face = "bold"), family = element_text(size = 12)), by_layer_x = T)) +
  theme(panel.spacing = unit(0.1, "lines"), strip.background = element_rect(fill = "grey95")) 
combined_subfamily_enriched_barplot
ggsave(plot = combined_subfamily_enriched_barplot, "plots/te_subfamily_dm_window_enrichment_barplot.pdf", width = 20, height = 11.25, device = cairo_pdf)