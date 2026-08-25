# Make boxplots of Chip-seq signal for hypermethylated, hypomethylated and unchanged windows (outside GEVE)

# Load required packages and functions
library(dplyr)
library(ggplot2)
library(GenomicRanges)
source("../helper_scripts/bigwig_summarize_over_regions.R")
source("../helper_scripts/ggplot_functions.R")

# Load GRanges with GEVE and relicts and with centromeres
geve_and_relicts_gr = rtracklayer::import.bed("../genome_files/GEVE_and_relicts_T2T.bed")

# Get hypermethylated, hypomethylated and unchanged windows and remove regions overlapping GEVE or relicts
hypermethylated_windows = subsetByOverlaps(rtracklayer::import.bed("../long_read_analysis/deletion_hypermethylated_windows.bed"), 
  geve_and_relicts_gr, invert = T)
names(hypermethylated_windows) = paste("Hypermethylated", seq_along(hypermethylated_windows), sep = "_")
hypomethylated_windows = subsetByOverlaps(rtracklayer::import.bed("../long_read_analysis/deletion_hypomethylated_windows.bed"), 
  geve_and_relicts_gr, invert = T)
names(hypomethylated_windows) = paste("Hypomethylated", seq_along(hypomethylated_windows), sep = "_")
unchanged_windows = subsetByOverlaps(rtracklayer::import.bed("../long_read_analysis/deletion_unchanged_windows.bed"), 
  geve_and_relicts_gr, invert = T)
names(unchanged_windows) = paste("Unchanged", seq_along(unchanged_windows), sep = "_")

# Take a random sample of the unchanged windows for plotting
set.seed(123)
unchanged_windows = sample(unchanged_windows, 250)

# Combine GRanges for hypermethylated, hypomethylated and unchanged windows
combined_windows = c(hypermethylated_windows, hypomethylated_windows, unchanged_windows)
combined_windows$name = names(combined_windows)

# Get paths to all coverage bigwigs
coverage_bigwigs = list.files("coverage_bigwig_files", full.names = T, pattern = ".bw")

# Remove windows with no coverage across all coverage bigwigs
window_coverage_values = bigwig_summarize_over_regions(bw_filepaths = coverage_bigwigs, regions_granges = combined_windows, 
  parallel_files = 3, statistic = "mean0")
combined_windows = combined_windows[row.names(window_coverage_values)[rowMeans(window_coverage_values) > 0]]

# Get paths to all normalized bigWigs
normalized_bigwigs = list.files("normalized_bigwigs", full.names = T)

# Select just one H3K4me3 sample (SRR1521636)
normalized_bigwigs = normalized_bigwigs[c(grep("H3K4me3", normalized_bigwigs, invert = T), grep("SRR1521636", normalized_bigwigs))]
bigwig_names = gsub("_.*", "", basename(normalized_bigwigs))
bigwig_names[bigwig_names == "RNA"] = "RNA Pol II"
histone_to_effects = setNames(c("Activating", "Silencing", "Activating", "Silencing", "Activating", "Activating", "Activating"), bigwig_names)

# Get the mean histone mark values for the windows. Replacing missing regions with 0
window_chip_seq_values = bigwig_summarize_over_regions(bw_filepaths = normalized_bigwigs, regions_granges = combined_windows, 
  column_names = bigwig_names, parallel_files = 3, statistic = "mean0")
window_chip_seq_values$window = gsub("_.*", "", row.names(window_chip_seq_values))

# Convert to long format
window_chip_seq_values = tidyr::pivot_longer(window_chip_seq_values, cols = -window, values_to = "score", names_to = "Mark")

# Classify histone marks
window_chip_seq_values$classification = histone_to_effects[window_chip_seq_values$Mark]

# Remove RNA pol II
window_chip_seq_values = filter(window_chip_seq_values, Mark != "RNA Pol II")

# Make boxplots for histone mark singal across windows
chip_boxplots = ggplot(window_chip_seq_values, aes(x = Mark, y = score, fill = window)) +
  geom_boxplot()
chip_boxplots = customize_ggplot_theme(chip_boxplots, ylab = "log2(Chip Vs Input)",  fill_title = "5mC Change", 
  show_legend = T, fill_labels = c("Increase", "Decrease", "Unchanged"), fill_colors = c("#CD2626", "#53868B", "grey")) +
  facet_grid(. ~  classification, scales = "free_x", space = "free_x") +
    geom_hline(yintercept = 0, linetype = "dashed") +
  theme(strip.background = element_rect(fill = "grey95"), 
    legend.position = c(0.9, 0.15), legend.background = element_rect(colour = "black", fill = "white", linewidth = 0.5))
chip_boxplots 
ggsave(plot = chip_boxplots, "plots/chlamy_chip_boxplots.pdf", width = 16, height = 9)

# Make boxplots with significance for Wilcoxon tests added

# Define comparisons to test
comparisons <- list(
  c("Hypermethylated", "Hypomethylated"),
  c("Hypermethylated", "Unchanged"),
  c("Hypomethylated", "Unchanged")
)

# Make boxplots with significance
chip_boxplots_test = ggplot(window_chip_seq_values, aes(x = window, y = score, fill = window)) +
  geom_boxplot()
chip_boxplots_test = customize_ggplot_theme(chip_boxplots_test, 
  ylab = "log2(Chip Vs Input)",  fill_title = "5mC Change", show_legend = T,
  fill_colors = c("#CD2626", "#53868B", "grey")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  ggh4x::facet_grid2(~ classification + Mark, scales = "free_x", space = "free_x", strip = ggh4x::strip_nested(bleed = T)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
  theme(strip.background = element_rect(fill = "grey95")) +
  ggpubr::geom_signif(
    comparisons = comparisons,
    test = "wilcox.test",
    map_signif_level = TRUE,
    y_position = c(2.5, 3.25, 4),
    tip_length = 0.03
  )
chip_boxplots_test 
ggsave(plot = chip_boxplots_test, "plots/chlamy_chip_boxplots_wilcoxon_tests.pdf", width = 16, height = 9)