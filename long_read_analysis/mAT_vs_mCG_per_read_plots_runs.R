# Plot the per-read levels of 5mCG vs 6mAT

# Load required packages and functions
library(data.table)
library(dplyr)
library(ggplot2)
library(GenomicRanges)
source("../helper_scripts/ggplot_functions.R")
source("../helper_scripts/modkit_functions.R")

# Load GRanges with GEVE
geve_gr = rtracklayer::import.bed("../genome_files/GEVE_T2T.bed")

# Get CG mod calls for third_run and filter them for quality calls
third_run_5mCG_mod_calls = fread("modkit_calls/ThirdRun_GEVEs_5mCG_modkit_calls.tsv.gz"); gc()
third_run_5mCG_mod_calls = add_dinucleotide_context_to_modkit_table(third_run_5mCG_mod_calls); gc()
third_run_5mCG_mod_calls = filter_modkit_calls(third_run_5mCG_mod_calls, sequence_context = "CG"); gc()

# Convert calls to a GRanges and separated based on overlap with GEVE
third_run_5mCG_mod_calls = make_granges_from_mod_calls(third_run_5mCG_mod_calls)
third_run_5mCG_mod_calls_geve = subsetByOverlaps(third_run_5mCG_mod_calls, geve_gr)
third_run_5mCG_mod_calls_host = subsetByOverlaps(third_run_5mCG_mod_calls, geve_gr, invert = T)

# Get AT mod calls for third_run and filter them for quality calls
third_run_6mA_mod_calls = fread("modkit_calls/ThirdRun_GEVEs_6mAT_modkit_calls.tsv.gz")
third_run_6mA_mod_calls = add_dinucleotide_context_to_modkit_table(third_run_6mA_mod_calls)
third_run_6mA_mod_calls = filter_modkit_calls(modkit_table = third_run_6mA_mod_calls, sequence_context = "AT")

# Convert calls to a GRanges and separated based on overlap with GEVE
third_run_6mA_mod_calls = make_granges_from_mod_calls(third_run_6mA_mod_calls)
third_run_6mA_mod_calls_geve = subsetByOverlaps(third_run_6mA_mod_calls, geve_gr)
third_run_6mA_mod_calls_host = subsetByOverlaps(third_run_6mA_mod_calls, geve_gr, invert = T)

# Find reads with both AT and CG for GEVE and host 
common_reads_geve = intersect(third_run_6mA_mod_calls_geve$read_id, third_run_5mCG_mod_calls_geve$read_id)
common_reads_host = intersect(third_run_6mA_mod_calls_host$read_id, third_run_5mCG_mod_calls_host$read_id)

# Summarize 5mC calls per read
third_run_5mC_calls_per_read_summary_geve = summarize_mod_calls_per_read(filter(data.frame(third_run_5mCG_mod_calls_geve), 
  read_id %in% common_reads_geve))
third_run_5mC_calls_per_read_summary_host = summarize_mod_calls_per_read(filter(data.frame(third_run_5mCG_mod_calls_host),
  read_id %in% common_reads_host))

# Summarize 6mA calls per read
third_run_6mA_calls_per_read_summary_geve = summarize_mod_calls_per_read(filter(data.frame(third_run_6mA_mod_calls_geve), 
  read_id %in% common_reads_geve))
third_run_6mA_calls_per_read_summary_host = summarize_mod_calls_per_read(filter(data.frame(third_run_6mA_mod_calls_host),
  read_id %in% common_reads_host))

# Mask reads with less than 10 calls
third_run_5mC_calls_per_read_summary_host$m[third_run_5mC_calls_per_read_summary_host$total_calls < 10] = NA
third_run_6mA_calls_per_read_summary_host$m[third_run_6mA_calls_per_read_summary_host$total_calls < 10] = NA
third_run_5mC_calls_per_read_summary_geve$m[third_run_5mC_calls_per_read_summary_geve$total_calls < 10] = NA
third_run_6mA_calls_per_read_summary_geve$m[third_run_6mA_calls_per_read_summary_geve$total_calls < 10] = NA

# Combine 5mC and 6mA read values for host
host_read_combined = data.frame(
  m = third_run_5mC_calls_per_read_summary_host$m,
  a = third_run_6mA_calls_per_read_summary_host$a
)

# Combine 5mC and 6mA read values for GEVE
geve_read_combined = data.frame(
  m = third_run_5mC_calls_per_read_summary_geve$m,
  a = third_run_6mA_calls_per_read_summary_geve$a
)

# Combine read data for host and GEVE
all_reads = bind_rows(
  list(GEVE = geve_read_combined, Host = host_read_combined), 
  .id = "origin"
)
all_reads$origin = factor(all_reads$origin, levels = c("Host", "GEVE"))

# Make scatter plot for per-read levels of 5mCG vs 6mAT in host and GEVE
per_read_plot_geve_reads = ggplot(geve_read_combined, aes(x = m*100, y = a*100)) +
  geom_point(alpha = 0.25)
per_read_plot_geve_reads = customize_ggplot_theme(per_read_plot_geve_reads, xlab = "Read-Level mCG %", ylab = "Read-Level mAT %") +
  theme(strip.background = element_rect(fill = "grey95"), aspect.ratio = 1)
per_read_plot_geve_reads
ggsave(plot = per_read_plot_geve_reads, "plots/per_read_plot_geve_reads_erazo_3rd_run.pdf", height = 10, width = 10, device = cairo_pdf)