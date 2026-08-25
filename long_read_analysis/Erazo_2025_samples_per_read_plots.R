# Plot per-read CpG methylation for reads with and without deletion 

# Load required packages and functions
library(Rsamtools)
library(dplyr)
source("../auxillary_scripts/filter_bam_for_regions.R")
source("../auxillary_scripts/modkit_functions.R")
source("../auxillary_scripts/ggplot_functions.R")
source("../auxillary_scripts/methylartist_locus_plot.R")

# Set path to chlamy genome FASTA
chlamy_fasta = "../genome_files/CC2937_T2T.fa"

# Import GRanges with GEVE and deletion
geve_gr = rtracklayer::import.bed("../genome_files/GEVE_T2T.bed")
deletion_gr = rtracklayer::import.bed("../genome_files/GEVE_deletion_T2T.bed")

# Get GEVE regions upstream and downstream of the deletion
upstream_region = setdiff(geve_gr, deletion_gr)[1]
downstream_region = setdiff(geve_gr, deletion_gr)[2]

# Filter downstream_region to region most affected by 5mC change
end(downstream_region) = 1600000

# Get paths to BAM files with reads overlapping the GEVEs
first_run_bam = "bam_files/FirstRun_GEVEs_5mCG_calls.sorted_remap_CC2937_T2T.bam"
second_run_bam = "bam_files/SecondRun_GEVEs_5mCG_calls.sorted_remap_CC2937_T2T.bam"
third_run_bam = "bam_files/ThirdRun_GEVEs_5mCG_calls.sorted_remap_CC2937_T2T.bam"

# Get the names of reads overlapping upstream region, deletion and downstream region in first, second and third runs
first_run_upstream_reads = filter_bam_for_regions(first_run_bam, upstream_region, nthreads = 10)[[1]]$qname
first_run_deletion_reads = filter_bam_for_regions(first_run_bam, deletion_gr, nthreads = 10)[[1]]$qname
first_run_downstream_reads = filter_bam_for_regions(first_run_bam, downstream_region, nthreads = 10)[[1]]$qname
second_run_upstream_reads = filter_bam_for_regions(second_run_bam, upstream_region, nthreads = 10)[[1]]$qname
second_run_deletion_reads = filter_bam_for_regions(second_run_bam, deletion_gr, nthreads = 10)[[1]]$qname
second_run_downstream_reads = filter_bam_for_regions(second_run_bam, downstream_region, nthreads = 10)[[1]]$qname
third_run_upstream_reads = filter_bam_for_regions(third_run_bam, upstream_region, nthreads = 10)[[1]]$qname
third_run_deletion_reads = filter_bam_for_regions(third_run_bam, deletion_gr, nthreads = 10)[[1]]$qname
third_run_downstream_reads = filter_bam_for_regions(third_run_bam, downstream_region, nthreads = 10)[[1]]$qname

# Find reads indicating the presence of the deletion by overlapping both upstream and downstream regions but not the deletion
first_run_reads_supporting_deletion = setdiff(intersect(first_run_upstream_reads, first_run_downstream_reads), first_run_deletion_reads)
second_run_reads_supporting_deletion = setdiff(intersect(second_run_upstream_reads, second_run_downstream_reads), second_run_deletion_reads)
third_run_reads_supporting_deletion = setdiff(intersect(third_run_upstream_reads, third_run_downstream_reads), third_run_deletion_reads)

# Find reads indicating the absence of the deletion by overlapping the deletion and either of the upstream regions
first_run_reads_not_supporting_deletion = c(intersect(first_run_upstream_reads, first_run_deletion_reads), 
  intersect(first_run_downstream_reads, first_run_deletion_reads))
second_run_reads_not_supporting_deletion = c(intersect(second_run_upstream_reads, second_run_deletion_reads), 
  intersect(second_run_downstream_reads, second_run_deletion_reads))
third_run_reads_not_supporting_deletion = c(intersect(third_run_upstream_reads, third_run_deletion_reads), 
  intersect(third_run_downstream_reads, third_run_deletion_reads))

# Define a function that will summarize modification calls separately for reads supporting and not supporting the deletion
summarize_per_read_mod_calls_by_deletion = function(bam, del_present_reads, del_absent_reads){
  
  # Extract modification calls for reads supporting deletion and filter for high-confidence modifications
  del_present_reads_mod_calls = extract_mod_calls_for_reads(modBAM = bam, 
  read_names = del_present_reads, nthreads = 10, reference_fasta = chlamy_fasta)
  del_present_reads_mod_calls = filter_modkit_calls(modkit_table = del_present_reads_mod_calls, 
    sequence_context = "CG", sequence_context_matches_ref = T)
  
  # Summarize per-read methylation calls and filter for reads with at least 10 calls
  del_present_reads_mod_summary = summarize_mod_calls_per_read(del_present_reads_mod_calls)
  del_present_reads_mod_summary = filter(del_present_reads_mod_summary, total_calls >= 10)
  
  
  # Extract modification calls for reads not supporting deletion and filter for high-confidence modifications
  del_absent_reads_mod_calls = extract_mod_calls_for_reads(modBAM = bam, 
  read_names = del_absent_reads, nthreads = 10, reference_fasta = chlamy_fasta)
  del_absent_reads_mod_calls = filter_modkit_calls(modkit_table = del_absent_reads_mod_calls, 
    sequence_context = "CG", sequence_context_matches_ref = T)
  
  # Summarize per-read methylation calls and filter for reads with at least 10 calls
  del_absent_reads_mod_summary = summarize_mod_calls_per_read(del_absent_reads_mod_calls)
  del_absent_reads_mod_summary = filter(del_absent_reads_mod_summary, total_calls >= 10)
  
  # Combine tables for deletion and non-deletion reads
  combined_summary_table = bind_rows(list(
    deletion = del_present_reads_mod_summary, 
    no_deletion = del_absent_reads_mod_summary), .id = "deletion_status")
  return(combined_summary_table)
  
}

# Summarize per-read modification calls for first, second and third runs based on deletion 
first_run_mod_calls_by_deletion_summary = summarize_per_read_mod_calls_by_deletion(bam = first_run_bam, 
  del_present_reads = first_run_reads_supporting_deletion, del_absent_reads = first_run_reads_not_supporting_deletion)
second_run_mod_calls_by_deletion_summary = summarize_per_read_mod_calls_by_deletion(bam = second_run_bam, 
  del_present_reads = second_run_reads_supporting_deletion, del_absent_reads = second_run_reads_not_supporting_deletion)
third_run_mod_calls_by_deletion_summary = summarize_per_read_mod_calls_by_deletion(bam = third_run_bam, 
  del_present_reads = third_run_reads_supporting_deletion, del_absent_reads = third_run_reads_not_supporting_deletion)

# Combine all tables into a single overall table
combined_mod_calls_by_deletion_summary = bind_rows(list(
  first_run = first_run_mod_calls_by_deletion_summary,
  second_run = second_run_mod_calls_by_deletion_summary,
  third_run = third_run_mod_calls_by_deletion_summary), .id = "run_number")
data.table::fwrite(combined_mod_calls_by_deletion_summary, sep = "\t", "combined_mod_calls_by_deletion_summary.tsv.gz")
combined_mod_calls_by_deletion_summary = data.table::fread("combined_mod_calls_by_deletion_summary.tsv.gz")

# Create a data.frame giving the numbers of reads for each group
combined_mod_calls_by_deletion_summary_count = summarise(group_by(combined_mod_calls_by_deletion_summary, deletion_status, run_number), 
  count = paste("n =", n()), y = max(m)+0.1)

# Create boxplots of per-read methylation averages
# About 60% of the reads in the 2nd run have the deletion while only 1% in the first and 3rd runs do
set.seed(123)
deletion_vs_no_deletion_boxplots = ggplot(combined_mod_calls_by_deletion_summary, aes(y = m*100, x = deletion_status, fill = deletion_status)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.375, alpha = 0.5, shape = 21) +
  geom_text(data = combined_mod_calls_by_deletion_summary_count, mapping = aes(y = y*100, label = count))
deletion_vs_no_deletion_boxplots = customize_ggplot_theme(deletion_vs_no_deletion_boxplots,
  ylab = "Per-Read mCG %", x_labels = c("Deletion", "No Deletion"), fill_colors = c("#b97978", "#e1d7AA"), facet_labels = c("First Run", "Second Run", "Third Run"), 
  facet = "run_number", show_legend = F, facet_scales = "fixed") + theme(strip.background = NULL)
deletion_vs_no_deletion_boxplots
ggsave(plot = deletion_vs_no_deletion_boxplots, "plots/deletion_vs_no_deletion_boxplots.pdf", width = 9, height = 9)

### Make methylartist plots

# Set paths to BAM files deletion tagged reads
first_run_deletion_tagged_bam = "bam_files/FirstRun_GEVEs_5mCG_calls_downstream_deletion_tagged.bam"
second_run_deletion_tagged_bam = "bam_files/SecondRun_GEVEs_5mCG_calls_downstream_deletion_tagged.bam"
third_run_deletion_tagged_bam = "bam_files/ThirdRun_GEVEs_5mCG_calls_downstream_deletion_tagged.bam"

# Define a region 1 kb upstream and downstream of the 3' end of the deletion
deletion_downstream_end = resize(deletion_gr, 1, fix = "end")
deletion_downstream_end_flanking = resize(deletion_downstream_end, width = 2001, fix = "center")

# Create Methylartist plots for reads
methylartist_locus_plot(modBAM = first_run_deletion_tagged_bam, region = deletion_downstream_end_flanking, 
  genome_fasta = chlamy_fasta, 
  additional_options = "--height 9 --width 9 --nticks 5 -p 0,10,1,0,0 --phased --readmarkersize 3 --hidelegend --samplepalette pink --primary_only", 
  outfile = "plots/first_run_methylartist_plot_downstream_deletion_tagged.svg")

methylartist_locus_plot(modBAM = second_run_deletion_tagged_bam, region = deletion_downstream_end_flanking, 
  genome_fasta = chlamy_fasta, 
  additional_options = "--height 9 --width 9 --nticks 5 -p 0,10,1,0,0 --phased --readmarkersize 3 --hidelegend --samplepalette pink --primary_only", 
  outfile = "plots/second_run_methylartist_plot_downstream_deletion_tagged.svg")

methylartist_locus_plot(modBAM = third_run_deletion_tagged_bam, region = deletion_downstream_end_flanking, 
  genome_fasta = chlamy_fasta, 
  additional_options = "--height 9 --width 9 --nticks 5 -p 0,10,1,0,0 --phased --readmarkersize 3 --hidelegend --samplepalette pink --primary_only", 
  outfile = "plots/third_run_methylartist_plot_downstream_deletion_tagged.svg")