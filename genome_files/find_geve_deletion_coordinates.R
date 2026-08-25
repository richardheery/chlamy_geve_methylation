# Find location of deletion in GEVE

# Load required packages
library(dplyr)
library(SummarizedExperiment)

# Read table with GEVE coverage
geve_coverage = data.table::fread("geve_depth.txt")
names(geve_coverage) = c("seqnames", "start", "chlamy_7A", "chlamy_8A", "chlamy_15A", "chlamy_18A")

# Add a column with genomic coordinates
geve_coverage = mutate(geve_coverage, coordinates = paste0(seqnames, ":", start), .keep = "unused")

# Convert to a RangedSummarizedExperiment
geve_coverage = SummarizedExperiment(assays = list(cov = dplyr::select(geve_coverage, starts_with("chlamy"))), 
  rowRanges = GRanges(geve_coverage$coordinates))

# Find uncovered regions (with coverage < 10) for each sample
uncovered_geve_regions = lapply(colnames(geve_coverage), function(x) 
  reduce(rowRanges(geve_coverage)[assay(geve_coverage)[[x]] <= 10]))
names(uncovered_geve_regions) = colnames(geve_coverage)

# The deletion is present at chromosome_15:1219979-1499485 in both 7A and 18A
deletion_gr = intersect(uncovered_geve_regions$chlamy_7A, uncovered_geve_regions$chlamy_18A)
export.bed(deletion_gr, "GEVE_deletion_T2T.bed")