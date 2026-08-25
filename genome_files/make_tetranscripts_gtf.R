# Process CC2937_T2T.TEs_v3_8.bed to create a GTF file for use with TEtranscripts

# Load required packages
library(rtracklayer)

# Import TE annotation as a data.frame and set names of columns
chlamy_tes = data.frame(data.table::fread("CC2937_T2T.TEs_v3_8.bed"))
names(chlamy_tes) = c("seqnames", "start", "end", "name", "score", "strand", "classification", "divergence", "percent_of_consensus", "autonomy", "location")

# Extract class and family names from classification
chlamy_tes$class_id = gsub("/.*", "", chlamy_tes$classification)
chlamy_tes$family_id = gsub(".*/", "", chlamy_tes$classification)

# Convert to a GRanges
chlamy_tes_gr = makeGRangesFromDataFrame(chlamy_tes, starts.in.df.are.0based = T, keep.extra.columns = T)

# Set type as exon for all TEs and set source and phase as .
chlamy_tes_gr$type = "exon"
chlamy_tes_gr$source = "."
chlamy_tes_gr$phase = "."

# Add class, family, gene and transcript metadata columns
chlamy_tes_gr$gene_id = gsub("_.*", "", chlamy_tes_gr$name)
chlamy_tes_gr$gene_name = paste0(chlamy_tes_gr$gene_id, ":TE")
chlamy_tes_gr$transcript_id = chlamy_tes_gr$name

# Make transcript_id unique
chlamy_tes_gr$transcript_id = make.unique(chlamy_tes_gr$transcript_id)

# Select the necessary metadata columns and put them in correct order and export a GTF file suitable for use with the TEToolkit suite
mcols(chlamy_tes_gr) = mcols(chlamy_tes_gr)[c("source", "type", "score", "phase", "gene_id", "transcript_id", "family_id", "class_id", "gene_name")]
export.gff2(chlamy_tes_gr, "CC2937_T2T.TEs_v3_8.bed.gtf")