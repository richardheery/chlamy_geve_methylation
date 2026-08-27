# Standardize Crein2937_GEVE.merged.gff and export as a GFF3

# Load required packages
library(rtracklayer)

# Load GFF file
chlamy_gff = import.gff("Crein2937_GEVE.merged.gff")

# Convert group to a character
chlamy_gff$group = as.character(chlamy_gff$group)

# Change mRNA type to transcript 
chlamy_gff[chlamy_gff$type == "mRNA"]$type = "transcript"

# Create exons from CDS, 5' UTR and 3' UTR regions
chlamy_exons = chlamy_gff[chlamy_gff$type %in% c("CDS", "5'-UTR", "3'-UTR")]
chlamy_exons$type = "exon"
chlamy_gff = c(chlamy_gff, chlamy_exons)

# Separate GEVE and host sequences. Some GEVE genes have no transcripts
host_gff = chlamy_gff[grep("GEVE", chlamy_gff$group, invert = T)]
geve_gff = chlamy_gff[grep("GEVE", chlamy_gff$group)]

# Separate host gff into genes, transcripts, exons and CDS
host_genes = host_gff[host_gff$type == "gene"]
host_transcripts = host_gff[host_gff$type == "transcript"]
host_exons = host_gff[host_gff$type == "exon"]
host_cds = host_gff[host_gff$type == "CDS"]

# For genes, set ID as group and Parent as NA and add gene_id from ID
host_genes$ID = host_genes$group
host_genes$Parent = NA
host_genes$gene_id = host_genes$ID

# For transcripts, use group as ID, set Parent as gene and add gene_id
host_transcripts$ID = host_transcripts$group
host_transcripts$Parent = gsub("\\..*", "", host_transcripts$group)
host_transcripts$gene_id = host_transcripts$Parent

# For exons extract transcript from group and set it as Parent
host_exons$Parent = stringr::str_match(host_exons$group, "\"(.*?)\"")[, 2]

# Sort exons and append exon number to transcript name to make ID
host_exons = host_exons[order(seqnames(host_exons), strand(host_exons), 
  ifelse(strand(host_exons) == "+", start(host_exons), -start(host_exons)))]
host_exons$ID = paste0(host_exons$Parent, "_exon_", ave(seq_along(host_exons$Parent), host_exons$Parent, FUN = seq_along))

# Add gene ID to host exons
host_exons$gene_id = gsub("\\..*", "", host_exons$Parent)

# For CDS extract transcript from group and set it as Parent
host_cds$Parent = stringr::str_match(host_cds$group, "\"(.*?)\"")[, 2]

# Sort CDS and append CDS number to transcript name to make ID
host_cds = host_cds[order(seqnames(host_cds), strand(host_cds), 
  ifelse(strand(host_cds) == "+", start(host_cds), -start(host_cds)))]
host_cds$ID = paste0(host_cds$Parent, "_cds_", ave(seq_along(host_cds$Parent), host_cds$Parent, FUN = seq_along))

# Add gene ID to host exons
host_cds$gene_id = gsub("\\..*", "", host_cds$Parent)

# Separate geve gff into genes, transcripts and exons
geve_genes = geve_gff[geve_gff$type == "gene"]
geve_transcripts = geve_gff[geve_gff$type == "transcript"]
geve_exons = geve_gff[geve_gff$type == "exon"]
geve_cds = geve_gff[geve_gff$type == "CDS"]

# For genes, extract gene from group and set Parent to NA and set ID as gene ID
geve_genes$ID = stringr::str_match(geve_genes$group, "ID=(.*?)[;]")[, 2]
geve_genes$Parent = NA
geve_genes$gene_id = geve_genes$ID

# For transcripts, set Parent as transcript and use group as ID
geve_transcripts$ID = stringr::str_match(geve_transcripts$group, "ID=(.*?)[;]")[, 2]
geve_transcripts$Parent = stringr::str_match(geve_transcripts$group, "Parent=(.*?)[;]")[, 2]
geve_transcripts$gene_id = geve_transcripts$Parent

# The transcript with Parent GEVEv2_relic1_3|GEVEv2_19 is missing transcript number from its ID so .t1 is added to its ID
geve_transcripts[geve_transcripts$Parent == "GEVEv2_relic1_3|GEVEv2_19"]$ID = 
  paste0(geve_transcripts[geve_transcripts$Parent == "GEVEv2_relic1_3|GEVEv2_19"]$ID, ".t1")

# Some GEVE genes are missing transcripts so the gene sequences are used as transcripts
no_transcript_geve_genes = geve_genes[!geve_genes$ID %in% geve_transcripts$Parent]
no_transcript_geve_transcripts = no_transcript_geve_genes
no_transcript_geve_transcripts$type = "transcript"
no_transcript_geve_transcripts$Parent = no_transcript_geve_transcripts$ID
no_transcript_geve_transcripts$ID = paste0(no_transcript_geve_transcripts$ID, ".t1")
geve_transcripts = c(geve_transcripts, no_transcript_geve_transcripts)

# For exons extract exon from group and set transcript as Parent
geve_exons$ID = stringr::str_match(geve_exons$group, "ID=(.*?)[;]")[, 2]
geve_exons$Parent = stringr::str_match(geve_exons$group, "Parent=(.*?)[;]")[, 2]

# Add .t1 to parents of exons from genes originally lacking transcripts
geve_exons[geve_exons$Parent %in% no_transcript_geve_genes$ID]$Parent = 
  paste0(geve_exons[geve_exons$Parent %in% no_transcript_geve_genes$ID]$Parent, ".t1")

# Set ID by adding _1 to parents of exons from genes originally lacking transcripts
geve_exons[geve_exons$Parent %in% no_transcript_geve_transcripts$ID]$ID = 
  paste0(geve_exons[geve_exons$Parent %in% no_transcript_geve_transcripts$ID]$Parent, "_exon_1")

# Change CDS in IDs to exon
geve_exons$ID = gsub("CDS", "exon", geve_exons$ID)

# For CDS extract exon from group and set it as Parent
geve_cds$ID = stringr::str_match(geve_cds$group, "ID=(.*?)[;]")[, 2]
geve_cds$Parent = stringr::str_match(geve_cds$group, "Parent=(.*?)[;]")[, 2]

# Add .t1 to parents of cds from genes originally lacking transcripts
geve_cds[geve_cds$Parent %in% no_transcript_geve_genes$ID]$Parent = 
  paste0(geve_cds[geve_cds$Parent %in% no_transcript_geve_genes$ID]$Parent, ".t1")

# Set ID by adding _1 to parents of cds from genes originally lacking transcripts
geve_cds[geve_cds$Parent %in% no_transcript_geve_transcripts$ID]$ID = 
  paste0(geve_cds[geve_cds$Parent %in% no_transcript_geve_transcripts$ID]$Parent, "_cds_1")

# Change CDS in IDs to exon
geve_cds$ID = gsub("CDS", "cds", geve_cds$ID)

# Combine genes, transcript and exons for host and GEVE
updated_chlamy_annotation_gr = c(host_genes, host_transcripts, host_exons, host_cds, 
  geve_genes, geve_transcripts, geve_exons, geve_cds)
mcols(updated_chlamy_annotation_gr) = mcols(updated_chlamy_annotation_gr)[c("source", "type", "score", "phase", "ID", "Parent")]

# Sort ranges and export as a GFF3
updated_chlamy_annotation_gr = sort(updated_chlamy_annotation_gr, ignore.strand = T)
rtracklayer::export.gff3(updated_chlamy_annotation_gr, "updated_chlamy_annotation.gff3")