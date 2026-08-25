# Export annotation files for different classes of genes

# Load required packages
library(dplyr)
library(rtracklayer)
library(GenomicFeatures)
library(BSgenome)

# Load genes from CC2937-T2T GTF and export
chlamy_txdb = txdbmaker::makeTxDbFromGFF("updated_chlamy_annotation_CC2937_T2T_liftoff_polished.gtf")
chlamy_genes = genes(chlamy_txdb)
saveRDS(chlamy_genes, "chlamy_genes_gr.rds")

# Export the GEVE TET gene GEVEv2_173 as a BED file
export.bed(chlamy_genes["GEVEv2_173"], "TET_GEVEv2_173.bed")

# Filter for host, GEVE and relict genes and export as BED files
host_genes = chlamy_genes[grep("GEVE", names(chlamy_genes), invert = T)]
geve_genes = chlamy_genes[intersect(grep("GEVE", names(chlamy_genes)), grep("relic", names(chlamy_genes), invert = T))]
relict_genes = chlamy_genes[grep("relic", names(chlamy_genes))]
export.bed(geve_genes, "GEVE_all_genes.bed")

# Load annotation for genes. GEVEv2_204 is removed since it is missing from GTF file
geve_genes_annotation = xlsx::read.xlsx("../genome_files/GEVE_annotation.xlsx", sheetIndex = 1)[-1]
geve_genes_annotation = filter(geve_genes_annotation, gene != "GEVEv2_204")

# Indicate if GEVE genes are BVP, DVP or unassigned
geve_genes$cluster = "unclustered"
geve_genes[filter(geve_genes_annotation, cluster == "BVP")$gene]$cluster = "BVP"
geve_genes[filter(geve_genes_annotation, cluster == "DVP")$gene]$cluster = "DVP"

# Create a list with GEVE gene clusters and host genes
gene_class_list = list(
  Host = names(host_genes),
  BVP = names(geve_genes[geve_genes$cluster == "BVP"]),
  DVP = names(geve_genes[geve_genes$cluster == "DVP"]),
  GEVE_unclustered = names(geve_genes[geve_genes$cluster == "unclustered"]),
  Relict = names(relict_genes)
)
saveRDS(gene_class_list, "gene_class_list.rds")

# Export BVP and DVP genes classes as BED files
export.bed(chlamy_genes[gene_class_list$BVP], "GEVE_BVP_genes.bed")
export.bed(chlamy_genes[gene_class_list$DVP], "GEVE_DVP_genes.bed")

# Turn gene_class_list into a vector
gene_class_vec = unlist(with(reshape2::melt(gene_class_list), split(L1, value)))

### Export amino acid sequences for host, relict and GEVE proteins

# Get names of host, GEVE and relict genes
host_genes = names(host_genes)
geve_genes = names(geve_genes)
relict_genes = names(relict_genes)

# Get CDS and get sequences for CDS
cds_gr = cdsBy(chlamy_txdb, "gene")
cds_seqs = extractTranscriptSeqs(readDNAStringSet("CC2937_T2T.fa"), cds_gr)

# Remove sequences containing Ns and whose length is not a multiple of 3
cds_seqs = cds_seqs[grep("N", cds_seqs, invert = T)]
cds_seqs = cds_seqs[width(cds_seqs) %% 3 == 0]

# Translate to proteins and export as a FASTA
protein_seqs = translate(cds_seqs)

# Save FASTA files with host, G£VE and relict proteins
host_protein_seqs = protein_seqs[names(protein_seqs) %in% host_genes]
geve_protein_seqs = protein_seqs[names(protein_seqs) %in% geve_genes]
relict_protein_seqs = protein_seqs[names(protein_seqs) %in% relict_genes]
writeXStringSet(host_protein_seqs, "host_proteins.fa")
writeXStringSet(geve_protein_seqs, "geve_proteins.fa")
writeXStringSet(relict_protein_seqs, "relict_proteins.fa")

### Assign genes in the relicts to BVP or DVP classes based on sequence similarity

# Find most similar GEVE protein for each relict protein using blastp
system("makeblastdb -in geve_proteins.fa -dbtype prot -out geve_proteins_blastdb/geve_proteins")
system("blastp -query relict_proteins.fa -db geve_proteins_blastdb/geve_proteins -outfmt \"6 qseqid sseqid qlen length qcovs pident evalue\" -max_target_seqs 1 -out relict_vs_geve_blastp_hits.tsv")

# Load blastp results aligning relict proteins to GEVE proteins
relict_vs_geve_blastp_hits = data.table::fread("relict_vs_geve_blastp_hits.tsv", 
  col.names = c("relict_gene", "geve_gene", "qlen", "length", "qcovs", "pident", "evalue"))

# Add gene class of most similar GEVE protein
relict_vs_geve_blastp_hits$geve_gene_class = gene_class_vec[relict_vs_geve_blastp_hits$geve_gene]

# Get the longest hit for each relict gene and remove those with E values < 1e-10
relict_vs_geve_blastp_hits = filter(group_by(relict_vs_geve_blastp_hits, relict_gene), length == max(length))
relict_vs_geve_blastp_hits = filter(relict_vs_geve_blastp_hits , evalue < 1e-10, geve_gene_class != "GEVE_unclustered")

# Create a vector matching relict genes to BVP or DVP genes
relict_gene_class_list = split(relict_vs_geve_blastp_hits$relict_gene, relict_vs_geve_blastp_hits$geve_gene_class)
saveRDS(relict_gene_class_list, "relict_gene_class_list.rds")