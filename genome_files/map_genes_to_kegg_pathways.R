# Create a list linking T2T gene IDs to chlamydomoas-specific KEGG pathways

# Load required packages
library(dplyr)
library(rtracklayer)
library(Biostrings)

# Load GTF for Chalmy 5.5 and filter for CDS 
chlamy_5.5_gtf = import.gff2("GCF_000002595.2_Chlamydomonas_reinhardtii_v5.5_genomic.gtf.gz")
chlamy_5.5_cds = chlamy_5.5_gtf[chlamy_5.5_gtf$type == "CDS"]

# Make a vector connecting gene and protein IDs
gene_to_protein_df = data.frame(protein_id = chlamy_5.5_cds$protein_id, gene_id = chlamy_5.5_cds$locus_tag)
gene_to_protein_df = filter(gene_to_protein_df, !duplicated(protein_id))
gene_to_protein_vec = setNames(gene_to_protein_df$gene_id, gene_to_protein_df$protein_id)

# Load amino acid sequneces for proteins and rename using gene names
chlamy_v5.5_proteins = readAAStringSet("GCF_000002595.2_Chlamydomonas_reinhardtii_v5.5_protein.faa.gz")
names(chlamy_v5.5_proteins) = gsub(" .*", "", names(chlamy_v5.5_proteins))
names(chlamy_v5.5_proteins) = gene_to_protein_vec[names(chlamy_v5.5_proteins)]
writeXStringSet(chlamy_v5.5_proteins, "GCF_000002595.2_Chlamydomonas_reinhardtii_v5.5_protein.faa.gz", compress = TRUE)

# Load table matching Chlamy KEGG pathway IDs and names and shorten IDs and names
chlamy_kegg_pathway_name_to_id_df = data.table::fread("cre_pathway_names.txt", col.names = c("ID", "Name"), header = F)
chlamy_kegg_pathway_name_to_id_df$Name = gsub(" - Chlamydomonas reinhardtii", "", chlamy_kegg_pathway_name_to_id_df$Name)
chlamy_kegg_pathway_name_to_id_vec = setNames(chlamy_kegg_pathway_name_to_id_df$Name, chlamy_kegg_pathway_name_to_id_df$ID)

# Load table matching genes to pathway ID and shorten gene and pathway IDs
kegg_pathway_genes = data.table::fread("cre_genes_to_pathways.txt", col.names = c("gene_id", "pathway_id"), header = F)
kegg_pathway_genes$gene_id = gsub("cre:", "", kegg_pathway_genes$gene_id)
kegg_pathway_genes$pathway_id = gsub("path:", "", kegg_pathway_genes$pathway_id)

# Remove genes missing from gene_to_protein_vec
kegg_pathway_genes = filter(kegg_pathway_genes, gene_id %in% gene_to_protein_vec)

# Find most similar Cr v5.5 protein for each CC2937-T2T assembly protein using Diamond
system("diamond makedb --in GCF_000002595.2_Chlamydomonas_reinhardtii_v5.5_protein.faa.gz --db chlamy_v5.5_proteins")
system("diamond blastp --query host_proteins.fa --db chlamy_v5.5_proteins --out t2t_vs_v5.5_diamond.tsv --outfmt 6 qseqid sseqid qlen length qcovhsp pident evalue bitscore --max-target-seqs 1")

# Read in table of matches and filter for matches with E-value < 1e-10
t2t_vs_v5.5_diamond = data.table::fread("t2t_vs_v5.5_diamond.tsv", header = F, 
  col.names = c("t2t_gene", "chlamy_v5.5_gene", "qlen", "length", "qcovhsp", "pident", "evalue", "bitscore"))
t2t_vs_v5.5_diamond = filter(t2t_vs_v5.5_diamond, evalue < 1e-10)

# Make a list matching chlamy v5.5 genes to T2T genes
chlamy_v5.5_to_t2t_genes = setNames(t2t_vs_v5.5_diamond$t2t_gene, t2t_vs_v5.5_diamond$chlamy_v5.5_gene)

# Add T2T genes to kegg_pathway_genes and remove genes  with no matching T2T gene
kegg_pathway_genes$t2t_gene = chlamy_v5.5_to_t2t_genes[kegg_pathway_genes$gene_id]
kegg_pathway_genes = filter(kegg_pathway_genes, !is.na(t2t_gene))

# Create a list giving genes in each pathway 
kegg_pathway_genes_list = split(kegg_pathway_genes$t2t_gene, kegg_pathway_genes$pathway_id)
names(kegg_pathway_genes_list) = chlamy_kegg_pathway_name_to_id_vec[names(kegg_pathway_genes_list)]
saveRDS(kegg_pathway_genes_list, "kegg_pathway_genes_list.rds")
