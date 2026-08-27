# Make RangedSummarizedExperiments with mCG and mAT values

# Load required packages
library(methodical)
library(dplyr)

# Read in CC2937 T2T as a DNAStringSet
chlamy_genome_CC2937 = Biostrings::readDNAStringSet("../genome_files/CC2937_T2T.fa")[1:17]

### Make RSE for 5mC

# Get meth sites for CC2937 genome
CC2937_cpgs_gr = extractMethSitesFromGenome(chlamy_genome_CC2937, standard_seqs_only = F)

# Get paths to bedMethyl files for 5mC data remapped to CC2937_T2T.fa
bedmethyl_files_5mC = list.files("t2t_bedmethyl_files", pattern = "5mC", full.names = T)

# Create colData for samples
coldata_5mc = data.frame(
  status = c("Full GEVE", "Deletion", "Deletion", "Full GEVE"),
  row.names = gsub("_5mCG.*", "", basename(bedmethyl_files_5mC))
)

# Create an RSE for bedmethyl_files_5mC
chlamy_5mc_rse = makeMethRSEFromInputFiles(meth_files = bedmethyl_files_5mC, 
  seqnames_col = 1, start_col = 2, total_reads_col = 10, meth_reads_col = 12, 
  zero_based = T, collapse_strands = T, sequence_context = "CG", meth_sites = CC2937_cpgs_gr, 
  sample_metadata = coldata_5mc, hdf5_dir = "CC2937_T2T_5mC_rse")

### Make RSE for 6mA

# Get meth sites for CC2937 genome
CC2937_apts_gr = extractMethSitesFromGenome(chlamy_genome_CC2937, pattern = "AT", standard_seqs_only = F)

# Get paths to bedMethyl files for 6mA data remapped to CC2937_T2T.fa
bedmethyl_files_6mA = list.files("t2t_bedmethyl_files", pattern = "6mA", full.names = T)

# Create colData for samples
coldata_6ma = data.frame(
  status = c("Full GEVE", "Deletion", "Deletion", "Full GEVE"),
  row.names = gsub("_6mA.*", "", basename(bedmethyl_files_6mA))
)

# Create an RSE for bedmethyl_files_6mA
chlamy_6ma_rse = makeMethRSEFromInputFiles(meth_files = bedmethyl_files_6mA, 
  seqnames_col = 1, start_col = 2, total_reads_col = 10, meth_reads_col = 12, 
  zero_based = T, collapse_strands = T, sequence_context = "AT", meth_sites = CC2937_apts_gr, 
  sample_metadata = coldata_6ma, hdf5_dir = "CC2937_T2T_6mA_rse")

### Export bedGraphs for 5mC and 6mA

# Mask methylation of CpGs covered by < 10 reads or >= 350 reads
chlamy_5mc_rse = HDF5Array::loadHDF5SummarizedExperiment("CC2937_T2T_5mC_rse")
colnames(chlamy_5mc_rse) = paste0(colnames(chlamy_5mc_rse), "_5mC")
assay(chlamy_5mc_rse, 1)[assay(chlamy_5mc_rse, 2) < 10 | assay(chlamy_5mc_rse, 2) >= 350] = NA

# Export bigWigs for each 5mC sample
export_files_from_rse(meth_rse = chlamy_5mc_rse, filetype = "bigWig", output_dir = "bigWigs")

# Mask methylation of CpGs covered by < 10 reads or >= 350 reads
chlamy_6ma_rse = HDF5Array::loadHDF5SummarizedExperiment("CC2937_T2T_6mA_rse")
colnames(chlamy_6ma_rse) = paste0(colnames(chlamy_6ma_rse), "_6mA")
assay(chlamy_6ma_rse, 1)[assay(chlamy_6ma_rse, 2) < 10 | assay(chlamy_6ma_rse, 2) >= 350] = NA

# Export bigWigs for each 6ma sample
export_files_from_rse(meth_rse = chlamy_6ma_rse, filetype = "bigWig", output_dir = "bigWigs")

### Create meth RSE for mCG and mAT for Erazo-Garcia et al. 2025 samples

# Read in Erazo 2025 assembly as a DNAStringSet
chlamy_genome_Erazo_2025 = Biostrings::readDNAStringSet("../genome_files/Erazo_2025.fasta")
names(chlamy_genome_Erazo_2025) = gsub(" .*", "", names(chlamy_genome_Erazo_2025))

# Get CpG sites for Erazo 2025
Erazo_2025_cpgs_gr = extractMethSitesFromGenome(chlamy_genome_Erazo_2025, standard_seqs_only = F)

# Get paths to bedMethyl files for 5mC for Erazo 2025
erazo_bedmethyl_files_5mC = list.files("erazo_bedmethyl_files", pattern = "5mCG", full.names = T)

# Create colData for samples
Erazo_2025_coldata = data.frame(row.names = c("FirstRun",  "SecondRun", "ThirdRun"))

# Create an RSE for erazo_bedmethyl_files_5mC
Erazo_2025_5mc_rse_contigs = makeMethRSEFromInputFiles(meth_files = erazo_bedmethyl_files_5mC, 
  seqnames_col = 1, start_col = 2, total_reads_col = 10, meth_reads_col = 12, 
  zero_based = T, collapse_strands = T, sequence_context = "CG", meth_sites = Erazo_2025_cpgs_gr, 
  sample_metadata = Erazo_2025_coldata, hdf5_dir = "Erazo_2025_5mC_rse_contigs")

# Get CpG sites for Erazo 2025
Erazo_2025_apts_gr = extractMethSitesFromGenome(chlamy_genome_Erazo_2025, pattern = "AT", standard_seqs_only = F)

# Get paths to bedMethyl files for 5mC for Erazo 2025
erazo_bedmethyl_files_6mA = list.files("erazo_bedmethyl_files", pattern = "6mA", full.names = T)

# Create an RSE for erazo_bedmethyl_files_5mC
Erazo_2025_6mA_rse_contigs = makeMethRSEFromInputFiles(meth_files = erazo_bedmethyl_files_6mA, 
  seqnames_col = 1, start_col = 2, total_reads_col = 10, meth_reads_col = 12, 
  zero_based = T, collapse_strands = T, sequence_context = "AT", meth_sites = Erazo_2025_apts_gr, 
  sample_metadata = Erazo_2025_coldata[2:3, ], hdf5_dir = "Erazo_2025_6mA_rse_contigs")

# Liftover Erazo RSEs to CC9937-T2T
chain = rtracklayer::import.chain("../genome_files/Erazo_2025_to_T2T.chain", exclude = NA)
Erazo_2025_5mC_rse_T2T = liftoverMethRSE(meth_rse = Erazo_2025_5mc_rse_contigs, chain = chain)
Erazo_2025_6mA_rse_T2T = liftoverMethRSE(meth_rse = Erazo_2025_6mA_rse_contigs, chain = chain)
HDF5Array::saveHDF5SummarizedExperiment(Erazo_2025_5mC_rse_T2T, "Erazo_2025_5mC_rse")
HDF5Array::saveHDF5SummarizedExperiment(Erazo_2025_6mA_rse_T2T, "Erazo_2025_6mA_rse")