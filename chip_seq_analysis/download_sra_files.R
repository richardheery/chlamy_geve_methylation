# Get data for selected histone marks in Chlamydomonas from papers:
# "Lineage-specific chromatin signatures reveal a regulator of lipid metabolism in microalgae" and
# "DNA cytosine methylation suppresses meiotic recombination at the sex-determining region"

# Import SRA Toolkit functions
source("../helper_scripts/sratools.R")

### Download files for PRJNA255778 (Lineage-specific chromatin signatures reveal a regulator of lipid metabolism in microalgae)

# Read table for SRA metadata for PRJNA255778
sra_metadata_PRJNA255778 = data.table::fread("SraRunTable_PRJNA255778.csv")

# Shorten source_name to just the name of the histone mark
sra_metadata_PRJNA255778$histone_mark = gsub("-", "_", gsub(" ", "_", gsub(" \\(.*", "", sra_metadata_PRJNA255778$source_name)))

# Create a vector matching SRR accession to histone mark
sra_to_name_PRJNA255778 = setNames(sra_metadata_PRJNA255778$histone_mark, sra_metadata_PRJNA255778$Run)

# Prefetch files with SRA Toolkit
sra_prefetch(srr_accessions = sra_metadata_PRJNA255778$Run, output_directory = "sra_files", parallel_files = 8)

# Get paths to downloaded directories
sra_directories = list.files("sra_files", full.names = T)

# Extract FASTQ files from directories
system.time({sra_fastq_dump(srr_directory_paths = sra_directories, output_directory = "fastq_files/PRJNA255778", parallel_files = 8)})

# Remove SRA directories
unlink(sra_directories, recursive = T)

# Get paths to FASTQ files and rename them to include histone mark in their names
fastq_files = list.files("fastq_files/PRJNA255778", full.names = T, pattern = "fastq.gz")
names(fastq_files) = gsub(".fastq.gz", "", basename(fastq_files))
file.rename(fastq_files, paste0("fastq_files/PRJNA255778/", sra_to_name_PRJNA255778[names(fastq_files)], "_", names(fastq_files), ".fastq.gz"))

### Download files for PRJNA1029152 (DNA cytosine methylation suppresses meiotic recombination at the sex-determining region)

# Read table for SRA metadata for PRJNA1029152
sra_metadata_PRJNA1029152 = data.table::fread("SraRunTable_PRJNA1029152.csv")

# Set histone_mark as chip_antibody and if antibody is none, set it as input
sra_metadata_PRJNA1029152$histone_mark = sra_metadata_PRJNA1029152$chip_antibody
sra_metadata_PRJNA1029152$histone_mark[sra_metadata_PRJNA1029152$histone_mark == "none"] = "input"

# Shorten histone_mark to name of chip antibody and append SRR accession
sra_metadata_PRJNA1029152$histone_mark = paste(gsub(" .*", "", sra_metadata_PRJNA1029152$histone_mark), sra_metadata_PRJNA1029152$Run, sep = "_")

# Create a vector matching SRR accession to histone mark
sra_to_name_PRJNA1029152 = setNames(sra_metadata$chip_antibody, sra_metadata$Run)

# Prefetch files with SRA Toolkit
sra_prefetch(srr_accessions = sra_metadata$Run, output_directory = "sra_files", parallel_files = 4)

# Get paths to downloaded directories
sra_directories = list.files("sra_files", full.names = T)

# Extract FASTQ files from directories
system.time({sra_fastq_dump(srr_directory_paths = sra_directories[1], output_directory = "fastq_files/PRJNA1029152", parallel_files = 1)})

# Remove SRA directories
unlink(sra_directories, recursive = T)

# Get paths to FASTQ files and rename them to include histone mark in their names
fastq_files = list.files("fastq_files/PRJNA1029152", full.names = T, pattern = "SRR264*.fastq.gz")
names(fastq_files) = gsub(".fastq.gz", "", basename(fastq_files))
file.rename(fastq_files, paste0("fastq_files/PRJNA1029152/", sra_to_name_PRJNA1029152[gsub("_.", "", names(fastq_files))], gsub(".*_", "_", names(fastq_files)), ".fastq.gz"))