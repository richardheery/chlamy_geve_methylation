# Download FASTQ files for RNA-seq and align readsd

# Load required packages
library(parallel)

# Load required packages
source("../auxillary_scripts/sratools.R")

# Load SRA metadata
sra_metadata = data.table::fread("sra_metadata.txt")

# Get all run accessions for RNA-seq samples
rnaseq_accessions = sra_metadata$run_accession[1:10]
names(rnaseq_accessions) = sra_metadata$library_name[1:10]

# Prefetch files
system.time({sra_prefetch(srr_accessions = rnaseq_accessions, output_directory = "sra_files", parallel_files = 10)})

# Extract FASTQ files
srr_directories = list.files("sra_files", full.names = T)
system.time({sra_fastq_dump(srr_directory_paths = srr_directories, output_directory = "fastq_files", 
  compress_fastq_files = T, parallel_files = 10)})

### Align FASTQ files

# Get paths to FASTQ files
fastq_files = list.files("fastq_files", full.names = T)

# Get names of all forward FASTQ files
forward_fastqs = grep("_fwd", fastq_files, value = T)
reverse_fastqs = grep("_rev", fastq_files, value = T)

# Create names for BAM files
dir.create("bam_files")
bam_prefix_names = paste0("bam_files/", paste0(gsub("_paired.*", "", basename(forward_fastqs))))

# Set path to STAR index
index_path = "../genome_files/CC2937_T2T_star_index/"

# Define a function to align pairs of reads with STAR
star_align_pair = function(pair){
  
  print(paste("Starting pair", pair))
  alignment_command = paste("STAR", "--genomeDir", index_path, 
    "--readFilesIn", forward_fastqs[pair], reverse_fastqs[pair], 
    "--readFilesCommand zcat",
    "--outFilterMultimapNmax", "100",
    "--winAnchorMultimapNmax", "100",
    "--outSAMtype BAM SortedByCoordinate",
    "--runThreadN", "10",
    "--outFileNamePrefix", bam_prefix_names[pair]
  )
  system(alignment_command)
}

# Align pairs in parallel
mclapply(seq_along(forward_fastqs), star_align_pair, mc.cores = 4)