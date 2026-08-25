# Quantify TEs in BAM files with TEcount and TElocal

# Get paths to BAM files for full GEVE and deletion samples
bam_files = list.files("bam_files", pattern = ".bam$", full.names = T)
bam_sample_names = gsub("Aligned.sortedByCoord.out.bam", "", basename(bam_files))

# Set paths to gene and TE GTF files
genes_gtf = "../genome_files/updated_chlamy_annotation_CC2937_T2T_liftoff_polished.gtf"
te_index = "../genome_files/CC2937_T2T.TEs_v3_8.bed.gtf.locInd"

# Define a function call TEcount on BAM files
TEcount_quantify = function(bam_number){
  
  results_dir = paste0("TEcount_results/", bam_sample_names[bam_number])
  dir.create(results_dir)
  tecount_command = paste("TEcount", 
    "-b", bam_files[bam_number],
    "--format BAM",
    "--sortByPos",
    "--GTF", genes_gtf,
    "--TE", te_gtf,
    "--mode multi",
    "--outdir", results_dir
  )
  system(tecount_command)
}

# Run TEcount on BAM files
lapply(seq_along(bam_files), TEcount_quantify)
  
# Define a function call TElocal on BAM files
TElocal_quantify = function(bam_number){
  
  output = paste0("TElocal_results/", bam_sample_names[bam_number])
  telocal_command = paste("~/programs/TElocal/TElocal", 
    "-b", bam_files[bam_number],
    "--sortByPos",
    "--GTF", genes_gtf,
    "--TE", te_index,
    "--mode multi",
    "--project", output
  )
  system(telocal_command)
}

# Run TElocal on BAM files
lapply(seq_along(bam_files), TElocal_quantify)