#' Summarize values in bigWig files for a set of genomic regions
#' 
#' Calculates the mean (with or without non-covered bases counting as zeros) or the sum of values in bigWig files 
#' for provided genomic regions. Regions can be provided either as a GRanges object or as a path to a BED file.
#' Columns in the output file correspond to input bigWig files and rows correspond to genomic regions. 
#' Uses names of the provided BED file or GRanges as row names for the output if they are unique. 
#' Otherwise names regions region_1, region_2, etc. 
#'
#' @param bw_filepaths A vector of filepaths for bigWig files.
#' @param bed_filepath Filepath for a BED file with genomic regions of interest. One of either bed_filepath or regions_granges must be provided. 
#' @param regions_granges A GRanges object. One of either bed_filepath or regions_granges must be provided. 
#' @param bigwig_average_over_bed_path Optional path to bigWigAverageOverBed executable. 
#' Will attempt to search for bigWigAverageOverBed on the system PATH if not provided.
#' @param statistic Statistic to compute. One of "mean" (average in bigWig files over just covered bases in regions), 
#' "mean0" (average over bases with non-covered bases counting as zeroes) or "sum".
#' @param column_names A vector of names to use as column names of the output table. 
#' Default is to use the base name of bw_filepaths.
#' @param parallel_files The number of bigWig files to process at a time in parallel. Default is just 1.
#' @return A data.frame with the results.
#' @export
bigwig_summarize_over_regions = function(bw_filepaths, bed_filepath = NULL, bigwig_average_over_bed_path = NULL,
  regions_granges = NULL, statistic, column_names = NULL, parallel_files = 1){
  
  # If bigwig_average_over_bed_path is NULL, check if it is on the path
  if(is.null(bigwig_average_over_bed_path)){
    bigwig_average_over_bed_path = Sys.which("bigWigAverageOverBed")
  }
  if(bigwig_average_over_bed_path == ""){
    stop("Could not find bigWigAverageOverBed on system PATH and bigwig_average_over_bed_path is not provided")
  }
  
  # Check that bigwig_average_over_bed_path actually exists
  if(!file.exists(bigwig_average_over_bed_path)){stop("Provided bigwig_average_over_bed_path does not exist: ", bigwig_average_over_bed_path)}
  
  # Check that one and only one of bed_filepath or regions_granges is provided
  if(!is.null(regions_granges) && !is.null(bed_filepath)){
    stop("Either bed_filepath or regions_granges should be provided, but not both")
  } else if(is.null(regions_granges) && is.null(bed_filepath)){
    stop("One of bed_filepath or regions_granges must be provided")
  }

  # If regions_granges is provided add names if they are missing or are not unique
  if(!is.null(regions_granges) && is.null(bed_filepath)){
    
    # Name regions if they are not unique
    if(is.null(names(regions_granges)) || anyDuplicated(names(regions_granges))){
      message("names(regions_granges) are not unique. Will thus name them region_1, region_2, etc.")
      names(regions_granges) = paste0("region_", seq_along(regions_granges))
    }
    
    # Remove regions_granges metadata and export regions_granges as a temporary BED file
    mcols(regions_granges) = NULL
    bed_filepath = tempfile(pattern = "temp_bed")
    on.exit(unlink(bed_filepath), add = TRUE)
    rtracklayer::export.bed(regions_granges, bed_filepath)
    
  }
  
  # Check that one of the correct choices is provided for statistic
  match.arg(arg = statistic, choices = c("mean", "mean0", "sum"), several.ok = F)
  
  # Check that if column_names is provided, it has the same length as bw_filepaths
  if(!is.null(column_names)){if(length(bw_filepaths) != length(column_names)){stop("Length of bw_filepaths not equal to length of column_names")}}
  
  # Get names of regions from BED file and check that they are unique
  regions_granges = rtracklayer::import.bed(bed_filepath)
  region_names = regions_granges$name
  if(anyDuplicated(region_names)){
    message("Names of regions in BED are not unique. Will thus name them region_1, region_2, etc.")
    region_names = paste0("region_", seq_along(region_names))
    
    # Update region_names and export regions_granges to a temporary BED file
    regions_granges$name = region_names
    bed_filepath = tempfile(pattern = "temp_bed")
    on.exit(unlink(bed_filepath), add = TRUE)
    rtracklayer::export.bed(regions_granges, bed_filepath)
  }
  
  # Select the mean, mean0 or sum  columns from the bigWigAverageOverBed output files for downstream use
  statistic_col = ifelse(statistic == "sum", 4, ifelse(statistic == "mean0", 5, 6))
  
  # Create temporary output directory
  temp_dir = tempfile("bigwig_summarize_over_regions_tmp")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)
  
  # Create names for the output files that will be produced by bigWigAverageOverBed  
  summary_over_bed_file_names = paste(temp_dir, paste0(tools::file_path_sans_ext(basename(bw_filepaths)), "_summary_over_bed.txt"), sep ="/")
  
  # If no sample names provided, set to the basename of the input bigWig files
  if(is.null(column_names)){column_names = basename(tools::file_path_sans_ext(bw_filepaths))}
  
  # Create cluster if parallel_files greater than 1
  if(parallel_files > 1){
    cluster = parallel::makeCluster(parallel_files)
    doParallel::registerDoParallel(cluster, cores = parallel_files)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    `%dopar%` = foreach::`%dopar%`
  } else {
    `%dopar%` = foreach::`%do%`
  }
  
  # Use bigWigAverageOverBed to calculate mean value for each bigwig file for each region in the BED file and save the results in output directory
  foreach::foreach(bigwig_number = seq_along(bw_filepaths)) %dopar% {
    system2(bigwig_average_over_bed_path, args = c(bw_filepaths[bigwig_number], bed_filepath, 
      summary_over_bed_file_names[bigwig_number]))}
  
  # Create a matrix with rows corresponding to genomic regions and columns corresponding to genomic features
  region_values = data.frame(foreach::foreach(result_file = 
   summary_over_bed_file_names, .combine = "cbind") %dopar% {
     data.frame(data.table::fread(result_file, sep = "\t", header = F, select = c(1, statistic_col), nThread = 1), row.names = 1)[region_names, ]
  })
  
  # Assign row names and column names to data.frame and return it
  rownames(region_values) = region_names
  colnames(region_values) = column_names
  return(region_values)
  
}