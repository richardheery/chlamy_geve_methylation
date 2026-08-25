#' Test if elements of a query vector are enriched in a test vector using Fisher's exact test
#'
#' @param test A character vector. Should only contain unique elements. 
#' @param query A character vector. Should only contain unique elements. 
#' @param universe A vector containing all the elements that test and query could contain. Should only contain unique elements. 
#' @param alternative The alternative hypothesis for the test. Default is "greater".
#' @param return_overlap A logical value indicating whether to return the overlap between test and query as a character with elements separated by ";".
#' @return A data.frame with the results of Fisher's exact test
#' @export
fisher_test_vectors = function(test, query, universe, alternative = "greater", return_overlap = F){
  
  # Check if there are any duplicated elements in test, query or universe
  if(anyDuplicated(test)){stop("test cannot contain duplicated elements")}
  if(anyDuplicated(query)){stop("query cannot contain duplicated elements")}
  if(anyDuplicated(universe)){stop("universe cannot contain duplicated elements")}
  
  # Perform Fisher tests
  ft_result = fisher.test(
    x = factor(universe %in% test, levels = c("TRUE", "FALSE")), 
    y = factor(universe %in% query, levels = c("TRUE", "FALSE")), 
    alternative = alternative)
  
  # Calculate the size of the overlaps of the test and the universe with the query
  test_overlap_size = sum(test %in% query)
  universe_overlap_size = sum(universe %in% query)
  
  # Calcualte the proportion of the test vector which overlaps the query
  test_overlap_proportion = test_overlap_size/length(test)
  
  # Calculate the proportion of the universe which overlaps the query
  universe_overlap_proportion = universe_overlap_size/length(universe)
  
  # Calculate the enrichment of the query in the test relative to the universe
  relative_enrichment = test_overlap_proportion/universe_overlap_proportion
  
  # Create a data.frame with the results
  result_df = data.frame(test_size = length(test), query_size = length(query), universe_size = length(universe),
    test_overlap_size, universe_overlap_size, test_overlap_proportion, universe_overlap_proportion, relative_enrichment, p_value = ft_result$p.value)
  
  # Add the intersection of test and query if requested
  if(return_overlap){
    test_query_overlap = paste(intersect(test, query), collapse = "; ")
     if(length(test_query_overlap) == 0){test_query_overlap = NA}
    result_df$test_query_overlap = test_query_overlap
  }
  
  # Return results data.frame
  return(result_df)
}

#' Test if a list of query vectors is enriched in a test vector using Fisher's exact test
#'
#' @param test A character vector. Should only contain unique elements. 
#' @param query_list A list of character vectors to test for enrichment in test. Each vector should only contain unique elements. Names of list are included in the results.
#' @param universe A vector containing all the elements that test and query_list could contain. Should only contain unique elements. 
#' @param alternative The alternative hypothesis for the test. Default is "greater".
#' @param p_adjust_method Method to use to adjust p-values. Default is "fdr".
#' @param return_overlap A logical value indicating whether to return the overlap between test and query as a character with elements separated by ";". 
#' @return A data.frame with the resuls of Fisher's exact test
#' @export
fisher_test_apply = function(test, query_list, universe, alternative = "greater", p_adjust_method = "fdr", return_overlap = F, ncores = 1){
  
  # Create cluster if ncores greater than 1
  if(ncores > 1){
    cl = parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl, ncores)
    `%dopar%` = foreach::`%dopar%`
    on.exit(parallel::stopCluster(cl))
  } else {
    `%dopar%` = foreach::`%do%`
  }
  
  # Create a list with results for each vector in query_list
  results_list = foreach::foreach(query = query_list) %dopar% {
    
    fisher_test_vectors(test = test, query = query, universe = universe, 
      alternative = alternative, return_overlap = return_overlap)
    
  }
  
  # Name results_list with names of query_list
  names(results_list) = names(query_list)
    
  # Combine results_list into a single data.frame
  results_df = dplyr::bind_rows(results_list, .id = "query_name")
  results_df$q_value = p.adjust(results_df$p_value, method = p_adjust_method)
  
  # Swap q_value and test_query_overlap columns if return_overlap is set to TRUE
  if(return_overlap){
    results_df = dplyr::relocate(results_df, q_value, .before = test_query_overlap)
  }
  
  # Return results data.frame
  return(results_df)
}
