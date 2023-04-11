create_comparisons <- function(deseq_object, comparisons) {
  res = list()
  
  for (comparison_name in names(comparisons)) {
    res[[comparison_name]] = results(deseq_object, alpha = 0.05, contrast = list(comparisons[[comparison_name]]))
  }
  
  return(res)
}

exclude_variety <- function(hold_out_variety, ...) {
  result <- setdiff(intersect(...), hold_out_variety)
  result <- as.data.frame(result)
  colnames(result)[1] <- "gene"
  
  return(result)
}

only_in_tolerant <- function(degs, up_down_degs, varieties, treatment, only_p_value_filter=FALSE) {
  # indicesToRemove <- c(1, 3, 5, 6, 7)
  # c(2,4,5,6,9) -> potv4
  
  if (only_p_value_filter) {
    prefix = "DEG_p"
  } else {
    prefix = "DEG"
  }
  
  indices_to_remove <- c("transcript_id", "baseMean", "lfcSE", "stat", "pvalue")
  merge_suffixes <- c('_SMira', '_SShona')
  
  variety_names <- c("Duke", "SMira", "SShona")
  tolerant_variety_names <- c("SMira", "SShona")
 
  genes_up = lapply(variety_names, function (x) { up_down_degs[[paste(prefix, x, treatment, "UP", sep = ".")]]$gene })
  only_tol_up_genes <- do.call(exclude_variety, genes_up)
  
  genes_down = lapply(variety_names, function (x) { up_down_degs[[paste(prefix, x, treatment, "DOWN", sep = ".")]]$gene })
  only_tol_down_genes <- do.call(exclude_variety, genes_down)
  
  tolerant_varieties_annotations <- lapply(tolerant_variety_names, function (x) {
    name <- paste(prefix, x, treatment, sep = ".")
    degs[[name]]
  })
  
  only_tol = list(
    up = merge_with_annotations(
      only_tol_up_genes,
      tolerant_varieties_annotations,
      indices_to_remove,
      merge_suffixes),
    
    down = merge_with_annotations(
      only_tol_down_genes,
      tolerant_varieties_annotations,
      indices_to_remove,
      merge_suffixes)
  )
  
  return(only_tol)
}
