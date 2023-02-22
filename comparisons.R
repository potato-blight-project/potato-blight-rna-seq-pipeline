only_in_a_variety_type <- function(variety_a, variety_b, variety_c) {
  result <- setdiff(intersect(variety_a, variety_b), variety_c)
  result <- as.data.frame(result)
  colnames(result)[1] <- "gene"
  
  return(result)
}

only_in_tolerant <- function(degs, up_down_degs, varieties) {
  indices_to_remove <- c("transcript_id", "baseMean", "lfcSE", "stat", "pvalue")
  merge_suffixes <- c('_SMira', '_SShona')
  
  only_tol_up_genes <- only_in_a_variety_type(DEG.SMira.Out.UP$gene, DEG.SShona.Out.UP$gene, DEG.Duke.Out.UP$gene)
  only_tol_down_genes <- only_in_a_variety_type(DEG.SMira.Out.DOWN$gene, DEG.SShona.Out.DOWN$gene, DEG.Duke.Out.DOWN$gene)
  
  OnlyTol.Out.DEG.up.annotation <- merge_with_annotations(
    only_tol_up_genes,
    DEG.SMira.Out.annotation,
    DEG.SShona.Out.annotation,
    indices_to_remove,
    merge_suffixes)
  
  OnlyTol.Out.DEG.down.annotation <- merge_with_annotations(
    only_tol_down_genes,
    DEG.SMira.Out.annotation,
    DEG.SShona.Out.annotation,
    indices_to_remove,
    merge_suffixes)
}
