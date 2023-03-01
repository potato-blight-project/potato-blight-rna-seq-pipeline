to_tibble <- function (.data) {
  data.frame(.data) %>%
    rownames_to_column(var="gene") %>%
    as_tibble()
}

in_intersection <- function(a, b) {
  result <- as.data.frame(intersect(a, b))
  names(result) <- c("gene")
  
  return(result)
}

merge_with_annotations <- function(x, annotations, to_remove, suffixes) {
  unique_annotations = lapply(annotations, function (x) { unique(x[,!(names(x) %in% to_remove)]) })
  
  return(merge(
    merge(x,
          unique_annotations[[1]],
          by = "gene",
          all.x = T
    ),
    unique_annotations[[2]],
    by = "gene",
    all.x = T,
    suffixes = suffixes
  ))
}