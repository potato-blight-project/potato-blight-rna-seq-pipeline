#' Convert a dataframe to a tibble
#' 
#' @param .data The input data
to_tibble <- function (.data) {
  data.frame(.data) %>%
    rownames_to_column(var='gene') %>%
    as_tibble()
}


#' Create a dataframe as the intersection of two vectors and label the column as 'gene'.
#' 
#' @param a First vector
#' @param b Second vector
#' @return The intersected dataframe
in_intersection <- function(a, b) {
  result <- as.data.frame(intersect(a, b))
  names(result) <- c('gene')
  
  return(result)
}
