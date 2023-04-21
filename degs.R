#' Create the DEGs from a comparisons list.
#' 
#' @param res The comparisons list created with 'create_comparisons'
#' @return The DEGs
create_degs <- function(res) {
  degs = list()
  
  for (res_name in names(res)) {
    deg_name = str_replace(res_name, 'res', 'DEG')
    deg_p_name = str_replace(res_name, 'res', 'DEG_p')
    
    degs[[deg_name]] <- res[[res_name]] %>% select_degs()
    degs[[deg_p_name]] <- res[[res_name]] %>% select_degs_only_padj()
  }
  
  return(degs)
}


#' Create a list of DEGs up or down regulated.
#' 
#' @param degs The degs created using the 'create_degs' function
#' @return The DEGs separated as up or down regulated
create_up_down_degs <- function(degs) {
  up_down_degs = list()
  
  for (deg_name in names(degs)) {
    up_deg_name = paste(deg_name, "UP", sep = '.')
    down_deg_name = paste(deg_name, "DOWN", sep = '.')
    
    up_down_degs[[up_deg_name]] = subset(degs[[deg_name]], log2FoldChange > 1)
    up_down_degs[[down_deg_name]] = subset(degs[[deg_name]], log2FoldChange < -1)
  }
  
  return(up_down_degs)
}


#' Filter the data by adjusted p-value.
#' 
#' @param .data The data to filter
#' @return The filtered data
select_degs_only_padj <- function(.data) {
  to_tibble(.data) %>%
    filter(padj < 0.05)
}


#' Filter the data by adjusted p-value and log2 fold change.
#' 
#' @param .data The data to filter
#' @retun The filtered data
select_degs <- function(.data) {
  to_tibble(.data) %>%
    filter(padj < 0.05) %>%
    filter(abs(log2FoldChange) > 1)
}