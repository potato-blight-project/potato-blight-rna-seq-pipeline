create_degs <- function(res) {
  degs = list()
  
  for (res_name in names(res)) {
    deg_name = str_replace(res_name, "res", "DEG")
    deg_p_name = str_replace(res_name, "res", "DEG_p")
    
    degs[[deg_name]] <- res[[res_name]] %>% select_degs()
    degs[[deg_p_name]] <- res[[res_name]] %>% select_degs_only_padj()
  }
  
  return(degs)
}

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

select_degs_only_padj <- function(.data) {
  to_tibble(.data) %>%
    filter(padj < 0.05)
}

select_degs <- function(.data) {
  to_tibble(.data) %>%
    filter(padj < 0.05) %>%
    filter(abs(log2FoldChange) > 1)
}