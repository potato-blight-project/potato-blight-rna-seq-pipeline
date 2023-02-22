create_comparisons <- function(deseq_object, contrasts) {
  res = list()
  
  for (comparison_name in names(comparisons)) {
    res[[comparison_name]] = results(deseq_object, alpha = 0.05, contrast = list(comparisons[[comparison_name]]))
  }
  
  return(res)
}

volcano_plots_from_comparisons <- function(res, condition) {
  volcano_plots = list()
  
  for (res_name in names(res)) {
    if (condition(res_name)) {
      volcano_plots[[res_name]] <- volcano_plot(res[[res_name]], res_name)
    }
  }
  
  grid.arrange(grobs = volcano_plots,
               ncol = 2,
               nrow = 3,
               widths = unit(c(160,160), c("mm","mm")),
               heights = unit(c(150,150,150), c("mm","mm","mm")))
}

volcano_plot <- function(x, title) {
  EnhancedVolcano(x, 
                  lab = rownames(x),
                  x = 'log2FoldChange',
                  y = 'padj',
                  pCutoff = 0.05,
                  FCcutoff = 1,
                  xlim = c(-5, 8),
                  title = title)
}

venn_diagram_2way <- function(x, title) {
  venn.diagram(filename = NULL, x, 
               scaled = FALSE, col = "black", fill = c("#EFC000FF","#CD534CFF"), 
               euler.d = FALSE, lty = 1 , cex = 1, lwd = 1, cat.cex = 1, sigdigs = 1,
               print.mode=c("raw","percent"), height = 20, width = 20, main = title,
               rotation.degree = 180,
               cat.just=list(c(0.3,-1.5), c(0.7,-1.5)),cat.pos = c(-40,40),cat.dist = c(0.01, 0.01))
}

venn_diagram_3way <- function(x, title) {
  venn.diagram(filename = NULL, x, 
               scaled = FALSE, col = "black", fill = c("#E69F00","#56B4E9", "#009E73"),
               main = title,
               euler.d = FALSE,lty = 1 , cex = 0.9, lwd = 1, cat.cex = 1, sigdigs = 1,
               print.mode=c("raw","percent"), height = 20, width = 20, cat.pos = c(-20,20,-180), cat.dist = c(0.06, 0.06, 0.04))
}

to_tibble <- function (.data) {
  data.frame(.data) %>%
    rownames_to_column(var="gene") %>%
    as_tibble()
}

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

in_intersection <- function(a, b) {
  result <- as.data.frame(intersect(a, b))
  colnames(result)[1] <- "gene"
  
  return(result)
}

merge_with_annotations <- function(x, annotations_df_1, annotations_df_2, to_remove, suffixes) {
  unique_annotations_df_1 = unique(annotations_df_1[,!(names(annotations_df_1) %in% to_remove)])
  unique_annotations_df_2 = unique(annotations_df_2[,!(names(annotations_df_1) %in% to_remove)])
  
  return(merge(
    merge(x,
          unique_annotations_df_1,
          by = "gene",
          all.x = T
    ),
    unique_annotations_df_2,
    by = "gene",
    all.x = T,
    suffixes = suffixes
  ))
}