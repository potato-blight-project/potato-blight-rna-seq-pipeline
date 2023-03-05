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

volcano_plots_from_comparisons <- function(res, ncol, nrow, widths, heights, filename, condition = NULL) {
  volcano_plots = list()
  
  for (res_name in names(res)) {
    if (is.null(condition) || condition(res_name)) {
      volcano_plots[[res_name]] <- volcano_plot(res[[res_name]], res_name)
    }
  }
  
  g <- grid.arrange(grobs = volcano_plots,
               ncol = ncol,
               nrow = nrow,
               widths = widths,
               heights = heights)
  
  ggsave(filename, g)
}

get_venn_diagram <- function(x, title, categories, filename = NULL) {
  if (length(categories) == 2) {
    fill <- c("#EFC000FF","#CD534CFF")
    cat.dist <- c(0.055, 0.055)
  } else {
    fill <- c("#E69F00","#56B4E9", "#009E73")
    cat.dist <- c(0.055, 0.055, 0.055)
  }
  
  venn.diagram(filename = filename,
               x = x, 
               scaled = FALSE,
               col = "black",
               fill = fill, 
               euler.d = FALSE,
               main = title,
               category.names = categories,
               main.cex = 2,
               imagetype = 'png',
               lty = 1.5,
               cex = 2,
               lwd = 1,
               cat.cex = 1.5,
               cat.default.pos = "outer",
               # cat.pos = c(120, 120),
               cat.dist = cat.dist)
}

create_venn_diagrams_grid <- function(up_down_degs, ncol, nrow, widths, heights, filename, config) {
  venn_diagrams <- list()
  i = 1
  
  for (config_item in config) {
    deg_names = config_item[["x"]]
    title = config_item[["title"]]
    categories = config_item[["categories"]]
    
    genes <- lapply(deg_names, function(x) { up_down_degs[[x]]$gene })
    names(genes) <- lapply(deg_names, function (x) {
      a = str_split(x, "\\.")
      paste(a[[1]][2], a[[1]][3], sep = "_")
    })
    
    venn_diagrams[[i]] <- get_venn_diagram(genes, title, categories)
    
    i <- i + 1
  }
  
  g <- grid.arrange(grobs = lapply(venn_diagrams, grobTree),
               ncol = ncol,
               nrow = nrow,
               widths = widths,
               heights = heights)
  
  ggsave(filename, g)
}
