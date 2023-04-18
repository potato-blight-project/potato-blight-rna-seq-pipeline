filter_by_list <- function(
  annotated_degs,
  comparison,
  comparison_list,
  filter_by,
  use_only_func_ann = TRUE
) {
  annotated_degs[[comparison]] %>%
    filter(.data[[filter_by]] %in% comparison_list) %>%
    dplyr::select(-transcript_id, -pvalue, -stat, -baseMean, -TAIR, -GOterm) %>%
    dplyr::filter(if (use_only_func_ann) !is.na(name) else TRUE) %>%
    distinct(gene, .keep_all = TRUE)
}

create_block <- function(degs_list, suffixes) {
  block <- purrr::map2(degs_list, suffixes, ~{
    cols_to_rename <- setdiff(names(.x), 'gene')
    purrr::set_names(.x, c('gene', paste0(cols_to_rename, .y)))
  }) %>% purrr::reduce(merge, by = 'gene', all = TRUE)
}

create_heatmap_matrix <- function(
  block_list,
  suffixes,
  filter_fn,
  comparison_list,
  comparison_list_order,
  go_annotations
) {
  heatmap_data <- create_block(block_list, suffixes)
  
  heatmap_data_idx <- heatmap_data %>%
    dplyr::select(starts_with('log2FoldChange')) %>%
    mutate(id = row_number()) %>%
    filter(filter_fn(.)) %>%
    pull(id)
  
  heatmap_data <- heatmap_data %>% slice(heatmap_data_idx)
  
  func_names <- heatmap_data %>%
    dplyr::select(starts_with('name')) %>%
    rowwise() %>%
    mutate(func_name = list(first(na.omit(c_across(everything()))))) %>%
    dplyr::select(func_name)
  func_names <- unlist(lapply(func_names$func_name, function (x) { ifelse(length(x) == 0, NA, x) }))
  heatmap_data$label <- paste(heatmap_data$gene, func_names, sep = '_')
  heatmap_data <- heatmap_data %>% relocate(label, .after = gene)
  
  heatmap_data <- heatmap_data %>% dplyr::select(gene, label, starts_with('log2FoldChange'))
  
  if (!is.null(go_annotations)) {
    heatmap_data <- merge(heatmap_data, go_annotations, by="gene", all.x = T) %>%
      filter(GOterm %in% comparison_list) %>%
      dplyr::select(-transcript_id, -TAIR) %>%
      distinct(gene, .keep_all = TRUE) %>%
      relocate(GOterm, .after = label)
    
    if (!is.null(comparison_list_order)) {
      heatmap_data <- heatmap_data %>% arrange(match(GOterm, comparison_list_order))
    } else {
      heatmap_data <- heatmap_data %>% arrange(GOterm)
    }
    
    heatmap_data$label <- paste(heatmap_data$GOterm, heatmap_data$label, sep = '_')
  }
  
  return(list(
    heatmap_data = heatmap_data,
    heatmap_data_matrix = as.matrix(heatmap_data %>% dplyr::select(starts_with('log2FoldChange')))
  ))
}

plot_heatmap <- function(
  heatmap_definition,
  annotated_degs,
  comparison_list,
  filter_by,
  filter_fn,
  comparison_list_order = NULL,
  go_annotations = NULL,
  label_rows = TRUE,
  filename = 'heatmap.pdf',
  width = 8,
  height = 10
) {
  blocks <- list()
  
  for (block_name in names(heatmap_definition)) {
    block_def <- heatmap_definition[[block_name]]
    degs_list <- list()
    
    for (comparison in block_def$comparisons) {
      degs_list[[comparison]] <- filter_by_list(annotated_degs, comparison, comparison_list, filter_by)
    }
    
    blocks[[block_name]] <- create_block(degs_list, block_def$suffixes)
  }
  
  heatmap_data_result <- create_heatmap_matrix(
    block_list = blocks,
    suffixes = unlist(lapply(names(blocks), function(x) { paste0('_', x) })),
    filter_fn = filter_fn,
    comparison_list = comparison_list,
    comparison_list_order = comparison_list_order,
    go_annotations = go_annotations
  )
  
  if (!is.null(go_annotations)) {
    rle_go_term <- rle(heatmap_data_result$heatmap_data$GOterm)
    run_lengths <- rle_go_term$lengths
    rowsep <- cumsum(run_lengths)
    writeLines(unique(heatmap_data_result$heatmap_data$GOterm), file(paste(dirname(filename), 'go_terms.txt', sep = '/')))
  } else {
    rowsep <- NA
  }
  
  if (label_rows) {
    labels <- heatmap_data_result$heatmap_data$label
  } else {
    labels <- NA
  }
  
  pdf(filename, width = width, height = height)
  heatmap.2(heatmap_data_result$heatmap_data_matrix,
            scale = "none", 
            dendrogram = "none", 
            # key = FALSE,
            Colv = NA,
            Rowv = NA,
            col = colorRampPalette(c("blue","white","red")),
            trace = "none", 
            density.info = "none",
            labRow = labels,
            labCol = c(rep(c("D", "SM", "SS"), 3), rep(c("C", "M", "Inf", "Out"), 2)),
            cexCol = 0.5,
            cexRow = 0.4,
            offsetCol = 0.15,
            colsep = c(3,6,9,13),
            rowsep = rowsep,
            symkey = FALSE,
            lhei = c(1, 10),
            lwid = c(1, 6),
            keysize=0.1,
            key.par = list(cex=0.6),
            margins = c(3, 15),
            na.color = "darkgrey",
            key.title = "log2 FC")
  dev.off()
}
