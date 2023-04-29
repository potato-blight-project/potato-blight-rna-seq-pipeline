#' Plots a heatmap given a definition.
#' 
#' @param heatmap_definition A list containing the specification of the heatmap.
#' @param annotated_degs A list of annotated DEGs (see go_and_func_enrichment).
#' @param comparison_list A vector of items to filter the heatmap by.
#' @param filter_by The column in the DEGs to filter by (against the comparison_list).
#' @param filter_fn A function to filter rows.
#' @param comparison_list_order The order in which to plot the grouped items.
#' @param go_annotations A dataframe with GO annotations against transcript IDs.
#' @param label_rows Flag indicating whether to label rows with gene names.
#' @param filename The file path to save the PDF.
#' @param width The width of the heatmap.
#' @param height The height of the heatmap.
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
  
  # 1. Creates the column blocks from the heatmap definition list
  blocks <- list()
  
  for (block_name in names(heatmap_definition)) {
    block_def <- heatmap_definition[[block_name]]
    degs_list <- list()
    
    for (comparison in block_def$comparisons) {
      degs_list[[comparison]] <- filter_by_list(annotated_degs, comparison, comparison_list, filter_by)
    }
    
    blocks[[block_name]] <- create_block(degs_list, block_def$suffixes)
  }
  
  # 2. Creates the heatmap matrix
  heatmap_data_result <- create_heatmap_matrix(
    block_list = blocks,
    suffixes = unlist(lapply(names(blocks), function(x) { paste0('_', x) })),
    filter_fn = filter_fn,
    comparison_list = comparison_list,
    comparison_list_order = comparison_list_order,
    go_annotations = go_annotations
  )
  
  heatmap_filename <- basename(filename)
  heatmap_filename <- tools::file_path_sans_ext(heatmap_filename)
  
  # 3. Separate rows by GO terms if the annotations are provided
  if (!is.null(go_annotations)) {
    rle_go_term <- rle(heatmap_data_result$heatmap_data$GOterm)
    run_lengths <- rle_go_term$lengths
    rowsep <- cumsum(run_lengths)
    
    go_terms_filename <- paste(heatmap_filename, 'go_terms', sep = '_')
    writeLines(unique(heatmap_data_result$heatmap_data$GOterm), file(paste(dirname(filename), go_terms_filename, sep = '/')))
  } else {
    rowsep <- NA
  }
  
  # 4. Save the heatmap matrix as a CSV
  csv_filename <- paste(heatmap_filename, 'csv', sep = '.')
  write.csv(heatmap_data_result$heatmap_data, paste(dirname(filename), csv_filename, sep = '/'), row.names=FALSE)
  
  # 5. Plot the heatmap
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
            labCol = c(rep(c("D", "SM", "SS"), 2), rep(c("C", "M", "Inf", "Out"), 2)),
            cexCol = 0.5,
            cexRow = 0.25,
            offsetCol = 0.15,
            colsep = c(3, 6, 10),
            rowsep = rowsep,
            symkey = FALSE,
            lhei = c(1, 10),
            lwid = c(1, 6),
            keysize=0.1,
            key.par = list(cex=0.6),
            margins = c(3, 15),
            na.color = "darkgrey",
            key.title = "log2 FC")
  
  axis(
    side = 3,
    at = c(0.12, 0.25, 0.42, 0.6), #c(0.1, 0.22, 0.33, 0.47, 0.62),
    labels = unlist(lapply(heatmap_definition, function(x) {x$name})),
    tick = FALSE,
    cex.axis = 0.52)
  
  dev.off()
}


#' Filter the annotated DEGs
#' 
#' @param annotated_degs The input annotated DEGs.
#' @param comparison The comparison string in the DEGs (e.g. "DEG.Duke.Out").
#' @param comparison_list The comparison list (e.g. GO terms).
#' @param filter_by The column in the DEGs to filter by (against the comparison_list).
#' @return The filtered annotated DEGs
filter_by_list <- function(
  annotated_degs,
  comparison,
  comparison_list,
  filter_by
) {
  annotated_degs[[comparison]] %>%
    filter(.data[[filter_by]] %in% comparison_list) %>%
    dplyr::select(-transcript_id, -pvalue, -stat, -baseMean, -TAIR, -GOterm) %>%
    dplyr::filter(!is.na(name)) %>%
    distinct(gene, .keep_all = TRUE)
}


#' Given a list of DEGs, it merges them by gene
#' 
#' @param degs_list A list of DEG dataframes
#' @param suffixes The suffixes for each merged dataframe
#' @return The merged blocks
create_block <- function(degs_list, suffixes) {
  purrr::map2(degs_list, suffixes, ~{
    cols_to_rename <- setdiff(names(.x), 'gene')
    purrr::set_names(.x, c('gene', paste0(cols_to_rename, .y)))
  }) %>% purrr::reduce(merge, by = 'gene', all = TRUE)
}


#' Create the heatmap data matrix that will be plotted with 'heatmap.2'
#' 
#' @param block_list A list of blocks (i.e. merged DEG dataframes)
#' @param suffixes The suffixes for each merged dataframe
#' @param filter_fn A function to filter rows
#' @param comparison_list A vector of items to filter the heatmap by
#' @param comparison_list_order The order in which to plot the grouped items
#' @param go_annotations A dataframe with GO annotations against transcript IDs
#' @return A list with the heatmap data matrix and the plain numeric heatmap matrix
create_heatmap_matrix <- function(
  block_list,
  suffixes,
  filter_fn,
  comparison_list,
  comparison_list_order,
  go_annotations
) {
  # 1. Merge all blocks into one dataframe
  heatmap_data <- create_block(block_list, suffixes)
  
  # 2. Apply row filter (e.g. constitutive filter)
  heatmap_data_idx <- heatmap_data %>%
    dplyr::select(starts_with('log2FoldChange')) %>%
    mutate(id = row_number()) %>%
    filter(filter_fn(.)) %>%
    pull(id)
  
  heatmap_data <- heatmap_data %>% slice(heatmap_data_idx)
  
  # 3. Create label column
  func_names <- heatmap_data %>%
    dplyr::select(starts_with('name')) %>%
    rowwise() %>%
    mutate(func_name = list(first(na.omit(c_across(everything()))))) %>%
    dplyr::select(func_name)
  func_names <- unlist(lapply(func_names$func_name, function (x) { ifelse(length(x) == 0, NA, x) }))
  heatmap_data$label <- paste(heatmap_data$gene, func_names, sep = '_')
  heatmap_data <- heatmap_data %>% relocate(label, .after = gene)
  
  heatmap_data <- heatmap_data %>% dplyr::select(gene, label, starts_with('log2FoldChange'))
  
  # 4. Annotate heatmap data with GO terms and order
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
    
    # heatmap_data$label <- paste(heatmap_data$GOterm, heatmap_data$label, sep = '_')
  }
  
  return(list(
    heatmap_data = heatmap_data,
    heatmap_data_matrix = as.matrix(heatmap_data %>% dplyr::select(starts_with('log2FoldChange')))
  ))
}
