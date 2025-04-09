library(pheatmap)
library(dplyr)
library(stringr)  # for string manipulation

plot_single_pathway_heatmap <- function(
  gsea_obj,
  pathway_name,
  expression_data,
  sample_order,
  annotation_col = NULL,
  annotation_colors = NULL,
  output_prefix = "pathway_heatmap",
  scale_expr = "row"
) {
  # Convert GSEA results to a data frame
  gsea_df <- as.data.frame(gsea_obj@result)
  
  # Find matching pathway
  pathway_row <- gsea_df %>%
    dplyr::filter(
      grepl(pathway_name, Description, ignore.case = TRUE) | 
      grepl(pathway_name, ID, ignore.case = TRUE)
    )
  
  if (nrow(pathway_row) == 0) {
    message("No matching pathway found for: ", pathway_name)
    return(NULL)
  }
  
  # Extract the first matching row
  pathway_info <- pathway_row[1, ]
  
  # Extract core genes
  core_genes <- unlist(strsplit(pathway_info$core_enrichment, "/"))
  core_genes <- core_genes[core_genes %in% rownames(expression_data)]
  
  if (length(core_genes) == 0) {
    message("No genes in ", pathway_name, " found in expression data.")
    return(NULL)
  }
  
  # Subset expression matrix
  expr_sub <- expression_data[core_genes, sample_order, drop = FALSE]
  
  # Clean up the pathway name for the heatmap title
  # e.g. remove prefixes like "GOBP_", "KEGG_", etc., then wrap text
  cleaned_pathway_name <- pathway_info$Description %>%
    str_replace_all(c("^GOBP_" = "", "^KEGG_" = "", "^REACTOME_" = "", "^HALLMARK_" = "")) %>%
    str_replace_all("_", " ") %>% 
    str_wrap(width = 40)  # wrap text to 40 characters per line
  
  heatmap_title <- paste("Heatmap of", cleaned_pathway_name)
  
  # Plot with pheatmap and smaller fonts
  p <- pheatmap::pheatmap(
    mat = expr_sub,
    scale = scale_expr,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    color = colorRampPalette(c("navy", "white", "red"))(50),
    main = heatmap_title,
    fontsize = 9,           # overall font size
    fontsize_row = 7,       # row labels
    fontsize_col = 7        # column labels
  )
  
  # Save the plot
  out_filename <- paste0(
    output_prefix, "_", 
    gsub("[^A-Za-z0-9_]+", "_", pathway_info$Description),  # a safe file name
    ".pdf"
  )
  
  ggsave(out_filename, plot = p[[4]], width = 10, height = 8)
  # p[[4]] is the ggplot object for the heatmap in pheatmap 1.0+ 
  # (versions may vary—if in doubt, wrap the entire pheatmap() call in a pdf()... dev.off() block)
  
  return(invisible(p))
}
