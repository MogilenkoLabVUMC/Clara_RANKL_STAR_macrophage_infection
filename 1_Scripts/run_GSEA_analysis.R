source("/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/1_Scripts/GSEA_barplot.R")
source("/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/1_Scripts/runGSEA.R")
source("/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/1_Scripts/GSEA_dotplot.R")
source("/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/1_Scripts/GSEA_dotplot_facet.R")  # Add this line to source the facet dotplot script

# Function to run GSEA and create plots for one DE table
run_gsea_analysis <- function(
    de_table, 
    timepoint, 
    species = "Mus musculus",
    n_cat = 15,
    output_dir = "3_Results/imgs/GSEA/") {
  
  cat("\n========================================\n")
  cat("Starting GSEA analysis for timepoint:", timepoint, "\n")
  cat("========================================\n")
  
  # Define database-specific plot parameters
  db_plot_params <- list(
      HALLMARK = list(width = 6, height = 8, font.size = 8),
      REACTOME = list(width = 8, height = 8, font.size = 8),
      KEGG = list(width = 7, height = 8, font.size = 8),
      GOBP = list(width = 7, height = 8, font.size = 8),
      GOCC = list(width = 7, height = 8, font.size = 8),
      GOMF = list(width = 7, height = 8, font.size = 8),
      GRTD = list(width = 7, height = 8, font.size = 8),
      CGP = list(width = 7, height = 8, font.size = 8),
      PID = list(width = 7, height = 8, font.size = 8),
      WIKIPATHWAY = list(width = 7, height = 8, font.size = 8),
      TFT = list(width = 7, height = 8, font.size = 8),
      BIOCARTA = list(width = 7, height = 8, font.size = 8)
  )
  
  # Create the output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
      cat("Creating output directory:", output_dir, "\n")
      dir.create(output_dir, recursive = TRUE)
  } else {
      cat("Output directory already exists:", output_dir, "\n")
  }
  
  # A list to store all GSEA results by database
  result_list <- list()
  
  # Print DE table info
  cat("Input DE table dimensions:", nrow(de_table), "rows x", ncol(de_table), "columns\n")
  cat("Species:", species, "\n")
  cat("Number of categories to show:", n_cat, "\n\n")
  
  # Loop through each database config
  for (db_name in names(db_configs)) {
      cat("\n------------------------------------\n")
      cat("Processing database:", db_name, "\n")
      cat("------------------------------------\n")
      
      params <- db_plot_params[[db_name]]
      cat("Plot parameters - Width:", params$width, "Height:", params$height, "Font size:", params$font.size, "\n")
      
      # Run GSEA
      cat("Running GSEA analysis...\n")
      gsea_result <- runGSEA(
          DE_results = de_table,
          rank_metric = "t",
          species = species,
          category = db_configs[[db_name]]$category,
          subcategory = db_configs[[db_name]]$subcategory,
          padj_method = "fdr",
          nperm = 100000,
          pvalue_cutoff = 0.05
      )
      
      # Print GSEA result summary
      cat("GSEA completed. Result summary:\n")
      cat("  Total pathways analyzed:", nrow(gsea_result@result), "\n")
      sig_pathways <- sum(gsea_result@result$qvalue < 0.05)
      cat("  Significant pathways (q < 0.05):", sig_pathways, "\n")
      up_pathways <- sum(gsea_result@result$NES > 0 & gsea_result@result$qvalue < 0.05)
      down_pathways <- sum(gsea_result@result$NES < 0 & gsea_result@result$qvalue < 0.05)
      cat("  Upregulated pathways:", up_pathways, "\n")
      cat("  Downregulated pathways:", down_pathways, "\n")
      
      # Store the GSEA result in our output list
      result_list[[db_name]] <- gsea_result
      
      # Generate Upregulated (NES>0) dotplot
      cat("\nGenerating upregulated pathways dotplot...\n")
      up_title <- paste0(timepoint, " ", db_name, " Upregulated")
      up_filename <- file.path(output_dir, paste0(gsub(" ", "_", up_title), ".pdf"))
      cat("  Output file:", up_filename, "\n")
      
      GSEA_dotplot(
          gsea_result,
          filterBy = "NES_positive",
          sortBy = "GeneRatio",
          font.size = params$font.size,
          showCategory = n_cat,
          q_cut = 0.05,
          replace_ = TRUE,
          capitalize_1 = FALSE,
          capitalize_all = FALSE,
          strip_prefix = TRUE,
          min.dotSize = 2,
          title = up_title,
          output_dir = output_dir,
          width = params$width,
          height = params$height
      )
      
      # Generate Downregulated (NES<0) dotplot
      cat("\nGenerating downregulated pathways dotplot...\n")
      down_title <- paste0(timepoint, " ", db_name, " Downregulated")
      down_filename <- file.path(output_dir, paste0(gsub(" ", "_", down_title), ".pdf"))
      cat("  Output file:", down_filename, "\n")
      
      GSEA_dotplot(
          gsea_result,
          filterBy = "NES_negative",
          sortBy = "GeneRatio",
          font.size = params$font.size,
          showCategory = n_cat,
          q_cut = 0.05,
          replace_ = TRUE,
          capitalize_1 = FALSE,
          capitalize_all = FALSE,
          strip_prefix = TRUE,
          min.dotSize = 2,
          title = down_title,
          output_dir = output_dir,
          width = params$width,
          height = params$height
      )
      
      # Generate NES barplot
      cat("\nGenerating NES barplot...\n")
      bar_title <- paste0(timepoint, " ", db_name, " NES")
      bar_filename <- file.path(output_dir, paste0(gsub(" ", "_", bar_title), "_barplot.pdf"))
      cat("  Output file:", bar_filename, "\n")
      
      GSEA_barplot(
          gsea_result,
          top_n = n_cat * 3,
          title = bar_title,
          width = params$width,
          height = params$height,
          font.size = params$font.size,
          output_dir = output_dir,
          replace_ = TRUE,
          capitalize_1 = FALSE,
          capitalize_all = FALSE,
          strip_prefix = TRUE
      )
      
      # Generate faceted dotplot
      cat("\nGenerating faceted dotplot...\n")
      facet_title <- paste0(timepoint, " ", db_name, " Pathways")
      facet_filename <- file.path(output_dir, paste0(gsub(" ", "_", timepoint), "_", db_name, "_faceted_dotplot.pdf"))
      cat("  Output file:", facet_filename, "\n")
      
      facet_plot <- GSEA_faceted_dotplot(
          gsea_result,
          showCategory = n_cat,
          font.size = params$font.size,
          title = facet_title,
          replace_ = TRUE,
          capitalize_1 = FALSE,
          capitalize_all = FALSE,
          strip_prefix = TRUE
      )
      
      # Save the faceted dotplot
      ggsave(
          filename = facet_filename,
          plot = facet_plot,
          width = params$width,
          height = params$height,  # Slightly taller for the faceted plot
          dpi = 300
      )
      
      cat("\nCompleted processing for database:", db_name, "\n")
  }
  
  cat("\n========================================\n")
  cat("GSEA analysis completed for timepoint:", timepoint, "\n")
  cat("Total databases processed:", length(result_list), "\n")
  cat("========================================\n\n")
  
  # Return the list of GSEA results for this timepoint
  return(result_list)
}
