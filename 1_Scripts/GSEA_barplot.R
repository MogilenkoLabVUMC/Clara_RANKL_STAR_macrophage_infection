GSEA_barplot <- function(
    gsea_obj, 
    title = "GSEA NES Barplot", 
    q_cut = 0.05, top_n = 30,
    width = 10, height = 8,
    output_dir = "3_Results/imgs/GSEA/",
    font.size = 10,
    save_plot = TRUE,
    dpi = 300,
    replace_ = TRUE,
    capitalize_1 = TRUE,
    capitalize_all = FALSE,
    strip_prefix = TRUE) {
    
    # Extract significant pathways
    gsea_data <- as.data.frame(gsea_obj@result) %>%
        filter(qvalue < q_cut) %>%
        arrange(desc(abs(NES))) %>%
        head(top_n) %>%
        arrange(NES)
    
    # Modify Description field based on function arguments
    if (replace_) {
        gsea_data$Description <- gsea_data$Description %>%
            str_replace_all("_", " ") # Replace "_" with " " if 'replace' is TRUE
    }
    
    # Strip common prefixes if requested - BEFORE capitalization
    if (strip_prefix) {
        common_prefixes <- c(
            "HALLMARK ", "KEGG ", "REACTOME ", "BIOCARTA ", "GOBP ", "GOCC ", "GOMF ", "MEDICUS ", "GRTD",
            "PID ", "WIKIPATHWAY ", "GO ", "TFT", "CGP"
        )
        
        # Create a pattern for all prefixes at once for more efficient replacement
        prefix_pattern <- paste0("^(", paste(common_prefixes, collapse = "|"), ")")
        gsea_data$Description <- str_replace(gsea_data$Description, prefix_pattern, "")
    }
    
    if (capitalize_1) {
        gsea_data$Description <- gsea_data$Description %>%
            str_to_sentence() # Capitalize only the first word if 'capitalize_1' is TRUE
    }
    
    if (capitalize_all) {
        gsea_data$Description <- gsea_data$Description %>%
            str_to_title() # Capitalize all words if 'capitalize_all' is TRUE
    }
    
    # Add smart wrapping for long descriptions
    smart_wrap <- function(text, width = 40) {
        words <- unlist(strsplit(text, " "))
        total_chars <- nchar(text)
        
        if (total_chars > width) {
            mid_point <- length(words) %/% 2
            first_half <- paste(words[1:mid_point], collapse = " ")
            second_half <- paste(words[(mid_point + 1):length(words)], collapse = " ")
            return(paste(first_half, second_half, sep = "\n"))
        }
        return(text)
    }
    
    gsea_data$Description <- sapply(gsea_data$Description, smart_wrap)
    
    p <- ggplot(gsea_data, aes(x = reorder(Description, NES), y = NES)) +
        geom_bar(stat = "identity", aes(fill = NES > 0)) +
        scale_fill_manual(values = c("skyblue", "orange")) +
        coord_flip() +
        labs(title = title,
             x = NULL,
             y = "Normalized Enrichment Score") +
        custom_minimal_theme_with_grid() +
        theme(
            panel.background = element_rect(fill = "white", color = NA), # White background with black border
            plot.title = element_text(hjust = 0.5, size = font.size + 2),
            axis.text.x = element_text(size = font.size),
            axis.text.y = element_text(size = font.size, hjust = 1, lineheight = 0.8),
            legend.position = "none"
        )

    # Save the plot if save_plot is TRUE
    if (save_plot) {
        # Create directory if it doesn't exist
        if (!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }
        
        # Create filename from title
        filename <- file.path(output_dir, paste0(gsub(" ", "_", title), "_barplot.pdf"))
        
        ggsave(
            filename = filename,
            plot = p,
            width = width,
            height = height,
            dpi = dpi
        )
    }
    
    return(p)
}
