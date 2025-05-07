GSEA_faceted_dotplot <- function(gsea_obj, showCategory = 10, font.size = 10, title = "GSEA Faceted Dotplot",
                                replace_ = TRUE, capitalize_1 = TRUE, capitalize_all = FALSE,
                                strip_prefix = TRUE,
                                filterBy = "qvalue",  # Added parameter
                                sortBy = "GeneRatio", # Added parameter
                                q_cut = 0.05) {       # Added parameter
    
    # Process the result data
    gsea_data <- as.data.frame(gsea_obj@result)
    
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
    
    # Filter for significant pathways
    gsea_data_filtered <- gsea_data %>%
        filter(qvalue < q_cut)
    
    # Calculate GeneRatio if needed for sorting
    if (sortBy == "GeneRatio" && !"GeneRatio" %in% colnames(gsea_data_filtered)) {
        # Calculate the gene count from 'core_enrichment' by counting '/' and adding 1
        gene_count <- gsea_data_filtered %>%
            group_by(ID) %>%
            summarise(count = sum(str_count(core_enrichment, "/")) + 1)
        
        # Merge gene counts with the filtered data and calculate GeneRatio
        gsea_data_filtered <- left_join(gsea_data_filtered, gene_count, by = "ID") %>%
            mutate(GeneRatio = count / setSize)
    }
    
    # First filter by NES direction
    pos_data <- gsea_data_filtered %>% filter(NES > 0)
    neg_data <- gsea_data_filtered %>% filter(NES < 0)
    
    # Then apply sorting based on sortBy parameter
    if (sortBy == "GeneRatio") {
        pos_data <- pos_data %>% arrange((GeneRatio)) %>% head(showCategory)
        neg_data <- neg_data %>% arrange((GeneRatio)) %>% head(showCategory)
    } else if (sortBy == "qvalue") {
        pos_data <- pos_data %>% arrange(qvalue) %>% head(showCategory)
        neg_data <- neg_data %>% arrange(qvalue) %>% head(showCategory)
    } else if (sortBy == "NES") {
        pos_data <- pos_data %>% arrange(desc(NES)) %>% head(showCategory)
        neg_data <- neg_data %>% arrange(NES) %>% head(showCategory)
    } else {
        warning("Invalid sortBy parameter. Defaulting to GeneRatio.")
        # If GeneRatio doesn't exist, calculate it
        if (!"GeneRatio" %in% colnames(pos_data)) {
            gene_count <- rbind(pos_data, neg_data) %>%
                group_by(ID) %>%
                summarise(count = sum(str_count(core_enrichment, "/")) + 1)
            
            pos_data <- left_join(pos_data, gene_count, by = "ID") %>%
                mutate(GeneRatio = count / setSize)
            neg_data <- left_join(neg_data, gene_count, by = "ID") %>%
                mutate(GeneRatio = count / setSize)
        }
        
        pos_data <- pos_data %>% arrange(desc(GeneRatio)) %>% head(showCategory)
        neg_data <- neg_data %>% arrange(desc(GeneRatio)) %>% head(showCategory)
    }
    
    # Add a factor for ordering in the plot
    pos_data$order <- 1:nrow(pos_data)
    neg_data$order <- 1:nrow(neg_data)
    
    # Combine data
    plot_data <- rbind(pos_data, neg_data) %>%
        mutate(Direction = ifelse(NES > 0, "Upregulated", "Downregulated"))
    
    # Create faceted plot - use the order column instead of NES for reordering
    p <- ggplot(plot_data, aes(x = abs(NES), y = reorder(Description, order))) +
        geom_point(aes(size = -log10(qvalue), color = Direction)) +
        scale_color_manual(values = c("Upregulated" = "orange", "Downregulated" = "skyblue")) +
        facet_grid(Direction ~ ., scales = "free_y", space = "free") +
        labs(x = "Normalized Enrichment Score", y = NULL, title = title) +
        theme_minimal() +
        theme(
            text = element_text(size = font.size),
            plot.title = element_text(hjust = 0.5, size = font.size + 2)
        )
    
    return(p)
}
