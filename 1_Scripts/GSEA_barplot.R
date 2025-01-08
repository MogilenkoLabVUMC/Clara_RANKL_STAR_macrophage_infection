GSEA_barplot <- function(gsea_obj, title = "GSEA NES Barplot", q_cut = 0.05, top_n = 30, width = 10, height = 8) {
    # Extract significant pathways
    gsea_data <- as.data.frame(gsea_obj@result) %>%
        filter(qvalue < q_cut) %>%
        arrange(desc(abs(NES))) %>%
        head(top_n) %>%
        # Remove the smart_wrap call here since ggplot2 will handle text wrapping
        arrange(NES)
    
    p <- ggplot(gsea_data, aes(x = reorder(Description, NES), y = NES)) +
        geom_bar(stat = "identity", aes(fill = NES > 0)) +
        scale_fill_manual(values = c("skyblue", "orange")) +
        coord_flip() +
        labs(title = title,
             x = NULL,
             y = "Normalized Enrichment Score") +
        custom_minimal_theme_with_grid() +
        # Add text wrapping directly in theme
        theme(
            legend.position = "none",
            axis.text.y = element_text(lineheight = 0.8)
        )
    
    filename <- paste0("3_Results/imgs/GSEA/", gsub(" ", "_", title), "_barplot.pdf")
    ggsave(filename = filename, plot = p, width = width, height = height)
    
    return(p)
}
