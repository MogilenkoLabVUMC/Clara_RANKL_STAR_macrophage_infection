# Add this new function to combine volcano plots
create_combined_volcano_plots <- function(DE_rankl_4h, DE_rankl_24h, DE_rankl_time, 
                                        p_cutoff = 0.05, fc_cutoff = 2.0) {
    # Create all three plots but don't save them individually
    p1 <- create_volcano_plot(DE_rankl_4h,
        p_cutoff = p_cutoff, fc_cutoff = fc_cutoff,
        x_breaks = 2, max.overlaps = 10, label_method = "sig",
        title = "4 hours"
    )
    
    p2 <- create_volcano_plot(DE_rankl_24h,
        p_cutoff = p_cutoff, fc_cutoff = fc_cutoff,
        x_breaks = 2, max.overlaps = 8, label_method = "sig",
        title = "24 hours"
    )
    
    p3 <- create_volcano_plot(DE_rankl_time,
        p_cutoff = p_cutoff, fc_cutoff = fc_cutoff,
        x_breaks = 2, max.overlaps = 10, label_method = "sig",
        title = "24h vs 4h"
    )
    
    # Find common axis limits
    y_max <- max(c(
        max(-log10(DE_rankl_4h$P.Value), na.rm = TRUE),
        max(-log10(DE_rankl_24h$P.Value), na.rm = TRUE),
        max(-log10(DE_rankl_time$P.Value), na.rm = TRUE)
    ))
    
    x_max <- max(abs(c(
        DE_rankl_4h$logFC,
        DE_rankl_24h$logFC,
        DE_rankl_time$logFC
    )), na.rm = TRUE)
    x_max <- ceiling(x_max)
    
    # Apply common limits to all plots
    p1 <- p1 + coord_cartesian(xlim = c(-x_max, x_max), ylim = c(0, y_max))
    p2 <- p2 + coord_cartesian(xlim = c(-x_max, x_max), ylim = c(0, y_max))
    p3 <- p3 + coord_cartesian(xlim = c(-x_max, x_max), ylim = c(0, y_max))
    
    # Combine plots using patchwork
    combined_plot <- p1 + p2 + p3 + 
        patchwork::plot_layout(ncol = 3, guides = "collect")
    
    # Save the combined plot
    ggsave(
        filename = "3_Results/imgs/volcano/combined_volcano_plots.pdf",
        plot = combined_plot,
        width = 21,
        height = 7,
        dpi = 300
    )
    
    return(combined_plot)
}
