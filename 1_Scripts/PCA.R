source("1_Scripts/custom_minimal_theme.R")

create_pca_plot <- function(DGE_object, title = "PCA Plot") {
  # Calculate logCPM
  logCPM <- cpm(DGE_object, log = TRUE, prior.count = 1)

  # Calculate PCA
  pca <- prcomp(t(logCPM))

  # Prepare PCA data frame
  pca_data <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    condition = DGE_object$samples$group,
    sample_names = paste0(
      # Time point
      ifelse(grepl("4h", DGE_object$samples$group), "4h", "24h"),
      "_",
      # Cell type
      ifelse(grepl("STm", DGE_object$samples$group), "STm", "mock"),
      "_",
      # RANKL concentration
      ifelse(grepl("100", DGE_object$samples$group), "100", "0")
    )
  )

  # Calculate variance percentages
  percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100


  # Calculate axis limits with 10% buffer
  x_max <- max(abs(pca_data$PC1))
  y_max <- max(abs(pca_data$PC2))
  limit_buffer <- 1.5 # 50% buffer

  # Create and return the plot
  pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 5) +
    geom_text(aes(label = sample_names),
      vjust = -1, hjust = 1,
      size = 3,
      show.legend = FALSE
    ) +
    labs(color = "Condition", shape = "Condition") +
    xlab(paste0("PC1: ", round(percentVar[1], 1), "% variance")) +
    ylab(paste0("PC2: ", round(percentVar[2], 1), "% variance")) +
    scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                             "#0072B2", "#D55E00", "#CC79A7", "#000000")) +
    ggtitle(title) +
    custom_minimal_theme_with_grid() +
    xlim(-x_max * limit_buffer, x_max * limit_buffer) +
    ylim(-y_max * limit_buffer, y_max * limit_buffer) +
      theme(
      legend.position = "inside",
      legend.position.inside = c(0.7, 0.9),
      legend.justification = c(0, 1),
      legend.box.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
      legend.box.margin = margin(2, 2, 2, 2),
      plot.title = element_text(hjust = 0.2, vjust = 0.2, margin = margin(b = 10)),
      plot.margin = margin(20, 20, 20, 20)
    )

  # Save the plot
  filename <- paste0("3_Results/imgs/PCA/", gsub(" ", "_", title), ".pdf")
  ggsave(
    filename = filename,
    plot = pca_plot,
    width = 7,
    height = 5,
    dpi = 300
  )

  return(pca_plot)
}
