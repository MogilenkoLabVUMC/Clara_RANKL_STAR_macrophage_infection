# Function to install and load packages
install_and_load <- function(package_name, from_bioc = FALSE) {
    if (!require(package_name, character.only = TRUE)) {
        if (from_bioc) {
            if (!requireNamespace("BiocManager", quietly = TRUE)) {
                install.packages("BiocManager")
            }
            BiocManager::install(package_name)
        } else {
            install.packages(package_name)
        }
    }
    library(package_name, character.only = TRUE)
}



# Install and load necessary packages
packages <- list(
    list("ggplot2", FALSE),
    list("tidyverse", FALSE),
    list("pheatmap", FALSE),
    list("Cairo", FALSE),
    list("edgeR", TRUE),
    list("limma", TRUE),
    list("biomaRt", TRUE),
    list("EnhancedVolcano", TRUE),
    list("EGSEA", TRUE),
    list("clusterProfiler", TRUE)
)

# Iterate over the package list to install and load
for (pkg in packages) {
    install_and_load(pkg[[1]], from_bioc = pkg[[2]])
}
