process_rnaseq_data <- function(counts, sample_info, annotate = TRUE) {
    # Function to extract number from sample name
    extract_num <- function(x) {
        as.numeric(sub(".*\\.(\\d+)\\..*", "\\1", x))
    }
    
    # Sort sample names naturally by their numeric part
    sample_order <- order(sapply(colnames(counts), extract_num))
    counts <- counts[, sample_order]
    
    # Ensure sample info is in the same order
    sample_info <- sample_info[colnames(counts), ]
    
    # Create DGEList object
    DGE <- DGEList(counts = counts)
    
    # Add sample information
    DGE$samples$timepoint <- factor(sample_info$Time)
    DGE$samples$treatment <- factor(sample_info$Treatment)
    DGE$samples$rankl <- factor(sample_info$RANKL)
    # Create syntactically valid group names by adding 't' prefix to numbers
    DGE$samples$group <- factor(gsub("^(\\d+)h", "t\\1h", sample_info$Group))
    
    # Filter low expressed genes
    keep <- filterByExpr(DGE)
    DGE <- DGE[keep, , keep.lib.sizes = FALSE]
    
    # Normalize library sizes
    DGE <- calcNormFactors(DGE, method = "TMM")
    
    # Add gene annotations if requested
    if (annotate) {
        require(biomaRt)
        ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
        gene_ids <- rownames(DGE)
        
        annotations <- getBM(
            attributes = c("ensembl_gene_id", "external_gene_name", "description"),
            filters = "ensembl_gene_id",
            values = gene_ids,
            mart = ensembl
        )
        
        # Keep only unique gene names
        mapping_table <- annotations %>%
            distinct(external_gene_name, .keep_all = TRUE) %>%
            select(ensembl_gene_id, external_gene_name)
        
        # Match and rename genes
        matched_idx <- match(rownames(DGE), mapping_table$ensembl_gene_id)
        rownames(DGE) <- ifelse(!is.na(matched_idx), 
                               mapping_table$external_gene_name[matched_idx],
                               rownames(DGE))
    }
    
    return(DGE)
}
