process_rnaseq_data <- function(file_path, annotate = FALSE, normalize = TRUE) {
    # Read the data with check.names = FALSE to prevent X prefix
    raw_data <- read.csv(file_path, 
                        header = TRUE,
                        row.names = 1,
                        check.names = FALSE)
    
    # Create unique column names while preserving original names
    original_names <- unique(colnames(raw_data))
    groups <- factor(colnames(raw_data), levels = original_names)
    
    # Make column names unique without X prefix
    new_colnames <- character(length(colnames(raw_data)))

    for (name in original_names) {
        indices <- which(colnames(raw_data) == name)
        new_colnames[indices] <- paste0(name, "_", seq_along(indices))
    }

    colnames(raw_data) <- new_colnames
    
    # Extract timepoints cleanly
    groups <- factor(sub("_\\d+$", "", new_colnames))
    
    # Rest of the function remains the same...
    if (annotate) {
        require(biomaRt)
        ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
        gene_ids <- rownames(raw_data)
        
        annotations <- getBM(
            attributes = c("ensembl_transcript_id", "external_gene_name", "description"),
            filters = "ensembl_transcript_id",
            values = gene_ids,
            mart = ensembl
        )
        
        mapping_table <- annotations %>%
            distinct(external_gene_name, .keep_all = TRUE) %>%
            select(ensembl_transcript_id, external_gene_name)
        
        raw_data <- raw_data[match(mapping_table$ensembl_transcript_id, rownames(raw_data)), ]
        rownames(raw_data) <- mapping_table$external_gene_name
    }
    
    require(edgeR)
    DGE <- DGEList(counts = raw_data, group = groups)
   # DGE$samples$timepoint <- factor(timepoints)
    
    keep <- filterByExpr(DGE)
    DGE <- DGE[keep, , keep.lib.sizes = FALSE]
    
    if (normalize) {
        DGE <- normLibSizes(DGE, method = "TMM")
    }
    
    return(DGE)
}
