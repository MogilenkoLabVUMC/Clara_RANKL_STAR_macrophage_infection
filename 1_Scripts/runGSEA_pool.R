source("1_Scripts/get_pathway_genes.R")
source("1_Scripts/get_significant_pathways.R")
source("1_Scripts/get_pathway_genes_all.R")
source("1_Scripts/calculate_pathway_scores.R")
source("1_Scripts/get_significant_pathways.R")


run_pooled_gsea <- function(fit, contrasts, DGErankl, top_n = 30) {
    # Run GSEA for all pool contrasts
    gsea_results <- list()
    contrasts_to_analyze <- colnames(contrasts)
    
    # Run GSEA analysis for each contrast
    for(contrast in contrasts_to_analyze) {
        de_results <- topTable(fit, coef = contrast, 
                             sort.by = "t", adjust.method = "fdr", n = Inf)
        de_results <- de_results[rownames(de_results) != "", ]
        
        tryCatch({
            gsea_results[[contrast]][["hallmark"]] <- runGSEA(
                de_results, rank_metric = "t", species = "Mus musculus",
                category = "H", padj_method = "fdr", nperm = 100000,
                pvalue_cutoff = 0.01
            )
            
            gsea_results[[contrast]][["kegg"]] <- runGSEA(
                de_results, rank_metric = "t", species = "Mus musculus",
                category = "C2", subcategory = "CP:KEGG", padj_method = "fdr",
                nperm = 100000, pvalue_cutoff = 0.01
            )
            
            gsea_results[[contrast]][["gobp"]] <- runGSEA(
                de_results, rank_metric = "t", species = "Mus musculus",
                category = "C5", subcategory = "GO:BP", padj_method = "fdr",
                nperm = 100000, pvalue_cutoff = 0.01
            )
            
            gsea_results[[contrast]][["reactome"]] <- runGSEA(
                de_results, rank_metric = "t", species = "Mus musculus",
                category = "C2", subcategory = "CP:REACTOME", padj_method = "fdr",
                nperm = 100000, pvalue_cutoff = 0.01
            )
        }, error = function(e) {
            message(sprintf("Error in contrast %s: %s", contrast, e$message))
        })
    }
    
    # Get pools for each database
    pools <- list(
        hallmark = get_significant_pathways(lapply(gsea_results, `[[`, "hallmark")),
        kegg = get_significant_pathways(lapply(gsea_results, `[[`, "kegg")),
        gobp = get_significant_pathways(lapply(gsea_results, `[[`, "gobp")),
        reactome = get_significant_pathways(lapply(gsea_results, `[[`, "reactome"))
    )
    
    # Get pathway genes
    genes <- list(
        hallmark = get_pathway_genes_all(gsea_results, "hallmark", top = top_n),
        kegg = get_pathway_genes_all(gsea_results, "kegg", top = top_n),
        gobp = get_pathway_genes_all(gsea_results, "gobp", top = top_n),
        reactome = get_pathway_genes_all(gsea_results, "reactome", top = top_n)
    )
    
    # Calculate normalized expression scores
    norm_counts <- cpm(DGErankl, log = TRUE)
    
    scores <- list(
        hallmark = calculate_pathway_scores(norm_counts, genes$hallmark),
        kegg = calculate_pathway_scores(norm_counts, genes$kegg),
        gobp = calculate_pathway_scores(norm_counts, genes$gobp),
        reactome = calculate_pathway_scores(norm_counts, genes$reactome)
    )
    
    return(list(
        gsea_results = gsea_results,
        pools = pools,
        genes = genes,
        scores = scores
    ))
}
