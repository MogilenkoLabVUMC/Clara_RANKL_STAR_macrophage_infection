# Load libraries and functions
source("1_Scripts/Load_libraries.R")
source("1_Scripts/process_rnaseq_data.R")
source("1_Scripts/PCA.R")
source("1_Scripts/standard_volcano.R")
source("1_Scripts/GSEA_dotplot.R")
source("1_Scripts/runGSEA.R")
source("1_Scripts/runningSumGSEAplot.R")
source("1_Scripts/combined_volcano.R")

# Prepare data
## Read count data
counts <- read.delim("0_Data/STAR_Data/counts.txt", row.names = 1)

## Read sample information
sample_info <- read.delim("0_Data/STAR_Data/seq_reference.txt", row.names = 1)
colnames(sample_info) <- c("Group", "Time", "Treatment", "RANKL")
sample_info$Sample <- rownames(sample_info)

## Create DGEList object with sample information
DGErankl <- process_rnaseq_data(counts, sample_info)
# Display sample information
DT::datatable(DGErankl$samples)


# PCA
create_pca_plot(DGErankl, title = "PCA Plot")

# Set design matrix
## Set combinatorial design
design <- model.matrix(~ 0 + group, data = DGErankl$samples)
## Rename the columns to remove "group" prefix
colnames(design) <- gsub("group", "", colnames(design))
colnames(design) <- paste0("h", colnames(design))

## Set contrast matrix
### Select contrast desing
contrasts.m <- makeContrasts(
    # RANKL effect on infection response at 4h
    RANKL_effect_4h = (h4h_STm_100 - h4h_mock_100) - (h4h_STm_0 - h4h_mock_0),
    # RANKL effect on infection response at 24h
    RANKL_effect_24h = (h24h_STm_100 - h24h_mock_100) - (h24h_STm_0 - h24h_mock_0),
    # RANKL time-dependent effect 
    Time_RANKL_effect = ((h24h_STm_100 - h24h_mock_100) - (h24h_STm_0 - h24h_mock_0)) - ((h4h_STm_100 - h4h_mock_100) - (h4h_STm_0 - h4h_mock_0)),
    # Desing matrix
    levels = design
)

### Pool contrast design
contrasts.p <- makeContrasts(
    # STm effect without RANKL at 4h
    STm_4h_0 = h4h_STm_0 - h4h_mock_0,
    # STm effect with RANKL at 4h
    STm_4h_100 = h4h_STm_100 - h4h_mock_100,
    # RANKL effect without infection at 4h
    RANKL_4h_mock = h4h_mock_100 - h4h_mock_0,
    # RANKL effect with infection at 4h
    RANKL_4h_STm = h4h_STm_100 - h4h_STm_0,
    # STm effect without RANKL at 24h
    STm_24h_0 = h24h_STm_0 - h24h_mock_0,
    # STm effect with RANKL at 24h
    STm_24h_100 = h24h_STm_100 - h24h_mock_100,
    # RANKL effect without infection at 24h
    RANKL_24h_mock = h24h_mock_100 - h24h_mock_0,
    # RANKL effect with infection at 24h
    RANKL_24h_STm = h24h_STm_100 - h24h_STm_0,
    # Desing matrix 
    levels = design
)

#################################
#### Pooled contrast analsis ####
#################################

## Fit the model
fit <- edgeR::voomLmFit(
    counts = DGErankl,
    design = design,
    sample.weights = TRUE
)
## Compute contrasts from the linear model fit
fit <- contrasts.fit(fit, contrasts.p)
## Compute Bayes moderated t-statistics and p-values
fit <- eBayes(fit, robust = TRUE)

##############
#### GSEA ####
# Run GSEA for all pool contrasts
gsea_results <- list()
contrasts_to_analyze <- colnames(contrasts.p)

for(contrast in contrasts_to_analyze) {
    # Get DE results for each contrast
    de_results <- topTable(fit, coef = contrast, 
    sort.by = "t", adjust.method = "fdr", n = Inf)
    
    # Clean up gene names and remove any empty entries
    de_results <- de_results[rownames(de_results) != "", ]
    
    # Create ranked gene list
    ranked_genes <- de_results$t
    names(ranked_genes) <- rownames(de_results)
    ranked_genes <- ranked_genes[!is.na(ranked_genes)]
    
    tryCatch({
        # Run GSEA for each database
        gsea_results[[contrast]][["hallmark"]] <- runGSEA(
            de_results,
            rank_metric = "t",
            species = "Mus musculus",
            category = "H",
            padj_method = "fdr",
            nperm = 100000,
            pvalue_cutoff = 0.01
        )
        
        gsea_results[[contrast]][["kegg"]] <- runGSEA(
            de_results,
            rank_metric = "t",
            species = "Mus musculus",
            category = "C2",
            subcategory = "CP:KEGG",
            padj_method = "fdr",
            nperm = 100000,
            pvalue_cutoff = 0.01
        )
        
        gsea_results[[contrast]][["gobp"]] <- runGSEA(
            de_results,
            rank_metric = "t",
            species = "Mus musculus",
            category = "C5",
            subcategory = "GO:BP",
            padj_method = "fdr",
            nperm = 100000,
            pvalue_cutoff = 0.01
        )
    }, error = function(e) {
        message(sprintf("Error in contrast %s: %s", contrast, e$message))
    })
}

# Pool significant pathways
get_significant_pathways <- function(gsea_results, q_cutoff = 0.01) {
    significant_paths <- lapply(gsea_results, function(x) {
        if(!is.null(x)) {
            x$ID[x$p.adjust < q_cutoff]
        }
    })
    unique(unlist(significant_paths))
}

# Get pools for each database
hallmark_pool <- get_significant_pathways(lapply(gsea_results, `[[`, "hallmark"))
kegg_pool <- get_significant_pathways(lapply(gsea_results, `[[`, "kegg"))
gobp_pool <- get_significant_pathways(lapply(gsea_results, `[[`, "gobp"))


# Extract gene sets from GSEA results
get_pathway_genes <- function(gsea_results) {
    # Extract results data frame
    results_df <- gsea_results@result
    
    # Get significant pathways
    sig_results <- results_df[results_df$p.adjust < 0.01, ]
    
    # Create list of gene sets
    pathway_genes <- strsplit(sig_results$core_enrichment, "/")
    names(pathway_genes) <- sig_results$ID
    
    return(pathway_genes)
}

# Get pathway genes for each database
hallmark_genes <- get_pathway_genes(gsea_results$Time_RANKL_effect$hallmark)
kegg_genes <- get_pathway_genes(gsea_results$Time_RANKL_effect$kegg)
gobp_genes <- get_pathway_genes(gsea_results$Time_RANKL_effect$gobp)

# Calculate normalized expression scores
norm_counts <- cpm(DGErankl, log = TRUE)

calculate_pathway_scores <- function(expression_data, pathway_genes, method = "mean") {
    scores <- matrix(NA, 
                    nrow = ncol(expression_data),
                    ncol = length(pathway_genes),
                    dimnames = list(colnames(expression_data), 
                                  names(pathway_genes)))
    
    for(pathway in names(pathway_genes)) {
        genes <- pathway_genes[[pathway]]
        # Get genes that are present in expression data
        genes <- genes[genes %in% rownames(expression_data)]
        pathway_exp <- expression_data[genes, ]
        scores[, pathway] <- colMeans(pathway_exp, na.rm = TRUE)
    }
    return(scores)
}

# calculate normalized expression scores for each pathway database
hallmark_scores <- calculate_pathway_scores(norm_counts, hallmark_genes)
kegg_scores <- calculate_pathway_scores(norm_counts, kegg_genes)
gobp_scores <- calculate_pathway_scores(norm_counts, gobp_genes)


# Modified heatmap function
plot_pathway_heatmap <- function(scores, title, annotation_col) {
    pheatmap(t(scores),  # transpose the matrix
             scale = "row",
             clustering_method = "complete",
             clustering_distance_rows = "correlation",
             clustering_distance_cols = "correlation",
             show_rownames = TRUE,
             show_colnames = FALSE,
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             color = colorRampPalette(c("navy", "white", "red"))(50),
             main = title)
}

# Add sample annotations
annotation_col <- data.frame(
    Time = DGErankl$samples$timepoint,
    Treatment = DGErankl$samples$treatment,
    RANKL = DGErankl$samples$rankl,
    row.names = colnames(norm_counts)
)

# Create heatmaps with annotations
# Define annotation colors
ann_colors <- list(
    timepoint = c("4h" = "#1f77b4", "24h" = "#ff7f0e"),
    treatment = c("mock" = "#2ca02c", "STm" = "#d62728"),
    rankl = c("0" = "#9467bd", "100" = "#8c564b")
)


# Create heatmaps with annotations
p1 <- plot_pathway_heatmap(hallmark_scores, "Hallmark Pathways", annotation_col)
p2 <- plot_pathway_heatmap(kegg_scores, "KEGG Pathways", annotation_col)
p3 <- plot_pathway_heatmap(gobp_scores, "GO Biological Processes", annotation_col)
# Save the heatmaps
pdf("/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/Pooled/hallmark_heatmap.pdf", 
    width = 12, 
    height = 8)
p1
dev.off()

pdf("/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/Pooled/kegg_heatmap.pdf", 
    width = 12, 
    height = 8)
p2
dev.off()

pdf("/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/Pooled/gobp_heatmap.pdf", 
    width = 16, 
    height = 25)
p3
dev.off()


######################################
#### Independent contrast analsis ####
######################################

# Statistical estimation
## Fit the model
fit <- edgeR::voomLmFit(
    counts = DGErankl,
    design = design,
    sample.weights = TRUE
)
## Compute contrasts from the linear model fit
fit <- contrasts.fit(fit, contrasts)
## Compute Bayes moderated t-statistics and p-values
fit <- eBayes(fit, robust = TRUE)
## Extract contrasts coefficients
### Infection response with RANKL at 4h
DE_rankl_infection_4h <- topTable(fit,
    coef = "STm_response_100RANKL_4h",
    sort.by = "t",
    adjust.method = "fdr",
    n = Inf
)

# Infection response with RANKL at 24h
DE_rankl_infection_24h <- topTable(fit,
    coef = "STm_response_100RANKL_24h",
    sort.by = "t",
    adjust.method = "fdr",
    n = Inf
)

### Early infection RANKL modulation
DE_rankl_4h <- topTable(fit,
    coef = "RANKL_effect_4h",
    sort.by = "t",
    adjust.method = "fdr",
    n = Inf
)
### Late infection RANKL modulation
DE_rankl_24h <- topTable(fit,
    coef = "RANKL_effect_24h",
    sort.by = "t",
    adjust.method = "fdr",
    n = Inf
)
### Time-dependent RANKL modulation of the infection effect
DE_rankl_time <- topTable(fit,
    coef = "Time_RANKL_effect",
    sort.by = "t",
    adjust.method = "fdr",
    n = Inf
)
### Save results
write.csv(DE_rankl_4h, "3_Results/DE_tables/DE_rankl_4h.csv")
write.csv(DE_rankl_24h, "3_Results/DE_tables/DE_rankl_24h.csv")
write.csv(DE_rankl_time, "3_Results/DE_tables/DE_rankl_time.csv")



# Volcano plots 
## Early infection
create_volcano_plot(
  DE_rankl_4h,
  p_cutoff = 0.05,
  fc_cutoff = 2.0,
  x_breaks = 2,
  max.overlaps = 10,
  label_method = "sig",
  title = "RANKL+infection response change at 4 hours"
)

## Late infection
create_volcano_plot(
  DE_rankl_24h,
  p_cutoff = 0.05,
  fc_cutoff = 2.0,
  x_breaks = 2,
  max.overlaps = 8,
  label_method = "sig",
  title = "RANKL+infection response change at 24 hours"
)

## Time-dependent RANKL effect
create_volcano_plot(
  DE_rankl_time,
  p_cutoff = 0.05,
  fc_cutoff = 2.0,
  x_breaks = 2,
  max.overlaps = 10,
  label_method = "sig",
  title = "RANKL+infection DEGs changing over time: 24h vs 4h "
)

## Combined volcano 
create_combined_volcano_plots(
    DE_rankl_4h,
    DE_rankl_24h,
    DE_rankl_time,
    p_cutoff = 0.05,
    fc_cutoff = 2.0
)


# GSEA 
## Time-dependent RANKL effect 
### 1. REACTOME
gsea_reactome_time <- runGSEA(
  DE_rankl_time,
  rank_metric = "t",
  species = "Mus musculus",
  category = "C2",
  subcategory = "CP:REACTOME",
  padj_method = "fdr",
  nperm = 100000,
  pvalue_cutoff = 0.05
)


### 2. HALLMARK
gsea_hallmark_time <- runGSEA(
  DE_rankl_time,
  rank_metric = "t",
  species = "Mus musculus",
  category = "H",
  padj_method = "fdr",
  nperm = 100000,
  pvalue_cutoff = 0.05
)

# GSEA visualisation
## 1. REACTOME
GSEA_dotplot(
  gsea_obj = gsea_reactome_time,
  showCategory = 15,
  font.size = 10,
  title = "Top 15 KO enriched REACTOME time-dependent gene sets",
  sortBy = "GeneRatio",
  filterBy = "NES",
  q_cut = 0.05,
  replace_ = TRUE,
  capitalize_1 = FALSE,
  capitalize_all = FALSE,
  min.dotSize = 2
)

runSumGSEAplot(
  gsea_reactome_time,
  gene_set_ids = 1:5,
  title = "REACTOME running sum for top 5 pathways"
)

## 2. HALLMARK
GSEA_dotplot(
  gsea_obj = gsea_hallmark_time,
  showCategory = 15,
  font.size = 10,
  title = "Top 15 KO enriched HALLMARK time-dependent gene sets",
  sortBy = "GeneRatio",
  filterBy = "NES",
  q_cut = 0.05,
  replace_ = TRUE,
  capitalize_1 = FALSE,
  capitalize_all = FALSE,
  min.dotSize = 2
)

runSumGSEAplot(
  gsea_hallmark_time,
  gene_set_ids = 1:5,
  title = "HALLMARK running sum for top 5 pathways"
)
