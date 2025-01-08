##############################
### 1. Prepare environment ###
#############################

# Load libraries and functions
source("1_Scripts/Load_libraries.R")
source("1_Scripts/process_rnaseq_data.R")
source("1_Scripts/PCA.R")
source("1_Scripts/standard_volcano.R")
source("1_Scripts/GSEA_dotplot.R")
source("1_Scripts/runGSEA.R")
source("1_Scripts/runningSumGSEAplot.R")
source("1_Scripts/combined_volcano.R")
source("1_Scripts/runGSEA_pool.R")

#####################
### 2. Read data ###
#####################

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


#####################
### 3. Pool GSEA ###
#####################
# Set design matrix
## Set combinatorial design
design <- model.matrix(~ 0 + group, data = DGErankl$samples)
## Rename the columns to remove "group" prefix
colnames(design) <- gsub("group", "", colnames(design))

## Set contrast matrix
### Pool contrast design
contrasts.p <- makeContrasts(
    # STm effect without RANKL at 4h
    STm_4h_0 = t4h_STm_0 - t4h_mock_0,
    # STm effect with RANKL at 4h
    STm_4h_100 = t4h_STm_100 - t4h_mock_100,
    # RANKL effect without infection at 4h
    RANKL_4h_mock = t4h_mock_100 - t4h_mock_0,
    # RANKL effect with infection at 4h
    RANKL_4h_STm = t4h_STm_100 - t4h_STm_0,
    # STm effect without RANKL at 24h
    STm_24h_0 = t24h_STm_0 - t24h_mock_0,
    # STm effect with RANKL at 24h
    STm_24h_100 = t24h_STm_100 - t24h_mock_100,
    # RANKL effect without infection at 24h
    RANKL_24h_mock = t24h_mock_100 - t24h_mock_0,
    # RANKL effect with infection at 24h
    RANKL_24h_STm = t24h_STm_100 - t24h_STm_0,
    # Design
    levels = design
)

### Select contrast desing
contrasts.m <- makeContrasts(
    # RANKL effect on infection response at 4h
    RANKL_effect_4h = (t4h_STm_100 - t4h_mock_100) - (t4h_STm_0 - t4h_mock_0),
    # RANKL effect on infection response at 24h
    RANKL_effect_24h = (t24h_STm_100 - t24h_mock_100) - (t24h_STm_0 - t24h_mock_0),
    # RANKL time-dependent effect
    Time_RANKL_effect = ((t24h_STm_100 - t24h_mock_100) - (t24h_STm_0 - t24h_mock_0)) - ((t4h_STm_100 - t4h_mock_100) - (t4h_STm_0 - t4h_mock_0)),
    # Design
    levels = design
)

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
# Run pooled GSEA analysis
pooled_gsea_results <- run_pooled_gsea(fit, contrasts.p, DGErankl)


# Create ordered sample vector
sample_order <- with(DGErankl$samples, {
    # Create all combinations in desired order
    rankl_levels <- c("0", "100")
    treatment_levels <- c("mock", "STm")
    time_levels <- c("4h", "24h")
    
    # Generate ordered vector
    ordered_samples <- character()
    for(r in rankl_levels) {
        for(tr in treatment_levels) {
            for(t in time_levels) {
                ordered_samples <- c(ordered_samples, 
                    rownames(DGErankl$samples)[DGErankl$samples$rankl == r & 
                                             DGErankl$samples$treatment == tr & 
                                             DGErankl$samples$timepoint == t])
            }
        }
    }
    ordered_samples
})

# Modified heatmap function with controlled column order
plot_pathway_heatmap <- function(scores, title, annotation_col) {
    pheatmap(t(scores),  
             scale = "row",
             clustering_method = "complete",
             clustering_distance_rows = "correlation",
             cluster_cols = FALSE,  # Disable column clustering
             show_rownames = TRUE,
             show_colnames = FALSE,
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             color = colorRampPalette(c("navy", "white", "red"))(50),
             main = title,
             gaps_col = c(12, 24),  # Add visual separators between main groups
             colnames = sample_order)  # Set column order
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

# Create and save heatmaps
pdf("/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/Pooled/hallmark_heatmap.pdf", 
    width = 12, 
    height = 8)
plot_pathway_heatmap(pooled_gsea_results$scores$hallmark, "Hallmark Pathways", annotation_col)
dev.off()

pdf("/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/Pooled/kegg_heatmap.pdf", 
    width = 12, 
    height = 8)
plot_pathway_heatmap(pooled_gsea_results$scores$kegg, "KEGG Pathways", annotation_col)
dev.off()

pdf("/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/Pooled/gobp_heatmap.pdf", 
    width = 12, 
    height = 8)
plot_pathway_heatmap(pooled_gsea_results$scores$gobp, "GO Biological Processes", annotation_col)
dev.off()

pdf("/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/Pooled/reactome_heatmap.pdf", 
    width = 18, 
    height = 8)
plot_pathway_heatmap(pooled_gsea_results$scores$reactome, "REACTOME Pathways", annotation_col)
dev.off()

###################################
### 4. Specific-questions GSEA ###
###################################

### Select contrast desing
contrasts.m <- makeContrasts(
    # RANKL effect on infection response at 4h
    RANKL_effect_4h = (t4h_STm_100 - t4h_mock_100) - (t4h_STm_0 - t4h_mock_0),
    # RANKL effect on infection response at 24h
    RANKL_effect_24h = (t24h_STm_100 - t24h_mock_100) - (t24h_STm_0 - t24h_mock_0),
    # RANKL time-dependent effect
    Time_RANKL_effect = ((t24h_STm_100 - t24h_mock_100) - (t24h_STm_0 - t24h_mock_0)) - ((t4h_STm_100 - t4h_mock_100) - (t4h_STm_0 - t4h_mock_0)),
    # Design
    levels = design
)

# Statistical estimation
## Fit the model
fit <- edgeR::voomLmFit(
    counts = DGErankl,
    design = design,
    sample.weights = TRUE
)
## Compute contrasts from the linear model fit
fit <- contrasts.fit(fit, contrasts.m)
## Compute Bayes moderated t-statistics and p-values
fit <- eBayes(fit, robust = TRUE)
## Extract contrasts coefficients
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
p1 <- create_volcano_plot(
  DE_rankl_4h,
  p_cutoff = 0.05,
  fc_cutoff = 2.0,
  x_breaks = 2,
  max.overlaps = 10,
  label_method = "sig",
  title = "RANKL+infection response change at 4 hours"
)

## Late infection
p2 <- create_volcano_plot(
  DE_rankl_24h,
  p_cutoff = 0.05,
  fc_cutoff = 2.0,
  x_breaks = 2,
  max.overlaps = 8,
  label_method = "sig",
  title = "RANKL+infection response change at 24 hours"
)

## Time-dependent RANKL effect
p3 <- create_volcano_plot(
  DE_rankl_time,
  p_cutoff = 0.05,
  fc_cutoff = 2.0,
  x_breaks = 2,
  max.overlaps = 10,
  label_method = "sig",
  title = "RANKL+infection DEGs changing over time: 24h vs 4h "
)

## Combined volcano 
p4 <- create_combined_volcano_plots(
    DE_rankl_4h,
    DE_rankl_24h,
    DE_rankl_time,
    p_cutoff = 0.05,
    fc_cutoff = 2.0,
    max.overlaps = 15
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

### 3. GO:BP
gsea_gobp_time <- runGSEA(
  DE_rankl_time,
  rank_metric = "t",
  species = "Mus musculus",
  category = "C5",
  subcategory = "GO:BP",
  padj_method = "fdr",
  nperm = 100000,
  pvalue_cutoff = 0.05
)

### 4. KEGG
gsea_kegg_time <- runGSEA(
  DE_rankl_time,
  rank_metric = "t",
  species = "Mus musculus",
  category = "C2",
  subcategory = "CP:KEGG",
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

## 3. GO:BP
GSEA_dotplot(
  gsea_obj = gsea_gobp_time,
  showCategory = 15,
  font.size = 10,
  title = "Top 15 KO enriched GOBP time-dependent gene sets",
  sortBy = "GeneRatio",
  filterBy = "NES",
  q_cut = 0.05,
  replace_ = TRUE,
  capitalize_1 = FALSE,
  capitalize_all = FALSE,
  min.dotSize = 2
)

runSumGSEAplot(
  gsea_gobp_time,
  gene_set_ids = 1:5,
  title = "GOBP running sum for top 5 pathways"
)

## 4. KEGG 
GSEA_dotplot(
  gsea_obj = gsea_kegg_time,
  showCategory = 15,
  font.size = 10,
  title = "Top 15 KO enriched KEGG time-dependent gene sets",
  sortBy = "GeneRatio",
  filterBy = "NES",
  q_cut = 0.05,
  replace_ = TRUE,
  capitalize_1 = FALSE,
  capitalize_all = FALSE,
  min.dotSize = 2
)

runSumGSEAplot(
  gsea_kegg_time,
  gene_set_ids = 1:5,
  title = "KEGG running sum for top 5 pathways"
)

# Get pathway genes for each database
hallmark_genes <- get_pathway_genes(gsea_hallmark_time)
kegg_genes <- get_pathway_genes(gsea_kegg_time)
gobp_genes <- get_pathway_genes(gsea_gobp_time)
reactome_genes <- get_pathway_genes(gsea_reactome_time)

# Calculate normalized expression scores
norm_counts <- cpm(DGErankl, log = TRUE)

# Calculate pathway scores
hallmark_scores <- calculate_pathway_scores(norm_counts, hallmark_genes)
kegg_scores <- calculate_pathway_scores(norm_counts, kegg_genes)
gobp_scores <- calculate_pathway_scores(norm_counts, gobp_genes)
reactome_scores <- calculate_pathway_scores(norm_counts, reactome_genes)

# Create and save heatmaps
pdf("3_Results/imgs/GSEA/Time_dependent/hallmark_heatmap.pdf", width = 12, height = 8)
plot_pathway_heatmap(hallmark_scores, "Hallmark Pathways", annotation_col)
dev.off()

pdf("3_Results/imgs/GSEA/Time_dependent/kegg_heatmap.pdf", width = 12, height = 8)
plot_pathway_heatmap(kegg_scores, "KEGG Pathways", annotation_col)
dev.off()

pdf("3_Results/imgs/GSEA/Time_dependent/gobp_heatmap.pdf", width = 12, height = 8)
plot_pathway_heatmap(gobp_scores, "GO Biological Processes", annotation_col)
dev.off()

pdf("3_Results/imgs/GSEA/Time_dependent/reactome_heatmap.pdf", width = 20, height = 8)
plot_pathway_heatmap(reactome_scores, "REACTOME Pathways", annotation_col)
dev.off()


## tease apart which pathways RANKL is downregulating (or upregulating) in STm infected cells? i.e. RANKL effect on infection, 4h or 24h individually  

# Define databases to analyze
db_configs <- list(
    HALLMARK = list(category = "H", subcategory = NULL),
    REACTOME = list(category = "C2", subcategory = "CP:REACTOME"),
    KEGG = list(category = "C2", subcategory = "CP:KEGG"),
    GOBP = list(category = "C5", subcategory = "GO:BP")
)

# Function to run GSEA and create plots for one DE table
run_gsea_analysis <- function(de_table, timepoint) {
    # Define database-specific plot parameters
    db_plot_params <- list(
        HALLMARK = list(width = 10, height = 7, font.size = 10),
        REACTOME = list(width = 20, height = 11, font.size = 8),
        KEGG = list(width = 12, height = 8, font.size = 9),
        GOBP = list(width = 13, height = 8, font.size = 9)
    )
    
    for (db_name in names(db_configs)) {
        # Get database-specific parameters
        params <- db_plot_params[[db_name]]
        
        gsea_result <- runGSEA(
            DE_results = de_table,
            rank_metric = "t",
            category = db_configs[[db_name]]$category,
            subcategory = db_configs[[db_name]]$subcategory,
            padj_method = "fdr",
            nperm = 100000,
            pvalue_cutoff = 0.05
        )
        
        GSEA_dotplot(
            gsea_result,
            filterBy = "NES_positive",
            sortBy = "GeneRatio",
            font.size = params$font.size,
            showCategory = 15,
            q_cut = 0.05,
            replace_ = TRUE,
            capitalize_1 = FALSE,
            capitalize_all = FALSE,
            min.dotSize = 2,
            title = paste0(timepoint, " ", db_name, " Upregulated")
        )
        
        GSEA_dotplot(
            gsea_result,
            filterBy = "NES_negative",
            sortBy = "GeneRatio",
            font.size = params$font.size,
            showCategory = 15,
            q_cut = 0.05,
            replace_ = TRUE,
            capitalize_1 = FALSE,
            capitalize_all = FALSE,
            min.dotSize = 2,
            title = paste0(timepoint, " ", db_name, " Downregulated")
        )
        
        GSEA_barplot(
            gsea_result,
            title = paste0(timepoint, " ", db_name, " NES"),
            width = params$width,
            height = params$height
        )
    }
}

# Run analysis for both timepoints
run_gsea_analysis(DE_rankl_4h, "4h")
run_gsea_analysis(DE_rankl_24h, "24h")