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


# Set design matrix
## Set combinatorial design
design <- model.matrix(~ 0 + group, data = DGErankl$samples)
## Rename the columns to remove "group" prefix
colnames(design) <- gsub("group", "", colnames(design))
## Set contrast matrix
contrasts <- makeContrasts(
    RANKL_effect_4h = (`4h_STm_100` - `4h_mock_100`) - (`4h_STm_0` - `4h_mock_0`),
    RANKL_effect_24h = (`24h_STm_100` - `24h_mock_100`) - (`24h_STm_0` - `24h_mock_0`),
    Time_RANKL_effect = ((`24h_STm_100` - `24h_mock_100`) - (`24h_STm_0` - `24h_mock_0`)) -
        ((`4h_STm_100` - `4h_mock_100`) - (`4h_STm_0` - `4h_mock_0`)),
    levels = design
)

# PCA
create_pca_plot(DGErankl, title = "PCA Plot")


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
### Early infection
DE_rankl_4h <- topTable(fit,
    coef = "RANKL_effect_4h",
    sort.by = "t",
    adjust.method = "fdr",
    n = Inf
)
### Late infection
DE_rankl_24h <- topTable(fit,
    coef = "RANKL_effect_24h",
    sort.by = "t",
    adjust.method = "fdr",
    n = Inf
)
### Time-dependent RANKL effect
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
