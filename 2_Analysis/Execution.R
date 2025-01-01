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
colnames(design) <- paste0("h", colnames(design))
## Set contrast matrix
### My contrast Desing
contrasts <- makeContrasts(
    # RANKL effect in mock conditions at 4h
    RANKL_effect_mock_4h = h4h_mock_100 - h4h_mock_0,
    # RANKL effect in mock conditions at 24h
    RANKL_effect_mock_24h = h24h_mock_100 - h24h_mock_0,
    # RANKL effect in STm infection at 4h
    RANKL_effect_STm_4h = h4h_STm_100 - h4h_STm_0,
    # RANKL effect in STm infection at 24h
    RANKL_effect_STm_24h = h24h_STm_100 - h24h_STm_0,
    # Infection response without RANKL at 4h
    STm_response_0RANKL_4h = h4h_STm_0 - h4h_mock_0,
    # Infection response with RANKL at 4h
    STm_response_100RANKL_4h = h4h_STm_100 - h4h_mock_100,
    # Infection response without RANKL at 24h
    STm_response_0RANKL_24h = h24h_STm_0 - h24h_mock_0,
    # Infection response with RANKL at 24h
    STm_response_100RANKL_24h = h24h_STm_100 - h24h_mock_100,
    # RANKL effect on infection response at 4h
    RANKL_effect_4h = (h4h_STm_100 - h4h_mock_100) - (h4h_STm_0 - h4h_mock_0),
    # RANKL effect on infection response at 24h
    RANKL_effect_24h = (h24h_STm_100 - h24h_mock_100) - (h24h_STm_0 - h24h_mock_0),
    # RANKL time-dependent effect 
    Time_RANKL_effect = ((h24h_STm_100 - h24h_mock_100) - (h24h_STm_0 - h24h_mock_0)) - ((h4h_STm_100 - h4h_mock_100) - (h4h_STm_0 - h4h_mock_0)),
    # Desing matrix
    levels = design
)
### Denis contrast desing 
contrasts.d <- makeContrasts(
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
    # Time effect in mock without RANKL
    Time_mock_0 = h24h_mock_0 - h4h_mock_0,
    # Time effect in mock with RANKL
    Time_mock_100 = h24h_mock_100 - h4h_mock_100,
    # Time effect in STm without RANKL
    Time_STm_0 = h24h_STm_0 - h4h_STm_0,
    # Time effect in STm with RANKL
    Time_STm_100 = h24h_STm_100 - h4h_STm_100,
    # Desing matrix 
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
