# RANKL Signaling in Macrophages

This repository contains the bulk RNA-seq analysis pipeline used to explore the effects of **RANKL signaling** on the **innate immune response** of macrophages to *Salmonella* infection. The analysis is performed in **R**, with a focus on both **differential expression** and **gene set enrichment analysis (GSEA)** to uncover pathway-level responses.

## 🧪 Experimental Design

All samples are murine macrophages treated with either **no RANKL (0)** or **100 ng/mL RANKL** for 48 hours prior to mock infection or *Salmonella Typhimurium* (STm) infection. Cells were harvested at either **4 hours** or **24 hours** post-infection. Each condition was run in **triplicate**.

```mermaid
graph TD
    A[Macrophages] --> B[No RANKL]
    A --> C[RANKL]
    
    B --> D[Mock Infection]
    B --> E[Salmonella]
    C --> F[Mock Infection]
    C --> G[Salmonella]
    
    D --> H[4h Harvest]
    D --> I[24h Harvest]
    E --> J[4h Harvest]
    E --> K[24h Harvest]
    F --> L[4h Harvest]
    F --> M[24h Harvest]
    G --> N[4h Harvest]
    G --> O[24h Harvest]
```

### 📊 Experimental Conditions (n = 8 groups, 3 replicates each)

- `4h_mock_0`  
- `24h_mock_0`  
- `4h_STm_0`  
- `24h_STm_0`  
- `4h_mock_100`  
- `24h_mock_100`  
- `4h_STm_100`  
- `24h_STm_100`  

---


## 🧬 Biological Objective

The main question is: **How does RANKL modulate the macrophage response to infection?**

Previous biological findings suggest RANKL **dampens pro-inflammatory signaling** (e.g., TLR/NF-κB pathways). This pipeline quantifies those effects at the transcriptional and pathway level, across timepoints.

---

## 🧠 Statistical Design & Contrasts

We define **specific contrasts** to isolate the biological effects of interest:
## Statistical Modelling 
To model the experiment in mathematical terms, I`m going to set the following contrasts: 

#### 🧮 Defining Contrasts

```r
contrasts.p <- makeContrasts(
    STm_4h_0 = t4h_STm_0 - t4h_mock_0,
    STm_4h_100 = t4h_STm_100 - t4h_mock_100,
    RANKL_4h_mock = t4h_mock_100 - t4h_mock_0,
    RANKL_4h_STm = t4h_STm_100 - t4h_STm_0,
    STm_24h_0 = t24h_STm_0 - t24h_mock_0,
    STm_24h_100 = t24h_STm_100 - t24h_mock_100,
    RANKL_24h_mock = t24h_mock_100 - t24h_mock_0,
    RANKL_24h_STm = t24h_STm_100 - t24h_STm_0,
    levels = design
)
```
These contrasts look at **baseline** RANKL and infection effects for each timepoint.

```r
contrasts.m <- makeContrasts(
    RANKL_effect_4h = (t4h_STm_100 - t4h_mock_100) - (t4h_STm_0 - t4h_mock_0),
    RANKL_effect_24h = (t24h_STm_100 - t24h_mock_100) - (t24h_STm_0 - t24h_mock_0),
    Time_RANKL_effect = ((t24h_STm_100 - t24h_mock_100) - (t24h_STm_0 - t24h_mock_0)) -
                        ((t4h_STm_100 - t4h_mock_100) - (t4h_STm_0 - t4h_mock_0)),
    levels = design
)
```
These **granular contrasts** isolate how **RANKL** modulates infection specifically at 4h and 24h, and how that change evolves over time.

Basically speaking, the `-` sign sets the comparison between two conditions. This helps us to isolate specific biological effects we`re interested in. 

### Early infection RANKL effect 
`RANKL_effect_4h = (4h_STm_100 - 4h_mock_100) - (4h_STm_0 - 4h_mock_0)`
* This shows how RANKL changes the infection response at 4 hours
* Think of it as: "Does RANKL treatment make cells respond differently to infection (Salmonella) at 4h?"
* In mathematical terms:
    - First, we calculate infection effect with RANKL (4h_STm_100 - 4h_mock_100)
    - Then, infection effect without RANKL (4h_STm_0 - 4h_mock_0)
    - The difference between these tells us if RANKL modified the infection response
* It compares the infection response (STm vs mock) between RANKL-treated and untreated cells at 4 hours
* Identifies genes where RANKL changes the early infection response

### Late infection RANKL effect 
`RANKL_effect_24h = (24h_STm_100 - 24h_mock_100) - (24h_STm_0 - 24h_mock_0)`
* Same as above, but at 24 hours
* Answers: "Does RANKL treatment make cells respond differently to infection at 24h?"
Uses the same mathematical logic, just with 24h samples
* It compares the infection response (STm vs mock) between RANKL-treated and untreated cells at 24 hours
* Identifies genes where RANKL changes the late infection response

### RANKL time effect 
`Time_RANKL_effect = ((24h_STm_100 - 24h_mock_100) - (24h_STm_0 - 24h_mock_0)) - ((4h_STm_100 - 4h_mock_100) - (4h_STm_0 - 4h_mock_0))`
* This shows how the RANKL effect on infection changes between 4h and 24h
* Answers: "Does RANKL's impact on infection response change over time?"
* Mathematically:
    - Calculate RANKL effect at 24h
    - Calculate RANKL effect at 4h
    - Subtract to see if the effect changes over time
* Identifies genes where RANKL's modification of the infection response differs between early and late timepoints
* Positive values indicate genes where RANKL's effect on infection response is stronger at 24h
* Negative values indicate genes where RANKL's effect on infection response is stronger at 4h



---

## 🧰 Pipeline Overview

### 1. **Environment Setup**
Scripts are modularized under `1_Scripts/`:
- `Load_libraries.R`: Loads required packages
- `process_rnaseq_data.R`: Creates `DGEList` and normalizes counts
- `PCA.R`: Generates PCA plots
- GSEA and visualization tools: `runGSEA.R`, `GSEA_dotplot.R`, `runningSumGSEAplot.R`, `combined_volcano.R`, etc.

---

### 2. **Data Preparation**
- Reads **count matrix** and **sample metadata**
- Constructs a `DGEList` for downstream analysis
- Runs **PCA** to assess sample clustering

---

### 3. **Pooled GSEA Analysis**

#### 🧮 Defining Contrasts
```r
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

# Granular contrasts
contrasts.m <- makeContrasts(
    RANKL_effect_4h = (t4h_STm_100 - t4h_mock_100) - (t4h_STm_0 - t4h_mock_0),
    RANKL_effect_24h = (t24h_STm_100 - t24h_mock_100) - (t24h_STm_0 - t24h_mock_0),
    Time_RANKL_effect = ((t24h_STm_100 - t24h_mock_100) - (t24h_STm_0 - t24h_mock_0)) - ((t4h_STm_100 - t4h_mock_100) - (t4h_STm_0 - t4h_mock_0)),
    levels = design
)
```
### 📈 Running Pooled GSEA

```r
pooled_gsea_results <- run_pooled_gsea(fit, contrasts.p, DGErankl)
```

**What does `run_pooled_gsea()` do?**  
- Iterates over each contrast (e.g., `RANKL_effect_4h`, `RANKL_effect_24h`), then:
  1. **Extracts differentially expressed genes** using `topTable()`.
  2. **Runs GSEA** on the ranked gene list (t-statistic) via `runGSEA()` using MSigDB references (Hallmark, KEGG, GO:BP, Reactome).
  3. Collects significant pathways from each contrast and **merges** them into a pooled set.
  4. Retrieves the corresponding **core genes** for these pathways using `get_pathway_genes_all()`.
  5. **Normalizes** counts (`norm_counts`) and **calculates pathway scores** by averaging gene expression for each pathway.

**Key outputs**:
- **`gsea_results`**: Full GSEA results (by contrast and database).
- **`pools`**: List of significantly enriched pathways across all contrasts.
- **`genes`**: Genes in each enriched pathway, per database.
- **`scores`**: Sample-by-pathway score matrix for heatmaps.

**Goal**:  
Perform GSEA across multiple contrasts and **aggregate** significant pathways and genes.

**Function chain**:  
1. `topTable()` → pulls DE genes per contrast  
2. `runGSEA()` → runs enrichment with `clusterProfiler` & MSigDB  
3. `get_significant_pathways()` → extracts top pathways across contrasts  
4. `get_pathway_genes_all()` → collects core genes for top pathways  
5. `calculate_pathway_scores()` → computes average expression per pathway

---

### 🧠 Supporting Functions Breakdown

- **`run_pooled_gsea()`**  
  - Loops over contrasts → runs GSEA → retrieves genes → computes scores  
  - Returns GSEA results, significant pathways, and pathway scores  

- **`runGSEA()`**  
  - Runs GSEA using `clusterProfiler` on a ranked list (t-statistic)  
  - Supports multiple MSigDB categories (Hallmark, KEGG, GO:BP, Reactome)

- **`get_pathway_genes()`**  
  - Extracts core enrichment genes from top significant pathways (by adjusted p-value)  

- **`get_pathway_genes_all()`**  
  - Aggregates pathway genes across all contrasts, ranking by minimum adjusted p-value  

- **`get_significant_pathways()`**  
  - Returns all pathways passing a specified FDR cutoff

- **`calculate_pathway_scores()`**  
  - For each pathway, selects its genes and averages their expression per sample  
  - Outputs a “samples × pathways” score matrix

---

## 🧩 Key Functions Explained

| Function                          | Description                                                                      |
|----------------------------------|----------------------------------------------------------------------------------|
| **`runGSEA()`**                  | Runs GSEA using a ranked gene list (e.g., t-statistics) and MSigDB categories.   |
| **`run_pooled_gsea()`**          | Applies `runGSEA()` to multiple contrasts and compiles pooled results.           |
| **`get_pathway_genes()`**        | Retrieves top significant pathways and their core enrichment genes.              |
| **`calculate_pathway_scores()`** | Computes average expression of pathway genes per sample.                         |
| **`create_volcano_plot()`**      | Generates volcano plots with labeling and threshold options.                     |
| **`plot_pathway_heatmap()`**     | Plots expression heatmaps with metadata annotations.                             |
| **`create_combined_volcano_plots()`** | Merges multiple volcano plots into a single visualization.               |

---

## 🖼️ Output Files

All plots and tables are saved under:

```
3_Results/
  ├── DE_tables/                   # Differential expression results
  ├── imgs/
      ├── GSEA/
          ├── A.Pair-wise comparisons/  # Heatmaps from pooled GSEA
          ├── B.Granular_questions/     # Dotplots, heatmaps, and GSEA enrichment plots for specific contrasts
```

---
