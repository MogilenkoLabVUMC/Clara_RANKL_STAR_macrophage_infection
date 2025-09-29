## RANKL–Macrophage bulk RNA‑seq analysis: Methods and code map

This document serves both as a human‑readable Methods section and as a technical map of how the repository’s R scripts interoperate to produce the results. It explains the analysis philosophy, the statistical modeling, the GSEA pooling logic, and the plotting metrics. ASCII schematics show end‑to‑end data flow.

---

## Biological objective

We ask how pre‑exposure to RANKL modulates macrophage transcriptional responses to Salmonella infection and how those responses differ between 4 h and 24 h post‑infection.

Samples: murine macrophages with/without RANKL (0 vs 100 ng/mL, 48 h pre‑treatment), mock vs Salmonella (STm), harvested at 4 h or 24 h; n = 3 replicates per group.

---

## High‑level pipeline

```
┌─────────────────────────┐    ┌──────────────────────┐     ┌───────────────────────────┐
│ Raw counts (genes × N) │    │ Sample metadata      │     │ R packages                 │
│ 0_Data/STAR_Data/      │    │ seq_reference.txt    │     │ edgeR, limma, clusterProf.│
└─────────────┬──────────┘    └──────────┬───────────┘     └──────────┬────────────────┘
              │                           │                            │
              ▼                           ▼                            ▼
       process_rnaseq_data()  →  DGEList (filtered + TMM‑normalized; factors: timepoint, treatment, rankl, group)
              │
              ▼
   Design (~ 0 + group) + Contrasts  →  voomLmFit(sample.weights=TRUE) → contrasts.fit → eBayes(robust=TRUE)
              │
              ├─────────────────────────────────────────────────────────────────────────────────┐
              ▼                                                                                 ▼
      "Pooled" GSEA branch                                                           "Granular questions" branch
   (all contrasts; pool pathways)                                                  (specific contrasts/timepoints)
              │                                                                                 │
              ▼                                                                                 ▼
   topTable per contrast (t‑stat ranking)                                            DE tables per contrast
              │                                                                                 │
              ▼                                                                                 ▼
   runGSEA() over multiple MSigDB databases                                           runGSEA() per database
              │                                                                                 │
              ▼                                                                                 ▼
   get_significant_pathways() across contrasts                                         GSEA_dotplot (NES±),
   get_pathway_genes_all() (core genes by p.adjust min)                                GSEA_barplot (NES),
   calculate_pathway_scores() (mean log‑CPM)                                           GSEA_dotplot_facet (NES±)
              │                                                                                 │
              ▼                                                                                 ▼
   Pathway‑level heatmaps (scores)                                                    Optional gene‑level heatmaps
   (pheatmap; annotations; fixed column order)                                        for specific pathways
```

---

## Data preparation and normalization

- Input files: `0_Data/STAR_Data/counts.txt` (raw counts) and `0_Data/STAR_Data/seq_reference.txt` (metadata).
- Function: `process_rnaseq_data()` constructs an edgeR `DGEList`, filters lowly expressed genes (`filterByExpr`), and normalizes library sizes by TMM (`calcNormFactors`). It adds factors: `timepoint` (4/24), `treatment` (mock/STm), `rankl` (0/100), and a `group` factor combining these levels (e.g., `t4h_STm_100`). Optional Ensembl→gene symbol annotation is performed via `biomaRt`.

Key operations performed by `process_rnaseq_data`:

- Natural sort of columns by sample index embedded in column names
- Filter genes: `filterByExpr(DGE)`
- TMM normalization: `calcNormFactors(DGE, method = "TMM")`
- Factor construction for downstream modeling

---

## Statistical modeling and contrasts

- Design: no‑intercept model with one column per `group`: `design <- model.matrix(~ 0 + group, data = DGErankl$samples)`.
- Fitting: `edgeR::voomLmFit(..., sample.weights = TRUE)` followed by `limma::contrasts.fit` and `limma::eBayes(robust = TRUE)`.


Two sets of contrasts are used:

1) Pool contrasts (effects within each timepoint/treatment combination)

```49:68:2_Analysis/Analysis.R
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

2) Granular “modulation” contrasts (RANKL’s effect on infection at a timepoint and its time‑dependence)

```71:80:2_Analysis/Analysis.R
contrasts.m <- makeContrasts(
    RANKL_effect_4h = (t4h_STm_100 - t4h_mock_100) - (t4h_STm_0 - t4h_mock_0),
    RANKL_effect_24h = (t24h_STm_100 - t24h_mock_100) - (t24h_STm_0 - t24h_mock_0),
    Time_RANKL_effect = ((t24h_STm_100 - t24h_mock_100) - (t24h_STm_0 - t24h_mock_0)) -
                        ((t4h_STm_100 - t4h_mock_100) - (t4h_STm_0 - t4h_mock_0)),
    levels = design
)
```

Interpretation:

- `RANKL_effect_4h` and `RANKL_effect_24h` quantify how RANKL modifies infection‑induced changes at 4 h or 24 h respectively (difference‑of‑differences).
- `Time_RANKL_effect` captures how the RANKL modulation of infection differs between 24 h and 4 h.

Rationale for not using a single unified “time × treatment × RANKL” model as the primary inferential frame:

- PCA and heatmaps showed 24 h samples cluster distinctly in a way consistent with medium deprivation stress, a design‑coupled artifact not orthogonal to biological time. Because all 24 h samples share that stress, adding “time” as a single factor does not isolate this confounder; it risks attributing starvation‑driven shifts to temporal biology. We therefore treated 4 h and 24 h comparisons symmetrically via pre‑specified contrasts and emphasized interpretation within each timepoint. In the code I have left the explicit time‑interaction contrast (`Time_RANKL_effect`) but caution that it aggregates signal that likely includes medium‑deprivation effects.

---

## GSEA strategies and philosophy

We deploy two complementary GSEA tracks:

- Pooled GSEA (broad landscape): run GSEA for every contrast in `contrasts.p` across multiple databases; pool significant pathways across contrasts; extract core genes per pooled pathway; compute pathway scores (per sample) as average expression of core genes. This yields compact heatmaps capturing recurring pathway programs.

- Granular GSEA (focused questions): for targeted contrasts (`RANKL_effect_4h`, `RANKL_effect_24h`, `Time_RANKL_effect`), rank by moderated t, run GSEA per database, and visualize directionality (NES sign), magnitude (NES), significance (q‑value), and gene coverage (GeneRatio).

Why pool?

- Pooling stabilizes pathway selection against single‑contrast idiosyncrasies, highlighting pathways that recur or show strong evidence in at least one contrast. It supports downstream summarization (pathway scores) and interpretable heatmaps spanning all samples.

---

## Script interoperation (who calls whom and what flows where)

```
topTable() per contrast (moderated t)
   │  (rank by t)
   ▼
runGSEA()  [msigdbr gene sets; fgsea backend; nPermSimple = 100000]
   │  outputs gseaResult (with ID, Description, NES, p.adjust/qvalue, core_enrichment, setSize)
   ├─────────────────────────────┬───────────────────────────────────────────────────────────────────────┐
   ▼                             ▼                                                                       ▼
runGSEA_pool.R             get_significant_pathways.R                                                run_GSEA_analysis.R
   │                       (collect IDs with q < 0.01 across contrasts)                                  │
   │                                                                                                      │
   ▼                                                                                                      ▼
get_pathway_genes_all.R                                                                                    GSEA_dotplot.R
   ├─ parse core_enrichment to gene lists per pathway (p.adjust < 0.01)                                    • X = GeneRatio (= #core genes / setSize)
   ├─ deduplicate by pathway ID across contrasts                                                           • Y = pathway Description
   └─ rank by min p.adjust across contrasts; take top N                                                    • Color = NES sign (±), Size = -log10(q)
   │                                                                                                      • Filters: NES_positive / NES_negative / qvalue
   ▼                                                                                                      
calculate_pathway_scores.R                                                                                 GSEA_barplot.R
   • For each pooled pathway, subset expression (logCPM),                                                  • Bars = NES (signed)
     compute per‑sample mean across genes                                                                  • Order by NES (↑ or ↓), color by NES>0
   │                                                                                                      
   ▼                                                                                                      GSEA_dotplot_facet.R
Heatmaps of pathway scores (pheatmap)                                                                      • Facets by direction; X = |NES|; Size = -log10(q)
```

Key files and functions:

- `runGSEA.R`: prepares ranked vector (names = gene symbols; values = moderated t), collects MSigDB sets with `msigdbr(species, collection=category, subcollection=subcategory)`, runs `clusterProfiler::GSEA` with `by = "fgsea"`, `pAdjustMethod = "fdr"`, `nPermSimple = 100000`.
- `runGSEA_pool.R`: loops all pool contrasts, runs `runGSEA` for Hallmark, KEGG, GO:BP, and Reactome (pvalueCutoff 0.01), aggregates significant pathway IDs (`get_significant_pathways`), extracts pathway core genes across contrasts (`get_pathway_genes_all`), computes per‑sample pathway scores (`calculate_pathway_scores`) from `cpm(DGErankl, log=TRUE)`.
- `get_significant_pathways.R`: from a list of `gseaResult` objects, return unique pathway IDs with `p.adjust < q_cutoff` (default 0.01).
- `get_pathway_genes.R`: for a single `gseaResult`, select top N pathways by adjusted p (default q < 0.01), parse `core_enrichment` to gene symbols.
- `get_pathway_genes_all.R`: across contrasts, collect core genes for each pathway with `p.adjust < 0.01`, deduplicate by pathway ID, compute min adjusted p across contrasts per pathway, and optionally keep top N.
- `calculate_pathway_scores.R`: for each pathway gene set, intersect with expression matrix rows and compute the column‑wise mean (per sample). Output: samples × pathways matrix of average log‑CPM.
- `run_GSEA_analysis.R`: “granular” driver. Iterates a configured set of MSigDB categories (`db_configs`) and, for each: runs `runGSEA` on the DE table; emits three plots (dotplot positive, dotplot negative, NES barplot) and one faceted dotplot; returns a list of `gseaResult` objects by database. Parameters mirror those used in `Analysis.R` (rank by t, `nperm = 100000`, `q_cut = 0.05`).

---

## Plotting functions and what they visualize

- PCA (`PCA.R::create_pca_plot`): 2D PCA of samples over logCPM; point color = `group`; labels = group; axis labels display % variance explained. Saved to `3_Results/imgs/PCA/`.

- Volcano plots
  - `standard_volcano.R::create_volcano_plot`
    - X = log2 fold change (`logFC`), Y = −log10 p‑value (`P.Value`).
    - Thresholds: dashed lines at p‑value cutoff (default 0.05) and |log2FC| cutoff (default 2.0).
    - Color encodes significance categories: p‑value only, FC only, both, or NS; labels controlled by `label_method` and `max.overlaps`.
    - Used in `Analysis.R` for: 4 h, 24 h, and time‑interaction; saved to `3_Results/imgs/volcano/`.
  - `combined_volcano.R::create_combined_volcano_plots`
    - Builds 3 volcano plots (4 h, 24 h, 24 vs 4) with common axis limits for fair visual comparison; outputs a single PDF (landscape panel).

- GSEA dotplot (`GSEA_dotplot.R::GSEA_dotplot`)
  - Input: `gseaResult`.
  - Computes `GeneRatio = # core_enrichment genes / setSize` per term.
  - Filters by `qvalue < q_cut` (e.g., 0.05 in `Analysis.R`).
  - Encodings: X = `GeneRatio`; Y = term; dot size = −log10(qvalue); dot color = NES sign (positive/negative).
  - Supports `filterBy = "NES_positive"` and `"NES_negative"` to show up‑ or down‑regulated pathways separately (as used in `run_GSEA_analysis.R`).

- GSEA faceted dotplot (`GSEA_dotplot_facet.R::GSEA_faceted_dotplot`)
  - Splits significant pathways into Upregulated (NES > 0) and Downregulated (NES < 0) facets.
  - Encodings: X = |NES|; Y = term; dot size = −log10(qvalue).
  - Sorting configurable by `GeneRatio`, `qvalue`, or `NES`; `showCategory` limits items per facet.

- GSEA NES barplot (`GSEA_barplot.R::GSEA_barplot`)
  - Selects up to `top_n` significant pathways (q < `q_cut`) ordered by |NES|.
  - Encodings: Y (horizontal axis after flip) = term; bar length = NES; bar color indicates NES > 0 (up) vs NES < 0 (down).
  - Title and database‑specific figure sizes set by caller; used in `run_GSEA_analysis.R` with `top_n = 3 × showCategory`.

- GSEA running‑sum plot (`runningSumGSEAplot.R::runSumGSEAplot`)
  - Wrapper over `enrichplot::gseaplot2` to show the enrichment running‑sum curve, ranked list heatmap, and tick marks for selected gene sets.
  - Typical usage in `Analysis.R`: `subplots = c(1,2,3)` for composite view; saved to `3_Results/imgs/RunSum/`.

- Pathway heatmaps
  - Pooled pathway scores: in `Analysis.R`, `plot_pathway_heatmap(scores, ...)` (inline helper) draws heatmaps of transposed score matrices (`t(scores)`), `scale = "row"`, fixed sample order, column annotations (`Time`, `Treatment`, `RANKL`), optional visual gaps to separate groups. Color scale: blue‑white‑red.
  - Single‑pathway gene heatmaps: `plot_single_pathway_heatmap.R::plot_single_pathway_heatmap`
    - Finds a term in a `gseaResult` by `Description`/`ID` (partial match), parses `core_enrichment` genes, subsets the normalized expression matrix, and plots a gene‑level heatmap with optional annotations and sample order. Filename includes timepoint and pathway ID.

Parameterization (as called from `2_Analysis/Analysis.R`):

- Ranking metric for GSEA: moderated t‑statistic (`rank_metric = "t"`).
- GSEA engine: `clusterProfiler::GSEA(..., by="fgsea", nPermSimple=100000, pAdjustMethod="fdr")`.
- Significance: pooled GSEA uses `p.adjust < 0.01` for pooling; granular plots use `q_cut = 0.05`.
- Databases: Hallmark (H), KEGG (C2:CP:KEGG / KEGG_LEGACY), Reactome (C2:CP:REACTOME), GO (C5:BP/CC/MF), and additional collections (e.g., C3:TFT). Selection is specified in `db_configs` where used.

---

## Outputs

- DE tables: `3_Results/DE_tables/DE_rankl_4h.csv`, `DE_rankl_24h.csv`, `DE_rankl_time.csv`.
- Volcano plots: `3_Results/imgs/volcano/*.pdf` (single and combined triptych).
- Pooled GSEA heatmaps: `3_Results/imgs/GSEA/A.Pair‑wise comparisons/*_heatmap.pdf`.
- Granular GSEA plots (per DB/timepoint): `3_Results/imgs/GSEA/B.Granular_questions/` and `C.4_24h_RANKL_separate/`.
- Single‑pathway gene heatmaps: `3_Results/imgs/GSEA/.../Pathways/*.pdf`.
