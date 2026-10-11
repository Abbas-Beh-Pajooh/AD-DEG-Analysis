# AD-DEG-Analysis

RNA-seq differential gene expression analysis for Alzheimer's disease using DESeq2 in R.

## Overview

This repository contains an R script that compares Alzheimer's disease and control brain samples using RNA-seq gene counts. It uses DESeq2 to identify differentially expressed genes (DEGs) while adjusting for covariates. Finally, important plots will be generated.

## Methods

- **Model:** `~ age_at_death + rin_score + post_mortem_interval + sex + diagnosis` (Wald test, parametric dispersion fit)
- **Contrast:** Alzheimer vs. control (control is the reference level)
- **Covariates:** age at death, RNA integrity number (RIN), post-mortem interval (PMI) and sex. The numeric covariates are mean-imputed where missing and then scaled.
- **DEG thresholds:** adjusted p-value (FDR) < 0.05 and |log2 fold change| ≥ 1
- **Heatmap:** variance-stabilizing transformation (VST), row Z-scores, top 50 DEGs by adjusted p-value

## Repository Contents

- `AD-DEG-Analysis.R`: the analysis script

## Requirements

R with the following packages: DESeq2, dplyr, ggplot2, EnhancedVolcano, ComplexHeatmap, circlize. See `sessionInfo()` at the end of the script for the versions used.

## Input Data

Input data are **not** included in this repository. Place these files in the working directory:

- `metadata_filtered.txt`: tab-delimited sample metadata with columns `group` (sample name, matching the count-matrix column names), `diagnosis` (`control` / `Alzheimer`), `sex` (`female` / `male`), `age_at_death`, `rin_score` and `post_mortem_interval`. The first column is read as row names, so `group` must be a separate column.
- `gene_counts_brain.txt`: tab-delimited raw gene counts (genes in rows, samples in columns, Ensembl gene IDs as row names).

## Data source: *PRJNA377568*

## Usage

Open R in the folder containing the input files and run:

```r
source("AD-DEG-Analysis.R")
```

## Outputs

- `dds_dps.rds`: fitted DESeq2 object (reused on later runs; delete it if the inputs change)
- `result_dps_d_df_clean.txt`: full DESeq2 results table
- `MA_plot_DESeq2_ggplot.png`
- `Volcano_plot_diagnosis_with_covariates.png`
- `heatmap_top50_dps_d.png`

## License

MIT


## Publication

Beh-Pajooh A. Identification of candidate diagnostic biomarkers and gene networks for moderate stages of Alzheimer’s disease in fusiform gyrus exhibiting neurofibrillary tangles. Journal of Alzheimer’s Disease. 2026;112(2):667-689. doi:10.1177/13872877261452277
  
