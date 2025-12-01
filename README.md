# Metabolomics Project (LC–MS, NEG Mode)

Author: **Yuliia Makhno**
Host institution: **IGZ Großbeeren (Science Support Platform)**
Focus: Untargeted metabolomics (NEG mode): preprocessing, QC, normalisation, statistics, annotation, visualisation

This repository contains a complete and reproducible workflow for untargeted LC–MS metabolomics in negative ionisation mode. The project covers all major stages:

* preprocessing of raw LC–MS data
* QC and evaluation of normalisation strategies
* statistical analysis (multivariate and limma)
* annotation (Mummichog and KEGG)
* identification of stress-core metabolites
* visualisation (feature-level and pathway-level)
* optional Shiny application for interactive exploration

---

# Data Source

This project uses publicly available LC–MS metabolomics data from MetaboLights:

**MTBLS4823 – "The suberin transporter StABCG1 is required for barrier formation in potato leaves."**
Authors: Steffen Neumann, Elvio H. Benatto Perino, Ulrike Smolka, Karin Gorzelka, Ramona Gritzner, Sylvestre Marillonnet, Khabat Vahabi, Sabine Rosahl
Link: [https://www.ebi.ac.uk/metabolights/MTBLS4823](https://www.ebi.ac.uk/metabolights/MTBLS4823)

For this internship project:

* only NEG-mode LC–MS data were used
* comparison groups: Wounding vs Null (control)
* organism: *Solanum tuberosum* (potato)
* transcriptomic data were not analysed in this project

---

# Project Structure

```
metabolomics_project/
├── data/
│   ├── raw/              # raw .mzML files
│   ├── intermediate/     # RDS objects from xcms
│   └── processed/        # feature tables for statistics
│
├── scripts/
│   ├── preprocessing/    # xcms/MSnbase pipeline
│   ├── stats/            # cleaning → QC → normalisation → limma
│   ├── annotation/       # Mummichog, KEGG, core metabolites
│   └── viz/              # feature-level and pathway-level visualisation
│
├── results/
│   ├── qc/               # QA/QC plots
│   ├── annotation/       # Mummichog and KEGG outputs
│   └── plots/            # visualisation outputs
│
├── docs/                 # optional notes
├── quarto/               # reports
└── README.md
```

Each major folder contains its own README with detailed explanations.

---

# How to Use This Repository

## 1. Clone

```
git clone https://github.com/makhnojuly-code/metabolomics_project.git
cd metabolomics_project
```

## 2. Place raw LC–MS data

Place `.mzML` files into:

```
data/raw/neg_mode/
```

## 3. Run preprocessing

```
source("scripts/preprocessing/00_setup.R")
source("scripts/preprocessing/01_create_sample_sheet_neg.R")
source("scripts/preprocessing/01_preprocessing.R")
source("scripts/preprocessing/01a_qc_plots_fast.R")
source("scripts/preprocessing/01b_qc_plots_full.R")
```

Outputs:

* processed feature tables → `data/processed/`
* QC plots → `results/qc/`

---

## 4. Run statistical analysis

```
source("scripts/stats/01_cleaning_neg.R")
source("scripts/stats/02_before_normalization_neg.R")
source("scripts/stats/03_screen_normalization_methods_neg.R")
source("scripts/stats/04_normalization_neg.R")
source("scripts/stats/05_pca_after_normalization_neg.R")
source("scripts/stats/06_outlier_check_neg.R")
source("scripts/stats/07_filter_after_outlier_neg.R")
source("scripts/stats/08_stats_limma.R")
```

Main output: `results/stats/limma_results_neg.csv`

---

## 5. Run annotation pipeline

```
source("scripts/annotation/01_build_annotation_candidates_neg.R")
source("scripts/annotation/02A_neg_annotation_auto.R")
source("scripts/annotation/03_integrate_mummichog_neg.R")
source("scripts/annotation/04_add_kegg_info_neg.R")
source("scripts/annotation/05_stress_core_filter_neg.R")
source("scripts/annotation/06_validate_core_neg.R")
```

Outputs stored in `results/annotation/`.

---

# Full Workflow Overview

## A. LC–MS Preprocessing

Detailed in `scripts/preprocessing/README.md`.

* setup (packages, paths, colours)
* sample sheet generation
* mzML import (onDisk)
* CentWave peak picking
* Obiwarp retention-time alignment
* feature grouping
* missing value filling
* QC plots (TIC/BPC, RT shift)

Outputs: `data/processed/`, `results/qc/`.

## B. Statistical Pipeline

Detailed in `scripts/stats/README.md`.

1. Cleaning (remove technical samples, filter NAs, impute)
2. QC before normalisation (density, boxplots, PCA)
3. Screening of normalisation methods
4. Final normalisation: PQN → log2 → Pareto
5. PCA QC after normalisation
6. Outlier detection
7. Remove outliers and rebuild final dataset
8. Differential analysis (limma)

Output: `limma_results_neg.csv`.

## C. Annotation Pipeline

Detailed in `scripts/annotation/README.md`.

* build annotation candidates (VIP + limma)
* prepare input for Mummichog
* automatic Mummichog via MetaboAnalystR
* KEGG compound information
* stress-core metabolite identification
* validation checks

Outputs: annotated tables in `results/annotation/`.

## D. Visualisation Workflow

Detailed in `scripts/viz/README.md`.

### Pre-annotation visualisation

* volcano plot
* heatmap of significant features
* PLS-DA
* VIP filtering

### Post-annotation visualisation

* annotated volcano
* KEGG bubble plot
* pathway summaries
* chord diagrams (optional)

### Shiny app

Interactive exploration of:

* core metabolites
* KEGG enrichment
* pathway connections

App located in `scripts/shiny/`.

---

# Key R Packages Used

Preprocessing: `xcms`, `MSnbase`, `mzR`, `BiocParallel`
Statistics: `dplyr`, `limma`, `ggplot2`, `FactoMineR`, `mixOmics`
Annotation: `MetaboAnalystR`, `KEGGREST`, `jsonlite`
Visualisation: `ggplot2`, `pheatmap`, `circlize`, `ComplexHeatmap`, `tidyverse`, `shiny`, `highcharter`

---

# License

This project is released under the MIT License (see LICENSE file).
