# NEG-mode metabolomics visualisation workflow

This folder contains R scripts for visualisation of untargeted metabolomics data in **negative ionisation mode (NEG)**.
The workflow is organised into several analytical stages:

1. **Pre-annotation visualisation (feature level)**
2. **Post-annotation visualisation (metabolite and pathway level)**
   (Shiny application is documented separately and is not part of this folder)

This reflects the standard structure of untargeted metabolomics pipelines.

---

## 1. Why visualisation is performed *before* annotation

In untargeted LC–MS metabolomics, raw data first produce **features** (m/z–RT signals), not identified metabolites.
Before annotation can be meaningfully performed, it is essential to understand the dataset through visualisation.

### ✔ Data quality

* sample comparability
* effect of normalisation / transformation
* presence of outliers
* batch effects

### ✔ Biological structure

* whether stress vs control separation exists
* how many features change significantly
* general patterns in the dataset

### ✔ Statistical signal

* distribution of p-values
* fold-change patterns
* variance structure
* importance of features (VIP scores)

These questions do not require metabolite IDs — they rely only on:

* feature intensities
* statistical tests
* multivariate models (PCA/PLSDA)
* fold changes
* significance thresholds

Thus, **pre-annotation visualisation is standard practice** in untargeted metabolomics and necessary before investing effort into annotation.

---

## 2. Pre-annotation visualisation scripts (feature level)

These scripts operate on statistically processed feature tables (log2FC, p-values, VIP scores).
They provide exploratory analysis and QC-like inspection *before annotation*.

They are modular and can be run independently.

### `01_volcano_plot.R`

Volcano plot of NEG-mode features based on log2 fold change and adjusted p-values.
Shows overall strength of the biological signal.

### `02_heatmap_sig_features_neg.R`

Heatmap of significant features, displaying expression patterns across samples and potential clustering.

### `03_plsda_neg.R`

PLS-DA model of NEG-mode data, useful for visualising sample separation and identifying potential outliers.

### `03a_plsda_filterVIP_neg.R`

Filters features using VIP scores from PLS-DA.
Helps identify features most responsible for group differences.

### Summary of this stage

* **Purpose:** QC, exploratory visualisation, evaluation of biological effect
* **Inputs:** Feature table + statistical outputs
* **Outputs:** Volcano, heatmap, PLS-DA plots, VIP-filtered feature lists
* **Annotation:** *Not required*

---

## 3. Post-annotation visualisation (metabolite & pathway level)

After statistical exploration, selected or significant features undergo metabolite annotation (MS/MS matching, KEGG mapping, Mummichog enrichment).

The script below uses annotated metabolites and pathways for biological interpretation.

### `visual_after_annotation_neg.R`

This integrated script includes:

* annotated volcano plots
* pathway-level visualisation
* KEGG and Mummichog enrichment plots
* summarised post-annotation figures
* preparation of data for Shiny

### Purpose of this stage

* integrate statistics with metabolite IDs
* interpret biological mechanisms
* explore pathways and enriched processes
* generate final biologically meaningful figures

---

## 4. Relationship to the Shiny application

The post-annotation script prepares tables and outputs that are used by the interactive Shiny app:

**Metabolomics visualisation (NEG mode)**

The Shiny app allows:

* dynamic volcano plots
* interactive KEGG bubble plots
* chord diagrams for pathways
* adjustable cutoffs and display options
* fast exploration during analysis or presentation

---

## 5. Full workflow overview

Below is the complete analytical workflow relevant to this visualisation pipeline:

```
Processed data (after statistics)
        ↓
Pre-annotation visualisation (feature-level):
    • volcano plot
    • heatmap of significant features
    • PLS-DA and VIP filtering
        ↓
Metabolite annotation (MS/MS, KEGG, Mummichog)
        ↓
Post-annotation visualisation (metabolite & pathway level):
    • annotated volcano
    • KEGG bubble plot
    • pathway-level figures

```

This workflow reflects the real analytical structure of the project:

* **Pre-annotation visualisation** focuses on feature-level statistical exploration before metabolite IDs are available.
* **Annotation** links statistically significant features to metabolites, pathways, and biological processes.
* **Post-annotation visualisation** provides biological interpretation using annotated metabolites and enrichment results.


---

## 6. Suggested usage

### To explore data before annotation:

```r
source("01_volcano_plot.R")
source("02_heatmap_sig_features_neg.R")
source("03_plsda_neg.R")
source("03a_plsda_filterVIP_neg.R")
```

### After metabolite annotation:

```r
source("visual_after_annotation_neg.R")
```

### For interactive exploration:

Run the Shiny app located in `scripts/shiny/`.

---

## 7. Notes

* All scripts are designed specifically for **NEG mode** data.
* Input file paths may need adjustment to match your project structure.
* Pre-annotation visualisation reflects **exploration after statistical analysis**.
* Post-annotation visualisation reflects **biological interpretation**.
* The combined workflow follows best practices in untargeted metabolomics research.
