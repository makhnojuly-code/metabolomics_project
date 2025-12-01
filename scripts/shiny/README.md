# Metabolomics visualisation (NEG mode) – Shiny app

This Shiny app shows an interactive visualisation for **untargeted metabolomics in NEG mode**.
It is designed as a **preview of possible plots** for stress vs. control comparisons.

The app includes:

* A **volcano plot** for core metabolites
* A **KEGG enrichment bubble plot**
* A **chord diagram** for the top 10 KEGG pathways (P1–P10) and regulation (Up / Down in stress)

> ⚠️ This is a **prototype / preview**, not a final polished tool.
> It shows what is possible and can be improved in the future.

---

## 1. Input data

The app uses **four CSV files** as input. The file paths are currently hard-coded in the script (not selected through the UI), 
but you can change them to match your own project structure.

All input files should:

* be in **CSV format**,
* contain the **required columns** described in sections 1.1–1.4,
* use **consistent metabolite identifiers** across tables (for example, `Empirical.Compound` should match the IDs listed in `cpd.hits`),
* correspond to **NEG mode** metabolomics data (stress vs. control comparison).

In short, these four tables provide:

* per-metabolite statistics (p-values and log2 fold changes),
* pathway enrichment results (KEGG / Mummichog / GSEA output),
* mapping between pathways and detected metabolites,
* the direction of regulation (Up / Down in stress).

The following subsections describe the expected format of each input file.


### 1.1 Core metabolites (for the volcano plot)

**Path in the script:**

```r
core_path <- "results/stress_core_neg/core_metabolites_neg.csv"
```
This file contains the list of core metabolites detected in NEG mode, together with their statistical values. 
The app uses this table to build the volcano plot.

**Required columns:**

* **min_adjP** — minimal adjusted p-value for the metabolite
* **mean_logFC** — mean log2 fold change (stress vs. control)

Additional columns (metabolite names, IDs, annotations) may be present, but they are not required for plotting.

These values are used to calculate:

* **neg_log10_p** = –log10(min_adjP)
* **regulation class** (Up in stress / Down in stress / Not significant), based on user-defined thresholds for p-value and |log2FC|

The volcano plot visualises all metabolites according to their fold change and statistical significance.

### 1.2 KEGG pathway enrichment (for the bubble plot)

**Path in the script:**

```r
enrich_path <- "results/annotation/mummicho_pathway_enrichment.csv"
```
This file contains KEGG/Mummichog/GSEA enrichment results. The app uses it to display enriched pathways in the bubble plot.

**Required columns:**

* **Total_Size** – total number of metabolites in the pathway
* **Hits** – detected metabolites mapped to that pathway
* **Sig_Hits** – significant metabolites
* at least one enrichment p-value column (e.g. **Combined_Pvals**)

The app automatically computes standard metrics (such as rich factor and –log10 p-values) and uses 
them for visualisation. Only the required columns above must be provided; the rest is handled internally.
