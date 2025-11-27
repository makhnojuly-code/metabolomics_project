# NEG Annotation Pipeline

This folder contains the **NEG-mode annotation pipeline** for LC–MS data.  
The pipeline performs:

- feature-level statistics (VIP + LIMMA)
- preparation of Mummichog input
- metabolite annotation (Mummichog, KEGG)
- selection of stress-core metabolites
- validation of the final tables

The workflow consists of **six main scripts** (+ one optional automatic Mummichog step).

---

## Overview of scripts

### 01_build_annotation_candidates_neg.R  
**Feature statistics & preparation of Mummichog input**

This script combines:

- PLS-DA VIP scores  
- LIMMA differential statistics  
- LC–MS feature metadata (m/z, RT)  

It produces:

- `annotation_candidates_neg.csv`  
- `metaboanalyst_input_neg.csv`  
- `metabo_mummichog_NEG.txt`  ← input file for Mummichog  

No annotation is performed here — this step only prepares and structures the input data.

---

### 02_export_mummichog_input_neg.R  
**(historical) Manual Mummichog preparation**

This script converts the simplified feature table into the text format required by the Mummichog module:

- columns: `m.z`, `p.value`, `t.score`, `r.t`  

The generated `.txt` file can be uploaded manually to the MetaboAnalyst web interface.

After running Mummichog on the web, MetaboAnalyst returns two key files:

- `mummichog_matched_compound_all.csv` — matched compounds  
- `mummichog_pathway_enrichment_integ.csv` — pathway enrichment  

These files can be placed into `results/annotation/`.  
The downstream annotation starts from this point.

> **Note:** this manual step is now optional. It is replaced by Script **02A**, which runs the same analysis automatically in R.

---

### 02A_neg_annotation_auto.R  
**Automatic Mummichog (MetaboAnalystR)**  
*(replaces the manual web step)*

This script:

1. Loads `results/annotation/metabo_mummichog_NEG.txt`  
2. Creates a MetaboAnalystR object for `"mass_all"` / `"mummichog"`  
3. Applies the same processing parameters as the online tool  
   - ppm tolerance  
   - ion mode (NEG)  
   - enrichment mode (`integ`, v2)  
   - p-value cutoff  
4. Runs Mummichog locally  
5. Detects the output CSV files in the working directory  

It automatically generates the same two key outputs as the web version:

- `mummicho_matched_compound_all.csv` — matched compounds  
- `mummicho_pathway_enrichment.csv` — pathway enrichment  

These files are saved into `results/annotation/` and can be used directly in Scripts 03–06.  
This makes the NEG annotation step **fully reproducible and automated**.

---

### 03_integrate_mummichog_neg.R  
**Integrate Mummichog annotation with feature statistics**

This script links LC–MS features with KEGG compounds using rounded m/z and retention time.

It:

- reads `annotation_candidates_neg.csv`  
- reads `mummicho_matched_compound_all.csv` *(or the manual file `mummichog_matched_compound_all.csv`)*  
- rounds m/z and RT in both tables  
- performs a left join on `(m/z, RT)`  

It produces:

- `features_with_mummichog_neg.csv` — full table with statistics + KEGG matches  
- `mummichog_hits_neg.csv` — only features that received KEGG annotation  

This is the first step where **raw LC–MS features become biological metabolites**.

---

### 04_add_kegg_info_neg.R  
**Add KEGG chemical information**

Using `KEGGREST`, this script retrieves for each KEGG ID:

- compound name  
- molecular formula  
- exact mass  

The KEGG information is merged into the annotated feature table.  
The resulting table contains:

- LIMMA statistics  
- VIP values  
- m/z and RT  
- KEGG IDs  
- KEGG chemical properties  

This creates a complete dataset for biological interpretation and pathway analysis.

---

### 05_stress_core_filter_neg.R  
**Identify the stress-core metabolome**

Core features are selected based on thresholds for:

- adjusted p-value  
- \|logFC\|  
- VIP  

From these, the script builds:

- `core_features_neg.csv` — core LC–MS features  
- `core_metabolites_neg.csv` — aggregated metabolite-level table (by KEGG ID)  

For each metabolite it computes:

- number of supporting features  
- minimum FDR  
- mean and maximum logFC  
- maximum VIP  

Additionally, each metabolite receives:

- a direction label: **“Up in stress”** / **“Down in stress”**  
- a short metabolite name  
- a simplified biological category in `stress_core_group`  
  (e.g. JA/oxylipins, flavonoids, aromatic amino acids, sugar phosphates, etc.)

---

### 06_validate_core_neg.R  
**Validation of core tables**

This script performs a series of checks on:

- `core_features_neg.csv`  
- `core_metabolites_neg.csv`  

It verifies:

- all required columns are present  
- no duplicated KEGG IDs  
- NA values are reasonable  
- direction labels are consistent with `mean_logFC`  
- all metabolites in the core table come from core features  
- distribution of `stress_core_group` categories  

Output:

- `results/annotation/validation_report_neg.txt`  

This confirms that the NEG annotation tables are **ready for visualization** and for downstream biological interpretation.

---

## Notes

- The pipeline is designed for **NEG ion-mode LC–MS data**.  
- Automatic Mummichog (Script 02A) uses the **Arabidopsis KEGG library** (`ath_kegg`).  
- Manual and automatic Mummichog outputs have been cross-checked and give **identical results**.
