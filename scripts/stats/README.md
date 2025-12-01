# Statistical Analysis Pipeline (NEG mode)

*Cleaning → QC → Normalization screening → Final normalization → PCA QC → Outliers → Final dataset → limma*

This section describes the complete statistical workflow used for the metabolomics project.
The goal is to ensure that differential analysis (limma) is performed on **clean, normalized, and biologically stable** data.

Below is a clear, step-by-step explanation of the entire pipeline.

---

## Full Statistical Pipeline — Simple Explanation

### **1. Cleaning the feature table**

*(script `01_cleaning_neg.R`)*

**Goal:** produce a clean and reliable data matrix.

**What happens:**

* remove technical samples (QC, Blank, MM8, ACN)
* compute missing values per feature and per sample
* save QC histograms
* remove features with **>30% NAs**
* impute remaining NAs using **half-minimum** (biologically correct for LC–MS)
* save cleaned matrix + summary

**Why:**
Unclean data makes later analysis unstable or misleading.

---

### **2. QC before normalization**

*(script `02_before_normalization_neg.R`)*

**Goal:** inspect raw data quality and check if normalization is needed.

**We evaluate:**

* density plots per feature
* overall log-density
* sample boxplots (log scale)
* PCA before normalization

**We answer:**

* Are sample intensities comparable?
* Are distributions similar?
* Are there outliers already?
* Do we see group separation before normalization?

---

### **3. Screening normalization methods**

*(script `03_screen_normalization_methods_neg.R`)*

**Goal:** objectively choose the best normalization method.

**Procedure:**

* test multiple methods: raw, log2, Pareto, autoscale, PQN, arcsinh, etc.
* for each method:

  * apply normalization
  * run PCA
  * compute a "separation score" (Null vs Wounding)
  * save PCA plots

**Outcome:** allows an unbiased choice.

**Best method in this dataset:**
 **PQN → log2 → Pareto scaling**

---

### **4. Final normalization**

*(script `04_normalization_neg.R`)*

**Pipeline:**

1. **PQN normalization** – corrects for sample-wide dilution / injection variability
2. **log2 transform** – stabilizes variance, improves normality
3. **Pareto scaling** – balances impact of strong and weak features

**Saved:**

* PQN factors
* reference medians
* epsilon used for log2
* QC summaries
* final normalized table

---

### **5. PCA QC after normalization**

*(script `05_pca_after_normalization_neg.R`)*

**Goal:** verify that normalization produced stable data.

**We check:**

* PCA shape
* sample grouping
* sample sheet alignment
* absence of technical injections
* absence of unexpected patterns

---

### **6. Outlier detection**

*(script `06_outlier_check_neg.R`)*

We identify problematic samples via:

1. **PC1/PC2 thresholds**
2. **Mahalanobis distance** (distance in multivariate space)

Samples flagged by multiple criteria are considered outliers.

Output: a table of flagged samples.

---

### **7. Filtering outliers + final QC**

*(script `07_filter_after_outlier_neg.R`)*

**Goal:** prepare the final dataset for statistics.

**Steps:**

* remove outliers and technical samples
* rebuild clean sample sheet
* PCA after filtering
* QC histograms
* save final filtered matrix and metadata

This dataset is used for differential analysis.

---

### **8. Differential expression analysis (limma)**

*(script `08_stats_limma.R`)*

**Steps:**

1. build expression matrix
2. align sample order with metadata
3. create design matrix (Wounding vs Null)
4. fit limma linear models
5. apply empirical Bayes moderation
6. export results: logFC, p-values, FDR

**Output:** `results/limma_results_neg.csv`

---
