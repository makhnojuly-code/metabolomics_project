# ================================================================
# 08_stats_limma.R
# Build expression matrix, align with metadata, run limma
# ================================================================

# --- 0. Load required packages ---
library(dplyr)
library(limma)

message("Packages loaded.")


# --- 1. Load processed data ---
# We load:
#   - feature table after preprocessing and filtering
#   - sample metadata with group information

feat_mat <- read.csv(
  "data/processed/feature_table_neg_processed_for_stats_filtered.csv",
  check.names = FALSE
)

meta <- read.csv(
  "data/processed/sample_sheet_neg_filtered.csv",
  check.names = FALSE
)

# Quick look at the data structure
head(feat_mat)
head(meta)


# --- 2. Extract feature IDs (FT0001, FT0002, ...) ---
# The first column of feat_mat contains feature IDs.
feature_ids <- feat_mat[[1]]


# --- 3. Build expression matrix ---
# Remove the first column (feature IDs). Keep only intensities.
expr <- feat_mat[, -1]

# Use feature IDs as row names
rownames(expr) <- feature_ids

# Convert to numeric matrix (required by limma and many multivariate methods)
expr <- as.matrix(expr)


# --- 4. Align expression matrix with metadata ---
# We want the columns of expr to be in the same order as rows in meta.
# We assume that meta$file_name matches the column names of expr.

# Check if sample order already matches
all(colnames(expr) == meta$file_name)

# If sample order does not match, reorder columns of expr:
if (!all(colnames(expr) == meta$file_name)) {
  expr <- expr[, meta$file_name]
}

# Final check — should be TRUE now
stopifnot(all(colnames(expr) == meta$file_name))

message("Expression matrix and metadata are aligned.")


# --- 5. Quick check of metadata ---
# meta is expected to contain a column 'group' with values:
#   "Null" and "Wounding"
colnames(meta)
table(meta$group)   # How many samples per group?


# --- 6. Create group factor for limma ---
# Convert the 'group' column into a factor that limma can use.
group <- factor(meta$group)

# Check the distribution of samples across groups
table(group)
levels(group)   # expected: "Null" "Wounding"


# --- 7. Build design matrix ---
# The design matrix encodes which samples belong to which group.
# Formula '~ 0 + group' means:
#   - no intercept term (0)
#   - one column per group (Null, Wounding)
design <- model.matrix(~ 0 + group)

# Name the columns after the group levels
colnames(design) <- levels(group)

# Inspect the first rows of the design matrix
head(design)


# --- 8. Fit the linear model with limma ---
# limma fits one linear model per feature (row of expr):
#   intensity_feature ~ group
fit <- limma::lmFit(expr, design)


# --- 9. Define contrast of interest ---
# We want to compare: Wounding vs Null
# Interpretation:
#   logFC > 0  → feature higher in Wounding
#   logFC < 0  → feature lower in Wounding
contrast_matrix <- limma::makeContrasts(
  Wounding_vs_Null = Wounding - Null,
  levels = design
)

# Apply contrast and empirical Bayes moderation
fit2 <- limma::contrasts.fit(fit, contrast_matrix)
fit2 <- limma::eBayes(fit2)


# --- 10. Extract full results table ---
# topTable() returns, for each feature:
#   - logFC        : log2 fold change (Wounding vs Null)
#   - AveExpr      : average intensity across all samples
#   - t            : moderated t-statistic
#   - P.Value      : raw p-value
#   - adj.P.Val    : FDR-adjusted p-value (Benjamini–Hochberg)
#   - B            : log-odds of differential expression
res <- limma::topTable(
  fit2,
  coef = "Wounding_vs_Null",  # name of the contrast
  number = Inf,               # return all features
  adjust.method = "fdr"
)

# Quick look at the top rows
head(res)


# --- 11. Save results to CSV ---
if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}

write.csv(
  res,
  file = "results/limma_results_neg.csv",
  row.names = TRUE
)

message("Limma differential analysis finished.")
message("Results saved to 'results/limma_results_neg.csv'.")



