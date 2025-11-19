# ================================================================
# 03_plsda_neg.R
# Supervised multivariate analysis (PLS-DA) for NEG mode
#
# Goal:
#   Build and validate a supervised model that separates Null vs Wounding
#   based on the full metabolite profile, and extract VIP scores and
#   basic performance metrics.
#
# Inputs:
#   data/processed/feature_table_neg_processed_for_stats_filtered.csv
#   data/processed/sample_sheet_neg_filtered.csv
#
# Outputs (tables):
#   results/plsda_scores_neg.csv         - sample scores on PLS-DA components
#   results/plsda_vip_neg.csv            - VIP scores per feature
#   results/plsda_confusion_neg.csv      - confusion matrix
#   results/plsda_cv_error_neg.csv       - CV error per component
#   results/plsda_permutation_neg.csv    - permutation CV errors
#   results/plsda_perf_summary_neg.csv   - summary: accuracy, sens, spec, p_perm
#
# Outputs (figures):
#   results/viz/plsda_scores_neg.png     - PLS-DA score plot (comp 1 vs 2)
#   results/viz/plsda_cv_error_neg.png   - CV error vs number of components
#   results/viz/plsda_permutation_neg.png- permutation-test histogram
#
# High-level steps:
#   1) Load processed data and build feature x sample matrix.
#   2) Align expression matrix with metadata (sample order).
#   3) Fit a PLS-DA model (2 components).
#   4) Extract sample scores + VIP scores and save as CSV.
#   5) Compute confusion matrix, accuracy, sensitivity, specificity.
#   6) Perform repeated M-fold cross-validation (perf).
#   7) Run a simple permutation test and estimate empirical p-value.
#   8) Generate three plots (scores, CV, permutation).
# ================================================================

# --- 0. Load required packages ----------------------------------

# mixOmics must be installed beforehand, for example:
# if (!requireNamespace("mixOmics", quietly = TRUE)) {
#   BiocManager::install("mixOmics")
# }

library(mixOmics)
library(dplyr)

message("Packages for PLS-DA loaded.")


# --- 1. Define input/output paths -------------------------------

expr_path <- "data/processed/feature_table_neg_processed_for_stats_filtered.csv"
meta_path <- "data/processed/sample_sheet_neg_filtered.csv"

out_dir_results <- "results"
out_dir_viz     <- "results/viz"

dir.create(out_dir_results, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_viz, recursive = TRUE, showWarnings = FALSE)


# --- 2. Load expression matrix and metadata ---------------------

feat_mat <- read.csv(expr_path, check.names = FALSE)
meta     <- read.csv(meta_path, check.names = FALSE)

message("Expression table rows (features): ", nrow(feat_mat),
        " | columns (1 ID + samples): ", ncol(feat_mat))
message("Metadata rows (samples): ", nrow(meta))


# --- 3. Build expression matrix (features x samples) ------------

# First column = feature IDs (e.g. FT0001, FT0002, ...)
feature_ids <- feat_mat[[1]]

# Drop ID column, keep only numeric intensities
expr <- as.matrix(feat_mat[, -1])
rownames(expr) <- feature_ids   # rows = features, columns = samples


# --- 4. Align expression matrix with metadata -------------------

# We expect a column "file_name" in metadata that matches colnames(expr)
stopifnot("file_name" %in% names(meta))

if (!all(colnames(expr) == meta$file_name)) {
  expr <- expr[, meta$file_name]
}

# Safety check: columns must now be in the same order as rows in meta
stopifnot(all(colnames(expr) == meta$file_name))

message("Expression matrix and metadata are aligned.")


# --- 5. Define response (class labels) --------------------------

# PLS-DA needs a factor with group labels, e.g. "Null" / "Wounding"
group <- factor(meta$group)

message("Class distribution:")
print(table(group))


# --- 6. Fit PLS-DA model (2 components) -------------------------

ncomp <- 2          # number of PLS-DA components
set.seed(123)       # for reproducibility

plsda_model <- plsda(
  X     = t(expr),  # mixOmics expects samples in rows → transpose
  Y     = group,
  ncomp = ncomp
)

message("PLS-DA model fitted with ", ncomp, " components.")


# --- 7. Sample scores (for plotting / downstream use) -----------

scores <- plsda_model$variates$X    # samples x components

scores_df <- data.frame(
  sample = rownames(scores),
  group  = group,
  scores
)

write.csv(
  scores_df,
  file = file.path(out_dir_results, "plsda_scores_neg.csv"),
  row.names = FALSE
)

message("PLS-DA scores saved to: results/plsda_scores_neg.csv")


# --- 8. Score plot (component 1 vs 2) ---------------------------

score_plot_file <- file.path(out_dir_viz, "plsda_scores_neg.png")

png(score_plot_file, width = 900, height = 700, res = 120)
plotIndiv(
  plsda_model,
  comp       = c(1, 2),
  group      = group,
  ind.names  = FALSE,
  legend     = TRUE,
  title      = "PLS-DA (NEG mode): Wounding vs Null"
)
dev.off()

message("PLS-DA score plot saved to: ", score_plot_file)


# --- 9. VIP scores (feature importance) -------------------------

# VIP scores: importance of each feature for the PLS-DA model
vip_scores <- vip(plsda_model)   # matrix: features x components

# Mean VIP across components (common summary measure)
vip_mean <- rowMeans(vip_scores, na.rm = TRUE)

vip_df <- data.frame(
  feature_id = rownames(vip_scores),
  vip_comp1  = vip_scores[, 1],
  vip_comp2  = if (ncomp >= 2) vip_scores[, 2] else NA,
  vip_mean   = vip_mean
) %>%
  arrange(desc(vip_mean))

write.csv(
  vip_df,
  file = file.path(out_dir_results, "plsda_vip_neg.csv"),
  row.names = FALSE
)

n_vip_gt1 <- sum(vip_df$vip_mean > 1, na.rm = TRUE)
message("VIP table saved to: results/plsda_vip_neg.csv")
message("Number of features with VIP_mean > 1: ", n_vip_gt1)


# --- 10. Confusion matrix, accuracy, sensitivity, specificity ---

# Predict class labels for the same samples used to fit the model
pred <- predict(plsda_model, newdata = t(expr), dist = "max.dist")

# mixOmics stores predictions per distance and component;
# we take the final component (ncomp) and "max.dist" distance
pred_class <- pred$class$max.dist[, ncomp]

conf_mat <- table(
  Truth     = group,
  Predicted = pred_class
)

message("Confusion matrix (training data):")
print(conf_mat)

accuracy <- sum(diag(conf_mat)) / sum(conf_mat)

# Define "positive" = Wounding, "negative" = Null
tp <- conf_mat["Wounding", "Wounding"]
tn <- conf_mat["Null",     "Null"]
fp <- conf_mat["Null",     "Wounding"]
fn <- conf_mat["Wounding", "Null"]

sensitivity <- tp / (tp + fn)  # recall for Wounding
specificity <- tn / (tn + fp)  # correctly predicted Null

perf_summary <- data.frame(
  metric = c("accuracy", "sensitivity_Wounding", "specificity_Null"),
  value  = c(accuracy, sensitivity, specificity)
)

write.csv(
  as.data.frame.matrix(conf_mat),
  file = file.path(out_dir_results, "plsda_confusion_neg.csv"),
  row.names = TRUE
)

write.csv(
  perf_summary,
  file = file.path(out_dir_results, "plsda_perf_summary_neg.csv"),
  row.names = FALSE
)

message("Accuracy:      ", round(accuracy    * 100, 2), "%")
message("Sensitivity (Wounding): ", round(sensitivity * 100, 2), "%")
message("Specificity (Null):     ", round(specificity * 100, 2), "%")
message("Confusion matrix and performance summary saved.")


# --- 11. Cross-validation (M-fold, repeated) -------------------

# Here we estimate how stable the model is by repeated M-fold CV.
# We compute overall CV error for 1..ncomp components.

cv_res <- perf(
  plsda_model,
  validation  = "Mfold",
  folds       = 5,
  nrepeat     = 20,
  dist        = "max.dist",
  progressBar = FALSE
)

# overall error rate per component (for "max.dist" distance)
# rows = components, column = "max.dist"
cv_overall <- cv_res$error.rate$overall[, "max.dist"]
cv_overall <- as.numeric(cv_overall)
comp_index <- seq_along(cv_overall)

# Save CV errors to CSV
cv_df <- data.frame(
  component      = comp_index,
  overall_cv_err = cv_overall
)

write.csv(
  cv_df,
  file = file.path(out_dir_results, "plsda_cv_error_neg.csv"),
  row.names = FALSE
)

# Plot CV error vs number of components
cv_plot_file <- file.path(out_dir_viz, "plsda_cv_error_neg.png")

png(cv_plot_file, width = 900, height = 700, res = 120)
plot(
  x    = comp_index,
  y    = cv_overall,
  type = "b",
  xlab = "Number of PLS-DA components",
  ylab = "Overall CV error (max.dist)",
  main = "Cross-validation error by number of components"
)
abline(v = ncomp, col = "red", lty = 2)  # vertical line at chosen ncomp
dev.off()

message("Cross-validation error plot saved to: ", cv_plot_file)

# Overall CV error for the chosen model (final component)
overall_err <- cv_overall[ncomp]


# --- 12. Permutation test (simple implementation) --------------

# Idea:
#   1) Randomly shuffle group labels (break the real structure).
#   2) Refit PLS-DA + run CV, collect CV error.
#   3) Repeat many times (n_permutations).
#   4) Compare real error to distribution of permuted errors.
#      → empirical p-value tells how likely such separation is by chance.

n_permutations <- 100 # increase for more stable p_perm, decrease if runtime is long
perm_errors    <- numeric(n_permutations)

set.seed(456)  # separate seed for permutation test

for (i in seq_len(n_permutations)) {
  # 1) shuffle labels
  group_perm <- sample(group)
  
  # 2) fit PLS-DA on permuted labels
  plsda_perm <- plsda(
    X     = t(expr),
    Y     = group_perm,
    ncomp = ncomp
  )
  
  # 3) CV performance for permuted model
  cv_perm <- perf(
    plsda_perm,
    validation  = "Mfold",
    folds       = 5,
    nrepeat     = 20,
    dist        = "max.dist",
    progressBar = FALSE
  )
  
  # Error for the final component
  perm_errors[i] <- cv_perm$error.rate$overall[ncomp, "max.dist"]
}

# Empirical p-value: proportion of permuted errors
# that are less than or equal to the real error
p_perm <- (sum(perm_errors <= overall_err) + 1) / (n_permutations + 1)

message("Permutation test completed with ", n_permutations, " permutations.")
message("Real overall CV error:        ", round(overall_err, 4))
message("Mean permuted overall CV error: ", round(mean(perm_errors), 4))
message("Empirical permutation p-value: ", signif(p_perm, 3))

# Save permutation errors to CSV
perm_df <- data.frame(
  perm_id       = seq_len(n_permutations),
  overall_error = perm_errors
)

write.csv(
  perm_df,
  file = file.path(out_dir_results, "plsda_permutation_neg.csv"),
  row.names = FALSE
)

# Add permutation summary into performance table
perm_summary <- data.frame(
  metric = c("overall_cv_error_real", "mean_cv_error_permuted", "p_perm"),
  value  = c(overall_err, mean(perm_errors), p_perm)
)

perf_all <- rbind(
  perf_summary,
  perm_summary
)

write.csv(
  perf_all,
  file = file.path(out_dir_results, "plsda_perf_summary_neg.csv"),
  row.names = FALSE
)

# --- 13. Permutation plot with REAL error line -----------------

perm_plot_file <- file.path(out_dir_viz, "plsda_permutation_neg.png")

# Define pretty x-limits so the histogram is close to the red line
x_min <- min(c(perm_errors, overall_err)) - 0.02
x_max <- max(c(perm_errors, overall_err)) + 0.02

png(perm_plot_file, width = 900, height = 700, res = 120)

hist(
  perm_errors,
  breaks = 20,
  col    = rgb(110, 160, 220, 120, maxColorValue=255),   # soft blue
  border = "white",
  xlab   = "Overall CV error (permuted labels)",
  main   = paste0("Permutation Test (", n_permutations, " permutations)"),
  xlim   = c(x_min, x_max)
)

abline(v = overall_err, col = "red", lwd = 3)

legend(
  "topleft",
  inset = c(0.1,0),
  legend = c(
    paste0("Real CV error = ", round(overall_err, 4)),
    paste0("Permutation p-value = ", signif(p_perm, 3))
  ),
  bty = "n"
)

dev.off()


message("Permutation histogram saved to: ", perm_plot_file)
message("PLS-DA analysis (model + CV + permutation) completed.")