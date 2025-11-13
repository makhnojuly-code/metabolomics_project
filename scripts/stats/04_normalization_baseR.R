# -------------------------------------------------------------
# 
# Purpose:
#   Normalize cleaned NEG LC–MS feature table for downstream
#   statistics (PCA, univariate tests, etc.).
#
#   Chosen pipeline (based on 02_screen_normalization_methods_neg.R):
#     1) PQN   – correct for total concentration / injection differences
#     2) log2  – compress large values and stabilize variance
#     3) Pareto – balance feature influence for multivariate models
#
# Input:
#   data/processed/feature_table_neg_clean.csv   (features x samples)
#     - technical injections already removed
#     - missing values already imputed (half-min per feature)
#
# Output:
#   data/processed/feature_table_neg_processed_for_stats.csv
#
# QC / logs:
#   results/qc/normalization_neg/
#     - pqn_factors.csv
#     - pqn_reference_row_medians.csv
#     - eps_used.csv
#     - normalization_summary.csv
#     - normalization_params.csv
#   results/logs/sessionInfo_normalization_<date>.txt
# -------------------------------------------------------------

source("scripts/preprocessing/00_setup.R")

# -------------------- Settings ---------------------------------

IN_CSV   <- here::here("data","processed","feature_table_neg_clean.csv")
OUT_CSV  <- here::here("data","processed","feature_table_neg_processed_for_stats.csv")

QC_DIR   <- here::here("results","qc","normalization_neg")
LOGS_DIR <- here::here("results","logs")
dir.create(QC_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(LOGS_DIR, recursive = TRUE, showWarnings = FALSE)

# How to handle tiny/zero values before log2
EPS_MODE   <- "fixed"      # "fixed" or "adaptive"
EPS_FIXED  <- 1e-12        # used when EPS_MODE = "fixed"

# Pareto type
PARETO_MODE <- "classic"   # "classic" (mean & SD) or "robust" (median & MAD)

# -------------------- Load data --------------------------------

stopifnot(file.exists(IN_CSV))
X <- read.csv(IN_CSV, row.names = 1, check.names = FALSE)

# Convert safely to numeric (CSV may read as character)
X[] <- lapply(X, function(x) suppressWarnings(as.numeric(x)))
X <- as.matrix(X)

cat("Input features:", nrow(X), "\n")
cat("Input samples :", ncol(X), "\n")

# Optional check: technical samples should already be removed
if (any(grepl("(?i)QC|MM8|ACN|Blank", colnames(X)))) {
  warning("Technical samples (QC/MM8/ACN/Blank) detected. ",
          "They should be removed in the cleaning step.")
}

# -------------------- Step 1: PQN -------------------------------
# Correct for different total concentrations between samples

# 1) Row-level reference (median per feature)
ref <- apply(X, 1, median, na.rm = TRUE)
ref[!is.finite(ref) | ref <= 0] <- 1e-12  # avoid division by zero / Inf

# 2) Divide each feature by its reference
Xn <- sweep(X, 1, ref, "/")

# 3) PQN factor per sample = median of quotients (finite only)
pqn_factor <- apply(
  Xn, 2,
  function(v) median(v[is.finite(v)], na.rm = TRUE)
)
pqn_factor[!is.finite(pqn_factor) | pqn_factor == 0] <- 1

# 4) Divide each sample by its PQN factor
Xpqn <- sweep(Xn, 2, pqn_factor, "/")

# Sanity check
if (!all(is.finite(pqn_factor))) {
  stop("Non-finite values in PQN factors. Check input/reference medians.")
}

# Save PQN info (QC)
write.csv(
  data.frame(sample = colnames(Xpqn), pqn_factor = as.numeric(pqn_factor)),
  file.path(QC_DIR, "pqn_factors.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(feature = rownames(Xpqn), ref_median = as.numeric(ref)),
  file.path(QC_DIR, "pqn_reference_row_medians.csv"),
  row.names = FALSE
)

# -------------------- Step 2: log2 transform --------------------
# Reduce dominance of very large peaks and stabilize variance

if (EPS_MODE == "adaptive") {
  min_pos <- suppressWarnings(
    min(Xpqn[Xpqn > 0 & is.finite(Xpqn)], na.rm = TRUE)
  )
  if (!is.finite(min_pos)) min_pos <- EPS_FIXED
  eps <- min(EPS_FIXED, min_pos / 2)
} else {
  eps <- EPS_FIXED
}

# Replace zeros / non-finite with tiny positive number
Xpqn[!is.finite(Xpqn) | Xpqn <= 0] <- eps
logX <- log2(Xpqn)

# Sanity check
if (!all(is.finite(logX))) {
  stop("Non-finite values after log2. Check eps handling.")
}

# Save eps (QC)
write.csv(
  data.frame(metric = "eps_used", value = eps),
  file.path(QC_DIR, "eps_used.csv"),
  row.names = FALSE
)

# -------------------- Step 3: Pareto scaling --------------------
# Make features contribute more equally in PCA / PLS-DA

if (PARETO_MODE == "robust") {
  # Robust: center = median, scale = MAD
  row_center <- apply(logX, 1, median, na.rm = TRUE)
  row_scale  <- apply(logX, 1, mad,    na.rm = TRUE, constant = 1)
} else {
  # Classic: center = mean, scale = SD
  row_center <- rowMeans(logX, na.rm = TRUE)
  row_scale  <- apply(logX, 1, sd, na.rm = TRUE)
}
row_scale[!is.finite(row_scale) | row_scale == 0] <- 1

# Apply Pareto: (value - center) / sqrt(scale)
pareto <- sweep(
  sweep(logX, 1, row_center, "-"),
  1, sqrt(row_scale), "/"
)

# Sanity check
if (!all(is.finite(pareto))) {
  stop("Non-finite values after Pareto scaling.")
}

# -------------------- Save results ------------------------------

# QC summary
norm_summary <- data.frame(
  metric = c("n_features_in", "n_samples_in",
             "pareto_mode", "eps_mode", "eps_value"),
  value  = c(nrow(X), ncol(X), PARETO_MODE, EPS_MODE, eps)
)

write.csv(
  norm_summary,
  file.path(QC_DIR, "normalization_summary.csv"),
  row.names = FALSE
)

# Save run parameters (useful for reproducibility)
write.csv(
  data.frame(
    param = c("EPS_MODE","EPS_FIXED","PARETO_MODE"),
    value = c(EPS_MODE, EPS_FIXED, PARETO_MODE)
  ),
  file.path(QC_DIR, "normalization_params.csv"),
  row.names = FALSE
)

# Keep original sample order (just in case)
pareto <- pareto[, colnames(X), drop = FALSE]

# Save final normalized table
write.csv(pareto, OUT_CSV, row.names = TRUE)
cat("DONE normalization → saved:", OUT_CSV, "\n")

# Session info for reproducibility
writeLines(
  capture.output(sessionInfo()),
  file.path(LOGS_DIR, paste0("sessionInfo_normalization_", Sys.Date(), ".txt"))
)

# Tip for PCA on this matrix:
message("PCA tip: prcomp(t(pareto), center = FALSE, scale. = FALSE)")