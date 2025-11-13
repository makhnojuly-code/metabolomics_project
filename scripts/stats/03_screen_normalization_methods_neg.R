###############################################
# 
# Goal:
#   Universal screening of normalization / transformation /
#   scaling methods on NEG LC–MS feature table.
#
#   Methods included:
#     1) raw
#     2) log2
#     3) log2_pareto
#     4) log2_auto
#     5) log2_robust
#     6) median_log2_pareto
#     7) pqn_log2_pareto
#     8) arcsinh_pareto
#
#   For each method:
#     - apply normalization / transform / scaling
#     - run PCA (PC1 vs PC2)
#     - calculate a PCA-based separation score (Null vs Wounding)
#     - save a PCA plot
#     - collect scores into a summary table
###############################################

source("scripts/preprocessing/00_setup.R")

library(tidyverse)
library(ggplot2)

# ---------- Paths ------------------------------------------------

feat_path   <- here("data", "processed", "feature_table_neg_clean.csv")
sample_path <- here("data", "processed", "sample_sheet_neg_clean.csv")
out_dir     <- here("results", "qc", "screen_normalization_methods_neg")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 1. Load feature table & sample sheet -----------------

feature_df  <- read.csv(feat_path, check.names = FALSE)
feature_mat <- feature_df %>% column_to_rownames("feature_id")

sample_sheet <- read.csv(sample_path)

cat("Features:", nrow(feature_mat), "\n")
cat("Samples :", ncol(feature_mat), "\n")

# ---------- 2. Match sample IDs & keep biological samples --------

# Map column names in the feature table to sample_key in the sample sheet
col_annot <- tibble(
  colname    = colnames(feature_mat),
  sample_key = colnames(feature_mat) %>%
    tolower() %>%
    gsub("\\.mzml$", "", ., ignore.case = TRUE)
)

# Keep only biological samples (Null, Wounding)
sample_bio <- sample_sheet %>%
  filter(group %in% c("Null", "Wounding")) %>%
  filter(!is.na(sample_key))

annot_join <- col_annot %>%
  inner_join(sample_bio, by = "sample_key")

# Subset feature table to matched biological samples
feature_bio <- feature_mat[, annot_join$colname, drop = FALSE]

# IMPORTANT: coerce to numeric matrix (avoid data.frame/list issues)
feature_bio[] <- lapply(feature_bio, function(x) suppressWarnings(as.numeric(x)))
feature_bio <- as.matrix(feature_bio)

meta <- annot_join %>%
  select(sample_key, group, colname)

stopifnot(ncol(feature_bio) == nrow(meta))

# ---------- 3. Basic building blocks -----------------------------

# Sample-wise median normalization:
# divide each sample by its median intensity.
median_normalize_samples <- function(mat) {
  mat <- as.matrix(mat)
  sample_medians <- apply(mat, 2, median, na.rm = TRUE)
  sample_medians[!is.finite(sample_medians) | sample_medians == 0] <- 1
  sweep(mat, 2, sample_medians, "/")
}

# Sample-wise PQN normalization:
# 1) normalize features to their median across samples
# 2) compute a median-based dilution factor per sample
# 3) divide each sample by its PQN factor
pqn_normalize_samples <- function(mat) {
  mat <- as.matrix(mat)
  ref <- apply(mat, 1, median, na.rm = TRUE)
  ref[!is.finite(ref) | ref <= 0] <- 1e-12
  Xn <- sweep(mat, 1, ref, "/")
  pqn_factor <- apply(Xn, 2, function(v) median(v[is.finite(v)], na.rm = TRUE))
  pqn_factor[!is.finite(pqn_factor) | pqn_factor == 0] <- 1
  sweep(Xn, 2, pqn_factor, "/")
}

# Log2 transformation with small epsilon to avoid log(0)
transform_log2 <- function(mat) {
  mat <- as.matrix(mat)
  mat[!is.finite(mat) | mat <= 0] <- 1e-12
  log2(mat)
}

# Arcsinh transformation (alternative variance-stabilizing transform)
transform_arcsinh <- function(mat) {
  mat <- as.matrix(mat)
  asinh(mat)
}

# Feature-wise Pareto scaling:
# (value - mean) / sqrt(sd)
scale_pareto <- function(mat) {
  mat <- as.matrix(mat)
  sds <- apply(mat, 1, sd, na.rm = TRUE)
  sds[!is.finite(sds) | sds == 0] <- 1
  sweep(mat, 1, sqrt(sds), "/")
}

# Feature-wise autoscaling (Z-score):
# (value - mean) / sd
scale_auto <- function(mat) {
  mat <- as.matrix(mat)
  means <- rowMeans(mat, na.rm = TRUE)
  sds   <- apply(mat, 1, sd, na.rm = TRUE)
  sds[!is.finite(sds) | sds == 0] <- 1
  x <- sweep(mat, 1, means, "-")
  sweep(x, 1, sds, "/")
}

# Feature-wise robust scaling:
# (value - median) / MAD
scale_robust <- function(mat) {
  mat <- as.matrix(mat)
  med <- apply(mat, 1, median, na.rm = TRUE)
  mad <- apply(mat, 1, mad,    na.rm = TRUE, constant = 1)
  mad[!is.finite(mad) | mad == 0] <- 1
  x <- sweep(mat, 1, med, "-")
  sweep(x, 1, mad, "/")
}

# ---------- 4. Master normalizer: combine steps by method --------

normalize_matrix <- function(mat, method) {
  x <- as.matrix(mat)
  
  if (method == "raw") {
    return(x)
  }
  
  if (method == "log2") {
    x <- transform_log2(x)
    return(x)
  }
  
  if (method == "log2_pareto") {
    x <- transform_log2(x)
    x <- scale_pareto(x)
    return(x)
  }
  
  if (method == "log2_auto") {
    x <- transform_log2(x)
    x <- scale_auto(x)
    return(x)
  }
  
  if (method == "log2_robust") {
    x <- transform_log2(x)
    x <- scale_robust(x)
    return(x)
  }
  
  if (method == "median_log2_pareto") {
    x <- median_normalize_samples(x)
    x <- transform_log2(x)
    x <- scale_pareto(x)
    return(x)
  }
  
  if (method == "pqn_log2_pareto") {
    x <- pqn_normalize_samples(x)
    x <- transform_log2(x)
    x <- scale_pareto(x)
    return(x)
  }
  
  if (method == "arcsinh_pareto") {
    x <- transform_arcsinh(x)
    x <- scale_pareto(x)
    return(x)
  }
  
  stop("Unknown method: ", method)
}

# ---------- 5. PCA-based separation score ------------------------

# Compute how well PCA separates Null vs Wounding in PC1/PC2.
# Score = distance between group centroids / average within-group variance.
evaluate_pca_separation <- function(mat_norm, meta) {
  mat_norm <- as.matrix(mat_norm)
  
  # prcomp expects rows = observations (samples), columns = variables (features)
  # here we have rows = features, columns = samples → transpose
  pca <- prcomp(t(mat_norm), center = TRUE, scale. = FALSE)
  scores <- as.data.frame(pca$x[, 1:2])
  scores$group <- meta$group
  
  # group centroids in PC1–PC2 space
  centroids <- scores %>%
    group_by(group) %>%
    summarise(across(starts_with("PC"), mean), .groups = "drop")
  
  if (nrow(centroids) < 2) {
    return(list(score = NA, scores = scores, pca = pca))
  }
  
  # Euclidean distance between group centroids
  dist_centroids <- sqrt(
    (centroids$PC1[1] - centroids$PC1[2])^2 +
      (centroids$PC2[1] - centroids$PC2[2])^2
  )
  
  # Average within-group variance on PC1 & PC2
  within_var <- scores %>%
    group_by(group) %>%
    summarise(across(starts_with("PC"), var), .groups = "drop") %>%
    summarise(across(starts_with("PC"), mean), .groups = "drop") %>%
    unlist() %>%
    mean(na.rm = TRUE)
  
  score <- dist_centroids / (within_var + 1e-8)
  
  list(score = score, scores = scores, pca = pca)
}

# ---------- 6. Loop over methods --------------------------------

methods <- c(
  "raw",
  "log2",
  "log2_pareto",
  "log2_auto",
  "log2_robust",
  "median_log2_pareto",
  "pqn_log2_pareto",
  "arcsinh_pareto"
)

scores <- numeric(length(methods))
names(scores) <- methods

for (m in methods) {
  message("Processing method: ", m)
  
  mat_norm <- normalize_matrix(feature_bio, m)
  eval_res <- evaluate_pca_separation(mat_norm, meta)
  
  scores[m] <- eval_res$score
  
  # PCA plot for this method
  p <- ggplot(eval_res$scores,
              aes(x = PC1, y = PC2, colour = group)) +
    geom_point(size = 2.5, alpha = 0.8) +
    theme_bw() +
    labs(
      title    = paste("PCA - method:", m),
      subtitle = "Null vs Wounding",
      x        = "PC1",
      y        = "PC2"
    )
  
  ggsave(
    file.path(out_dir, paste0("neg_pca_", m, ".png")),
    p, width = 8, height = 6, dpi = 300
  )
}

# ---------- 7. Summary table ------------------------------------

score_tbl <- tibble(
  method = methods,
  pca_separation_score = as.numeric(scores)
) %>%
  arrange(desc(pca_separation_score))

readr::write_csv(
  score_tbl,
  file.path(out_dir, "neg_normalization_method_scores.csv")
)

print(score_tbl)
cat("\n[OK] Normalization method screening finished.\nResults in:\n", out_dir, "\n")