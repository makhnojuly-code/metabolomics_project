# -------------------------------------------------------------
# 
# Purpose:
#   Reproduce Khabat’s “before normalization” step for NEG data:
#     - look at cleaned raw data BEFORE any normalization
#     - density / histogram
#     - boxplot across samples
#     - PCA scatter plot (PC1 vs PC2)
#     - help decide:
#         * are there strong outliers?
#         * is the distribution approximately normal (after log)?
#         * is there visible group separation (Null vs Wounding)?
#         * is normalization needed, and what type?
#
# Input:
#   data/processed/feature_table_neg_clean.csv   (features x samples)
#   data/processed/sample_sheet_neg_clean.csv    (sample_key, group, ...)
#
# Output (QC plots):
#   results/qc/before_normalization_neg/
#     - neg_raw_density_by_feature.png
#     - neg_raw_density_overall_log10.png
#     - neg_raw_boxplot_samples_log10.png
#     - neg_raw_pca_pc1_pc2_clean.png
#     - neg_raw_pca_pc1_pc2_labels.png
# -------------------------------------------------------------

source("scripts/preprocessing/00_setup.R")


# ---------- Paths ------------------------------------------------

feat_path   <- here("data", "processed", "feature_table_neg_clean.csv")
sample_path <- here("data", "processed", "sample_sheet_neg_clean.csv")
qc_dir      <- here("results", "qc", "before_normalization_neg")

stopifnot(file.exists(feat_path), file.exists(sample_path))
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 1. Load data -----------------------------------------

# Feature table: rows = features, first col = feature_id, others = samples
feature_df <- read.csv(feat_path, check.names = FALSE)
feature_mat <- feature_df %>%
  column_to_rownames("feature_id")

cat("Features:", nrow(feature_mat), "\n")
cat("Samples :", ncol(feature_mat), "\n")

# Sample sheet: must contain sample_key + group (Null / Wounding)
sample_sheet <- read.csv(sample_path, check.names = FALSE)

# ---------- 2. Harmonise sample IDs & keep biological samples ----

# Match colnames(feature_mat) to sample_sheet$sample_key
col_annot <- tibble(
  colname    = colnames(feature_mat),
  sample_key = colnames(feature_mat) %>%
    tolower() %>%
    gsub("\\.mzml$", "", ., ignore.case = TRUE)
)

sample_bio <- sample_sheet %>%
  dplyr::filter(group %in% c("Null", "Wounding")) %>%
  dplyr::filter(!is.na(sample_key))

annot_join <- col_annot %>%
  inner_join(sample_bio, by = "sample_key")

cat("Matched biological samples:", nrow(annot_join), "\n")

feature_bio <- feature_mat[, annot_join$colname, drop = FALSE]
stopifnot(ncol(feature_bio) == nrow(annot_join))

meta <- annot_join %>%
  dplyr::select(sample_key, group, colname)

# ---------- 3. Long format for plotting --------------------------

feature_long <- feature_bio %>%
  as.data.frame() %>%
  rownames_to_column("feature_id") %>%
  pivot_longer(
    cols      = -feature_id,
    names_to  = "colname",
    values_to = "intensity"
  ) %>%
  left_join(meta, by = "colname")

# ========== A) Density / Histogram ===============================

# 1) Density per feature for a small subset (first 9 features)
n_features_for_plot <- 9
top_features <- unique(feature_long$feature_id)[1:n_features_for_plot]

feature_long_sub <- feature_long %>%
  filter(feature_id %in% top_features)

p_density_features <- ggplot(feature_long_sub,
                             aes(x = intensity, colour = group)) +
  geom_density() +
  facet_wrap(~ feature_id, scales = "free") +
  theme_bw() +
  labs(
    title    = "Raw NEG data: density per feature",
    subtitle = "Raw intensities before normalization (subset of features)",
    x        = "Intensity (raw)",
    y        = "Density"
  )

ggsave(file.path(qc_dir, "neg_raw_density_by_feature.png"),
       p_density_features, width = 10, height = 7, dpi = 300)

# 2) Overall density on log10 scale (for normality check)
p_density_overall <- ggplot(feature_long,
                            aes(x = log10(intensity + 1), colour = group)) +
  geom_density(alpha = 0.4) +
  theme_bw() +
  labs(
    title    = "Raw NEG data: overall log10 density",
    subtitle = "Used to check approximate normality after log-transform",
    x        = "log10(intensity + 1)",
    y        = "Density"
  )

ggsave(file.path(qc_dir, "neg_raw_density_overall_log10.png"),
       p_density_overall, width = 8, height = 6, dpi = 300)

# ========== B) Boxplot (sample-wise variation) ===================

# Boxplot with log10 y-scale → easier to see differences
p_box_log <- ggplot(feature_long,
                    aes(x = colname, y = intensity, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_log10() +
  coord_flip() +
  theme_bw() +
  labs(
    title    = "Raw NEG data: boxplot per sample (log10 scale)",
    subtitle = "Used to check sample-wise intensity range and potential outliers",
    x        = "Sample (column name)",
    y        = "Intensity (raw, log10 scale)"
  )

ggsave(file.path(qc_dir, "neg_raw_boxplot_samples_log10.png"),
       p_box_log, width = 10, height = 8, dpi = 300)

# ========== C) PCA (scatter plot) ================================

# For PCA: rows = samples, cols = features
pca_input <- t(feature_bio)

pca_raw <- prcomp(
  pca_input,
  center = TRUE,   # centering only
  scale. = FALSE   # no scaling / normalization
)

pca_df <- as.data.frame(pca_raw$x) %>%
  mutate(
    sample_key = meta$sample_key,
    group      = meta$group
  )

# 1) PCA without text labels (clean scatter)
p_pca_clean <- ggplot(pca_df,
                      aes(x = PC1, y = PC2, colour = group)) +
  geom_point(size = 2.5, alpha = 0.8) +
  theme_bw() +
  labs(
    title    = "PCA on raw NEG data",
    subtitle = "Null vs Wounding before any normalization",
    x        = "PC1",
    y        = "PC2"
  )

ggsave(file.path(qc_dir, "neg_raw_pca_pc1_pc2_clean.png"),
       p_pca_clean, width = 8, height = 6, dpi = 300)

# 2) PCA with text labels (may be cluttered but useful for inspection)
p_pca_labels <- ggplot(pca_df,
                       aes(x = PC1, y = PC2, colour = group, label = sample_key)) +
  geom_point(size = 2) +
  geom_text(vjust = -0.5, size = 2, show.legend = FALSE) +
  theme_bw() +
  labs(
    title    = "PCA on raw NEG data (with labels)",
    subtitle = "Only for detailed inspection; can be hard to read",
    x        = "PC1",
    y        = "PC2"
  )

ggsave(file.path(qc_dir, "neg_raw_pca_pc1_pc2_labels.png"),
       p_pca_labels, width = 10, height = 8, dpi = 300)

cat("\n[OK] Raw NEG data inspection finished.\nPlots saved in:\n", qc_dir, "\n")