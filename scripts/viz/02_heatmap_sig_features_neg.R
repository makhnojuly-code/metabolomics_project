# ================================================================
# 02_heatmap_sig_features_neg.R
# Heatmap of significant features (NEG mode, Wounding vs Null)
#
# This script produces a clustered heatmap of the top statistically
# and biologically relevant features identified by limma.
# The steps below describe the complete workflow for documentation.
#
# ------------------------- WHAT THIS SCRIPT DOES -------------------------
#
# Step 1–3 – Load input tables:
#   • the expression table (feature × sample intensity values)
#   • the sample metadata (sample names, experimental groups)
#   • the limma result table (logFC, p-values, FDR)
#
# Step 4 – Convert the feature table into a numeric expression matrix:
#   • rows = features (FT IDs)
#   • columns = samples
#
# Step 5 – Ensure that the column order of the expression matrix
#   exactly matches the row order of the metadata (via `file_name`).
#   This alignment is critical for all downstream analyses.
#
# Step 6 – Select biologically and statistically relevant features:
#   • FDR-adjusted p-value < fdr_cutoff (default 0.05)
#   • |log2FC| > logfc_cutoff (default 1, i.e. ≥ 2-fold change)
#   • keep only the first max_features (default 50) to avoid
#     overcrowding the heatmap and ensure interpretability
#
# Step 7 – Build a column annotation dataframe containing the sample group
#   (Null / Wounding), which is shown as a colored bar above the heatmap.
#
# Step 8 – Define a diverging color palette:
#   • blue  = low abundance
#   • white = medium abundance
#   • red   = high abundance
#
# Step 9 – Draw a clustered heatmap:
#   • hierarchical clustering of rows and columns (Euclidean distance)
#   • intensities scaled per feature (`scale = "row"`)
#     → highlights *patterns* rather than absolute intensities
#   • final heatmap saved to: results/viz/heatmap_sig_features_neg.png
#
# ========================================================================


# --- 0. Load required packages ---------------------------------

source("scripts/preprocessing/00_setup.R")

library(tidyverse)
library(pheatmap)

message("Packages for heatmap loaded.")


# --- 1. Define input/output paths ------------------------------

expr_path <- "data/processed/feature_table_neg_processed_for_stats_filtered.csv"
meta_path <- "data/processed/sample_sheet_neg_filtered.csv"
res_path  <- "results/limma_results_neg.csv"

out_dir  <- "results/viz"
out_file <- file.path(out_dir, "heatmap_sig_features_neg.png")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# --- 2. Analysis parameters ------------------------------------

fdr_cutoff    <- 0.05
logfc_cutoff  <- 1.0     # 2-fold change
max_features  <- 20

# --- 3. Load data ----------------------------------------------

feat_mat <- read.csv(expr_path, check.names = FALSE)
meta     <- read.csv(meta_path, check.names = FALSE)
res      <- read.csv(res_path, row.names = 1, check.names = FALSE)

message("Input tables loaded.")


# --- 4. Build expression matrix -------------------------------

feature_ids <- feat_mat[[1]]
expr <- as.matrix(feat_mat[, -1])  # drop FT ID column
rownames(expr) <- feature_ids


# --- 5. Align expression matrix with metadata -----------------

stopifnot("file_name" %in% names(meta))

if (!all(colnames(expr) == meta$file_name)) {
  expr <- expr[, meta$file_name]
}

stopifnot(all(colnames(expr) == meta$file_name))

message("Expression matrix and metadata are aligned.")


# --- 6. Select significant features ----------------------------

sig_feats <- res %>%
  dplyr::filter(adj.P.Val < fdr_cutoff & abs(logFC) > logfc_cutoff) %>%
  dplyr::arrange(adj.P.Val)

if (nrow(sig_feats) == 0) {
  stop("No significant features with current cutoffs.")
}

if (nrow(sig_feats) > max_features) {
  sig_feats <- head(sig_feats, max_features)
}

sel_ids  <- rownames(sig_feats)
heat_mat <- expr[sel_ids, , drop = FALSE]

message("Number of features in heatmap: ", nrow(heat_mat))


# --- 7. Annotation: sample groups ------------------------------

annotation_col <- data.frame(
  Group = factor(meta$group)
)
rownames(annotation_col) <- meta$file_name


# --- 8. Color palette ------------------------------------------

palette_length <- 100
my_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(palette_length)


# --- 9. Draw and save heatmap ---------------------------------

pheatmap(
  mat       = heat_mat,
  color     = my_colors,
  scale     = "row",
  breaks = seq(-2.5, 2.5, length.out = palette_length),
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "complete",
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = annotation_col,
  main = paste0("Top ", nrow(heat_mat),
                " significant features (FDR < ", fdr_cutoff,
                ", |log2FC| > ", logfc_cutoff, ")"),
  filename = out_file,
  width = 8,
  height = 6
)

message("Heatmap saved to: ", out_file)




