
# -------------------------------------------------------------
# Purpose:
#   Identify potential outlier samples in normalized NEG LC–MS
#   data using:
#     - hard cut on PC1 / PC2 scores
#     - Mahalanobis distance in PCA space
#
# Input:
#   data/processed/feature_table_neg_processed_for_stats.csv
#   data/metadata/sample_sheet_neg.csv
#     - must contain 'file' or 'file_name'
#
# Output:
#   results/qc/outlier_check_neg/outlier_flags.csv
#   results/qc/outlier_check_neg/pca_outliers.png
#   results/qc/outlier_check_neg/pca_scree.png
# -------------------------------------------------------------

source("scripts/preprocessing/00_setup.R")

# -------------------- INPUT / OUTPUT ---------------------------
mat_path  <- here::here("data",  "processed",
                        "feature_table_neg_processed_for_stats.csv")
ss_path   <- here::here("data",  "metadata", "sample_sheet_neg.csv")

qc_dir    <- here::here("results", "qc", "outlier_check_neg")
out_flags <- file.path(qc_dir, "outlier_flags.csv")
out_pca   <- file.path(qc_dir, "pca_outliers.png")
out_scree <- file.path(qc_dir, "pca_scree.png")

dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------- PARAMETERS (easy to change) --------------
# Rule 1: simple cut on PC scores (similar to manual inspection)
PC1_LIMIT <- 60     # flag if |PC1| > 60
PC2_LIMIT <- 40     # flag if |PC2| > 40

# Rule 2: distance in PC space (Mahalanobis on first K PCs)
N_PCS_FOR_DISTANCE <- 5
ALPHA_DISTANCE     <- 0.001  # ~ very extreme (0.1% tail)

# -------------------- LOAD DATA --------------------------------
X <- read.csv(mat_path, row.names = 1, check.names = FALSE)
X[] <- lapply(X, function(x) suppressWarnings(as.numeric(x)))
X <- as.matrix(X)

ss <- read.csv(ss_path, check.names = FALSE)

# Helper: normalize file names so they match matrix columns
norm_names <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\","/", x)
  x <- sub(".*/", "", x)
  tolower(
    sub("\\.(raw|mzml|mzxml|d|csv|txt)$","", x, ignore.case = TRUE)
  )
}

# Build comparable keys and align orders ------------------------
# sample sheet may contain 'file_name' or 'file'
if ("file_name" %in% names(ss)) {
  file_col <- "file_name"
} else if ("file" %in% names(ss)) {
  file_col <- "file"
} else {
  stop("Neither 'file_name' nor 'file' column found in sample_sheet_neg.csv")
}

ss$sample_key <- norm_names(ss[[file_col]])
x_cols_key    <- norm_names(colnames(X))

keep <- base::intersect(x_cols_key, ss$sample_key)

X2  <- X[, match(keep, x_cols_key), drop = FALSE]
ss2 <- ss[match(keep, ss$sample_key), , drop = FALSE]
stopifnot(ncol(X2) == nrow(ss2))

# -------------------- PCA (same settings as before) ------------
# IMPORTANT: matrix is already normalized & scaled (Pareto),
# so no extra centering/scaling
pc <- prcomp(t(X2), center = FALSE, scale. = FALSE)

scores <- as.data.frame(pc$x)
expl   <- (pc$sdev^2) / sum(pc$sdev^2)

# -------------------- OUTLIER RULE 1: PC1 / PC2 cut ------------
flag_rule1 <- (abs(scores$PC1) > PC1_LIMIT) | (abs(scores$PC2) > PC2_LIMIT)

# -------------------- OUTLIER RULE 2: distance in PC space -----
# Use first K PCs (or as many as available)
K <- min(N_PCS_FOR_DISTANCE, ncol(scores))
S <- as.matrix(scores[, seq_len(K), drop = FALSE])

# Mahalanobis distance on normalized PCA scores
# If covariance is singular, add tiny jitter to diagonal
covS <- cov(S)
if (any(!is.finite(covS)) || det(covS) == 0) {
  diag(covS) <- diag(covS) + 1e-8
}
md <- mahalanobis(S, center = colMeans(S), cov = covS)

# Cutoff from chi-square with K degrees of freedom
cut_md    <- qchisq(1 - ALPHA_DISTANCE, df = K, lower.tail = TRUE)
flag_rule2 <- md > cut_md

# -------------------- COLLECT FLAGS ----------------------------
flags_df <- data.frame(
  sample_file = colnames(X2),
  PC1         = round(scores$PC1, 3),
  PC2         = round(scores$PC2, 3),
  MD          = round(md,        3),
  rule_PC_cut  = flag_rule1,
  rule_MD_tail = flag_rule2,
  any_flag     = flag_rule1 | flag_rule2,
  stringsAsFactors = FALSE
)

# Attach simple group if present (from previous steps)
if ("group" %in% names(ss2)) {
  flags_df$group <- ss2$group
}

# Save CSV
write.csv(flags_df, out_flags, row.names = FALSE)

# -------------------- PLOTS ------------------------------------
png(out_pca, width = 900, height = 650, res = 120)
par(mar = c(4,4,1,1))

grp  <- factor(if ("group" %in% names(ss2)) ss2$group else "All")
cols <- as.integer(grp) + 1   # simple palette
pch_all <- ifelse(flags_df$any_flag, 17, 19)  # triangles for flagged

plot(scores$PC1, scores$PC2,
     col = cols, pch = pch_all,
     xlab = paste0("PC1 (", round(100 * expl[1], 1), "%)"),
     ylab = paste0("PC2 (", round(100 * expl[2], 1), "%)"),
     main = "PCA with outlier flags")

abline(v = c(-PC1_LIMIT, PC1_LIMIT), lty = 3)
abline(h = c(-PC2_LIMIT, PC2_LIMIT), lty = 3)

# Label only flagged points (to avoid clutter)
lab_idx <- which(flags_df$any_flag)
if (length(lab_idx) > 0) {
  text(scores$PC1[lab_idx], scores$PC2[lab_idx],
       labels = lab_idx, pos = 3, cex = 0.8)
  legend("topright",
         legend = c(levels(grp), "flagged (triangle)"),
         col    = c(seq_len(nlevels(grp)) + 1, "black"),
         pch    = c(rep(19, nlevels(grp)), 17),
         bty    = "n")
} else {
  legend("topright", legend = levels(grp),
         col = seq_len(nlevels(grp)) + 1, pch = 19, bty = "n")
}
dev.off()

png(out_scree, width = 700, height = 450, res = 120)
barplot(100 * expl[1:min(10, length(expl))],
        xlab = "PC", ylab = "Explained variance (%)",
        main = "Scree (first 10 PCs)")
dev.off()

cat("QC outlier check saved:\n",
    " - Flags table: ", out_flags, "\n",
    " - PCA plot   : ", out_pca, "\n",
    " - Scree plot : ", out_scree, "\n", sep = "")

