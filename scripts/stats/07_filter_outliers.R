# -------------------------------------------------------------

# Purpose:
#   Remove flagged / technical outlier samples from normalized
#   NEG LC–MS data, align metadata, and check PCA after filtering.
#
# Input:
#   data/processed/feature_table_neg_processed_for_stats.csv
#   data/metadata/sample_sheet_neg.csv
#   results/qc/outlier_check_neg/outlier_flags.csv
#
# Output:
#   data/processed/feature_table_neg_processed_for_stats_filtered.csv
#   data/processed/sample_sheet_neg_filtered.csv
#
#   results/qc/after_filter_neg/removed_samples.csv
#   results/qc/after_filter_neg/pca_after_filter_colored.png
#   results/qc/after_filter_neg/pca_after_filter_scree.png
#   results/qc/after_filter_neg/hist_example_feature.png
#   results/qc/after_filter_neg/hist_mean_intensity.png
# -------------------------------------------------------------


source("scripts/preprocessing/00_setup.R")


# -------------------------- PATHS -----------------------------

mat_path   <- here::here("data","processed",
                         "feature_table_neg_processed_for_stats.csv")
ss_path    <- here::here("data","metadata","sample_sheet_neg.csv")
flags_path <- here::here("results","qc","outlier_check_neg","outlier_flags.csv")

out_mat     <- here::here("data","processed",
                          "feature_table_neg_processed_for_stats_filtered.csv")
out_ss      <- here::here("data","processed","sample_sheet_neg_filtered.csv")

qc_dir      <- here::here("results","qc","after_filter_neg")
out_removed <- file.path(qc_dir, "removed_samples.csv")
out_pca_png <- file.path(qc_dir, "pca_after_filter_colored.png")
out_scree   <- file.path(qc_dir, "pca_after_filter_scree.png")
out_hist1   <- file.path(qc_dir, "hist_example_feature.png")
out_hist2   <- file.path(qc_dir, "hist_mean_intensity.png")
out_info    <- file.path(qc_dir, "median_feature_used.csv")

dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------- Helper: clean sample names ---------------------

norm_names <- function(x){
  x <- as.character(x)
  x <- gsub("\\\\","/", x)
  x <- sub(".*/", "", x)
  tolower(sub("\\.(raw|mzml|mzxml|d|csv|txt)$","", x, ignore.case = TRUE))
}

# -------------------- Load input data --------------------

stopifnot(file.exists(mat_path), file.exists(flags_path), file.exists(ss_path))

X <- read.csv(mat_path, row.names = 1, check.names = FALSE)
X[] <- lapply(X, function(x) suppressWarnings(as.numeric(x)))
X <- as.matrix(X)

ss <- read.csv(ss_path, check.names = FALSE)
flags <- read.csv(flags_path, check.names = FALSE)

stopifnot("file_name" %in% names(ss))
stopifnot(all(c("sample_file","any_flag") %in% names(flags)))

# Align names
ss$sample_key    <- norm_names(ss$file_name)
x_cols_key       <- norm_names(colnames(X))
flags$sample_key <- norm_names(flags$sample_file)

keep_keys <- base::intersect(x_cols_key, ss$sample_key)

X2     <- X[, match(keep_keys, x_cols_key), drop = FALSE]
ss2    <- ss[match(keep_keys, ss$sample_key), , drop = FALSE]
flags2 <- flags[match(keep_keys, flags$sample_key), , drop = FALSE]

stopifnot(ncol(X2) == nrow(ss2))
stopifnot(nrow(flags2) == ncol(X2))

# --------------------- Filtering rules ---------------------

tech_re <- "(blank|qc|mm8|msms|ue)"
rm_flag <- (flags2$any_flag == TRUE) |
  grepl(tech_re, ss2$sample_key, ignore.case = TRUE)

rm_flag[is.na(rm_flag)] <- FALSE

removed_df <- data.frame(
  sample_file = colnames(X2)[rm_flag],
  sample_key  = ss2$sample_key[rm_flag],
  reason      = ifelse(flags2$any_flag[rm_flag], "any_flag=TRUE", "pattern_tech"),
  stringsAsFactors = FALSE
)

# Keep clean samples
Xf  <- X2[, !rm_flag, drop = FALSE]
ssf <- ss2[!rm_flag, , drop = FALSE]

# Add group if missing
if (!"group" %in% names(ssf)) {
  ssf$group <- ifelse(grepl("null", ssf$sample_key, ignore.case = TRUE), "Null",
                      ifelse(grepl("wounding", ssf$sample_key, ignore.case = TRUE),
                             "Wounding", "Other"))
}

# Save results
write.csv(Xf, out_mat, row.names = TRUE)
write.csv(ssf[, c("file_name","sample_key","group")],
          out_ss, row.names = FALSE)
write.csv(removed_df, out_removed, row.names = FALSE)

cat("Original samples:", ncol(X), "\n")
cat("Kept after filter:", ncol(Xf), "\n")
cat("Removed samples   :", nrow(removed_df), "\n\n")


# ------------------------- PCA after filtering -------------------------

pc <- prcomp(t(Xf), center = FALSE, scale. = FALSE)
sc <- as.data.frame(pc$x[, 1:2, drop = FALSE])
ve <- (pc$sdev^2) / sum(pc$sdev^2)

grp <- factor(ssf$group)
cols <- as.integer(grp) + 1

png(out_pca_png, width = 900, height = 650, res = 120)
par(mar = c(4,4,2,1))
plot(sc$PC1, sc$PC2,
     col = cols, pch = 19,
     xlab = paste0("PC1 (", round(100*ve[1], 1), "%)"),
     ylab = paste0("PC2 (", round(100*ve[2], 1), "%)"),
     main = "PCA after filtering (colored by group)")
legend("topright", legend = levels(grp),
       col = seq_len(nlevels(grp)) + 1, pch = 19, bty = "n")
dev.off()

png(out_scree, width = 700, height = 450, res = 120)
barplot(100 * ve[1:min(10, length(ve))],
        xlab = "PC", ylab = "Explained variance (%)",
        main = "Scree (first 10 PCs)")
dev.off()


# ------------------------- Histograms -------------------------

# 1) Choose median feature (robust)
feature_medians <- apply(Xf, 1, median, na.rm = TRUE)
feature_idx <- which.min(abs(feature_medians - median(feature_medians)))
median_feature_name <- rownames(Xf)[feature_idx]

df_feat <- data.frame(intensity = Xf[feature_idx, ])

median_info <- data.frame(
  feature = median_feature_name,
  q05     = quantile(df_feat$intensity, 0.05),
  q50     = quantile(df_feat$intensity, 0.50),
  q95     = quantile(df_feat$intensity, 0.95)
)

write.csv(median_info, out_info, row.names = FALSE)

p_hist1 <- ggplot(df_feat, aes(x = intensity)) +
  geom_histogram(bins = 30, fill = "grey80", colour = "black") +
  geom_vline(xintercept = median_info$q50,
             linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = c(median_info$q05, median_info$q95),
             linetype = "dotted", linewidth = 0.7) +
  theme_bw(base_size = 14) +
  labs(
    title = "Distribution of intensities (median feature)",
    subtitle = paste("Feature:", median_feature_name),
    x = "Intensity (log2, PQN Pareto)",
    y = "Frequency"
  )

ggsave(out_hist1, p_hist1, width = 8, height = 6, dpi = 300)


# 2) Histogram of mean intensity across all features
mean_values <- apply(Xf, 1, mean, na.rm = TRUE)
df_mean <- data.frame(mean_intensity = mean_values)

median_mean <- median(mean_values)
q_mean <- quantile(mean_values, c(0.05, 0.95))

p_hist2 <- ggplot(df_mean, aes(x = mean_intensity)) +
  geom_histogram(bins = 30, fill = "grey80", colour = "black") +
  geom_vline(xintercept = median_mean,
             linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = q_mean,
             linetype = "dotted", linewidth = 0.7) +
  theme_bw(base_size = 14) +
  labs(
    title = "Average intensity across all features",
    x = "Mean intensity (log2, PQN Pareto)",
    y = "Frequency"
  )

ggsave(out_hist2, p_hist2, width = 8, height = 6, dpi = 300)

cat("QC histograms saved to:", qc_dir, "\n")
