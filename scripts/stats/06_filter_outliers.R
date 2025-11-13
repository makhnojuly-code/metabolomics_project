
source("scripts/preprocessing/00_setup.R")



# ======================================================
# Remove flagged and technical outlier samples,
# align metadata, and perform PCA after filtering.
# ======================================================

# --- Paths
mat_path    <- here::here("data","processed","feature_table_neg_processed_for_stats.csv")
ss_path     <- here::here("data","metadata","sample_sheet_neg.csv")
flags_path  <- here::here("results","qc","outlier_flags.csv")

out_mat     <- here::here("data","processed","feature_table_neg_processed_for_stats_filtered.csv")
out_ss      <- here::here("data","processed","sample_sheet_neg_filtered.csv")
out_removed <- here::here("results","qc","removed_samples.csv")
out_pca_png <- here::here("results","qc","pca_after_filter_colored.png")
out_scree   <- here::here("results","qc","pca_after_filter_scree.png")

dir.create(here::here("results","qc"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("data","processed"), recursive = TRUE, showWarnings = FALSE)

#  Helper function to normalize sample names
# Converts file paths to lowercase clean identifiers
norm_names <- function(x){
  x <- as.character(x)
  x <- gsub("\\\\","/", x)
  x <- sub(".*/", "", x)
  tolower(sub("\\.(raw|mzml|mzxml|d|csv|txt)$","", x, ignore.case = TRUE))
}

#  Load input data
stopifnot(file.exists(mat_path), file.exists(flags_path), file.exists(ss_path))

X <- read.csv(mat_path, row.names = 1, check.names = FALSE)
X[] <- lapply(X, function(x) suppressWarnings(as.numeric(x)))  # Ensure numeric values
X <- as.matrix(X)

ss <- read.csv(ss_path, check.names = FALSE)
stopifnot("file_name" %in% names(ss))

flags <- read.csv(flags_path, check.names = FALSE)
stopifnot(all(c("sample_file","any_flag") %in% names(flags)))

# Align data by sample keys (use explicit base:: to avoid conflicts)
ss$sample_key     <- norm_names(ss$file_name)
x_cols_key        <- norm_names(colnames(X))
flags$sample_key  <- norm_names(flags$sample_file)

keep_keys   <- base::intersect(x_cols_key, ss$sample_key)
X2          <- X[, base::match(keep_keys, x_cols_key), drop = FALSE]
ss2         <- ss[base::match(keep_keys, ss$sample_key), , drop = FALSE]
flags2      <- flags[base::match(keep_keys, flags$sample_key), , drop = FALSE]
stopifnot(ncol(X2) == nrow(ss2), ncol(X2) == nrow(flags2))

#  Filter out flagged or technical samples (QC, Blank, MSMS, UE)
tech_re <- "(blank|qc|mm8|msms|ue)"
rm_flag  <- isTRUE(flags2$any_flag) | grepl(tech_re, ss2$sample_key, ignore.case = TRUE)
rm_flag[is.na(rm_flag)] <- FALSE

#  Prepare dataframe of removed samples for record keeping
removed_df <- data.frame(
  sample_file = colnames(X2)[rm_flag],
  sample_key  = ss2$sample_key[rm_flag],
  reason      = ifelse(isTRUE(flags2$any_flag[rm_flag]), "any_flag=TRUE", "pattern_tech"),
  stringsAsFactors = FALSE
)

# Keep only clean samples
Xf  <- X2[, !rm_flag, drop = FALSE]
ssf <- ss2[!rm_flag, , drop = FALSE]

#  Add group column if missing (based on file naming)
if (!("group" %in% names(ssf))) {
  ssf$group <- ifelse(grepl("null",     ssf$sample_key, ignore.case = TRUE), "Null",
                      ifelse(grepl("wounding", ssf$sample_key, ignore.case = TRUE), "Wounding", "Other"))
}

#  Save cleaned outputs
write.csv(Xf, out_mat, row.names = TRUE)
write.csv(ssf[, base::intersect(c("file_name","sample_key","group"), names(ssf))],
          out_ss, row.names = FALSE)
write.csv(removed_df, out_removed, row.names = FALSE)

cat("Original samples:", ncol(X),  "\n")
cat("Kept after filter:", ncol(Xf), "\n")
cat("Removed samples   :", nrow(removed_df), "\n\n")

#  PCA on filtered data
pc  <- prcomp(t(Xf), center = FALSE, scale. = FALSE)
sc  <- as.data.frame(pc$x[, 1:2, drop = FALSE])
ve  <- (pc$sdev^2) / sum(pc$sdev^2)

grp <- factor(ssf$group)
cols <- as.integer(grp) + 1  # Assign colors by group

#  PCA plot (colored by group)
png(out_pca_png, width = 900, height = 650, res = 120)
par(mar = c(4,4,2,1))
plot(sc$PC1, sc$PC2, col = cols, pch = 19,
     xlab = paste0("PC1 (", round(100*ve[1], 1), "%)"),
     ylab = paste0("PC2 (", round(100*ve[2], 1), "%)"),
     main = "PCA after filtering (colored by group)")
legend("topright", legend = levels(grp),
       col = seq_len(nlevels(grp)) + 1, pch = 19, bty = "n")
dev.off()

# Scree plot (first 10 PCs)
png(out_scree, width = 700, height = 450, res = 120)
barplot(100*ve[1:min(10,length(ve))],
        xlab = "PC", ylab = "Explained variance (%)",
        main = "Scree (first 10 PCs)")
dev.off()

cat("PCA plot saved to:", out_pca_png, "\n")
cat("Scree plot saved to:", out_scree, "\n")



X <- read.csv("data/processed/feature_table_neg_processed_for_stats_filtered.csv", row.names = 1)
hist(X[,1], main="Distribution of intensities (example feature)", xlab="Intensity (log2, PQN Pareto)")


mean_values <- apply(X, 1, mean, na.rm = TRUE)
hist(mean_values,
     main = "Average intensity across all features",
     xlab = "Mean intensity (log2, PQN Pareto)",
     col = "lightgray", border = "black")