source("scripts/preprocessing/00_setup.R")

# PCA after normalization (simple, robust)
# Input:
#   data/processed/feature_table_neg_processed_for_stats.csv
#   data/metadata/sample_sheet_neg.csv   (columns: file, ... )
# Output:
#   data/processed/sample_sheet_neg_clean.csv
#   results/qc/pca_after_norm.png

# Paths 
mat_path <- here::here("data","processed","feature_table_neg_processed_for_stats.csv")
ss_path  <- here::here("data","metadata","sample_sheet_neg.csv")
out_ss   <- here::here("data","processed","sample_sheet_neg_clean.csv")
out_png  <- here::here("results","pca_after_norm.png")

#  Load matrix (samples = columns) 
X <- read.csv(mat_path, row.names = 1, check.names = FALSE)
X[] <- lapply(X, function(x) suppressWarnings(as.numeric(x)))
X <- as.matrix(X)

# Load sample sheet 
ss <- read.csv(ss_path, check.names = FALSE)

#  Helper: normalize names (drop path + extension, case-insensitive) 
norm_names <- function(x){
  x <- as.character(x)
  x <- gsub("\\\\","/", x)
  x <- sub(".*/", "", x)
  x <- tolower(sub("\\.(raw|mzml|mzxml|d|csv|txt)$","", x, ignore.case=TRUE))
  x
}

# Build comparable name keys
stopifnot("file_name" %in% names(ss))
ss$sample_key <- norm_names(ss$file)
x_cols_key    <- norm_names(colnames(X))

# Drop technical injections by name pattern (both in sample sheet and X) 
tech_patterns <- c("blank","acn","mm8","qc")
tech_re <- paste(tech_patterns, collapse="|")
ss <- ss[ !grepl(tech_re, ss$sample_key, ignore.case = TRUE), , drop = FALSE]

# keep only intersection & align orders (X follows its current column order)
keep <- base::intersect(x_cols_key, ss$sample_key)
X2   <- X[, match(keep, x_cols_key), drop = FALSE]
ss2  <- ss[match(keep, ss$sample_key), , drop = FALSE]

# Safety: dimensions should match now
stopifnot(ncol(X2) == nrow(ss2))

#  Auto-assign simple groups from name (edit if needed) 
ss2$group <- ifelse(grepl("null",     ss2$sample_key, ignore.case = TRUE), "Null",
                    ifelse(grepl("wounding", ss2$sample_key, ignore.case = TRUE), "Wounding", NA))
if (any(is.na(ss2$group))) message("Note: some samples have group = NA.")

# Save clean sheet for stats (sample, group)
dir.create(dirname(out_ss), recursive = TRUE, showWarnings = FALSE)
write.csv(ss2[, c("sample_key","group")], out_ss, row.names = FALSE)

#  PCA (samples must be rows → transpose) 
pc  <- prcomp(t(X2), center = FALSE, scale. = FALSE)
sc  <- as.data.frame(pc$x[, 1:2, drop = FALSE])
ve  <- (pc$sdev^2) / sum(pc$sdev^2)  # variance explained
grp <- factor(ss2$group)
cols <- as.integer(grp) + 1  # simple palette

#  Plot & save 
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
png(out_png, width = 900, height = 650, res = 120)
par(mar = c(4,4,1,1))
plot(sc$PC1, sc$PC2,
     col = cols, pch = 19,
     xlab = paste0("PC1 (", round(100*ve[1], 1), "%)"),
     ylab = paste0("PC2 (", round(100*ve[2], 1), "%)"),
     main = "PCA after normalization")
if (nlevels(grp) > 0) {
  legend("topright", legend = levels(grp),
         col = seq_len(nlevels(grp)) + 1, pch = 19, bty = "n")
}
dev.off()

cat("Saved:\n  - clean sample sheet:", out_ss, "\n  - PCA plot:", out_png, "\n")
