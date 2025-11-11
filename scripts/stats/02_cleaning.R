source("scripts/preprocessing/00_setup.R")

# 1) Load LC–MS feature table (negative mode)
# 2) Remove technical samples (QC / MM8 / ACN blanks)
# 3) Check missing values (NAs) and save simple QC plots
# 4) Keep features with <=30% NAs
# 5) Replace remaining NAs with half of row minimum (half‑min)
# 6) Save cleaned table and summaries
# 7) Save session info for reproducibility

# • Technical samples are not biological → exclude them from statistics.
# • Too many NAs make features unreliable → filter by threshold.
# • NA ≈ “below detection limit” in LC–MS → half‑min keeps correct scale.
# • QC plots + logs = transparency and reproducibility.


# Setup 
source("scripts/preprocessing/00_setup.R")  

# Paths 
input_path  <- here("data", "processed", "feature_table_neg.csv")
output_path <- here("data", "processed", "feature_table_neg_clean.csv")
qc_dir      <- here("results", "qc")
logs_dir    <- here("results", "logs")

stopifnot(file.exists(input_path))
dir.create(qc_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

#  Load data 
# Expect: rows = features, columns = samples, first column = feature IDs
feature_table_raw <- read.csv(input_path, row.names = 1, check.names = FALSE)
ft <- feature_table_raw

cat("Initial features:", nrow(ft), "\n")
cat("Initial samples :", ncol(ft), "\n")

# Coerce all columns to numeric once (quietly). Non‑numeric → NA.
suppressWarnings({ ft[] <- lapply(ft, as.numeric) })

# Remove technical samples 
# Drop QC / blanks / pooled mixes by pattern
# They are not biological samples and bias results.
tech_pattern <- "(?i)MM8|QC|ACN"
technical <- grep(tech_pattern, colnames(ft), value = TRUE)

if (length(technical)) {
  tibble::tibble(removed_sample = technical) %>% 
    readr::write_csv(file.path(qc_dir, "technical_samples_removed.csv"))
  ft <- ft %>% dplyr::select(-all_of(technical))
}
cat("After removing technical samples:", ncol(ft), "samples\n")

# NA audit (before)
na_feature_before <- rowSums(is.na(ft))
na_sample_before  <- colSums(is.na(ft))

png(file.path(qc_dir, "NA_per_feature_before.png"), width = 900, height = 600)
hist(na_feature_before, main = "NAs per feature (before)", xlab = "count NAs")
dev.off()

png(file.path(qc_dir, "NA_per_sample_before.png"), width = 900, height = 600)
hist(na_sample_before, main = "NAs per sample (before)", xlab = "count NAs")
dev.off()

# Filter features by missingness 
# Keep features with <=30% missing values
# Too many NAs make data unstable for stats
na_threshold <- 0.30
keep <- rowMeans(is.na(ft)) <= na_threshold
ft <- ft[keep, , drop = FALSE]
cat("After missing filter:", nrow(ft), "features\n")

# Impute remaining NAs with half‑min
# Replace NA with half of row minimum
# NA ≈ below detection limit → use small non‑zero value
replaced_total <- 0L
for (i in seq_len(nrow(ft))) {
  vals <- ft[i, ]
  m <- suppressWarnings(min(as.numeric(vals), na.rm = TRUE))  # row minimum
  if (is.finite(m)) {
    na_idx <- is.na(vals)
    if (any(na_idx)) {
      ft[i, na_idx] <- m / 2
      replaced_total <- replaced_total + sum(na_idx)
    }
  }
}
cat("NAs replaced (half‑min):", replaced_total, "\n")

png(file.path(qc_dir, "NA_per_feature_after.png"), width = 900, height = 600)
hist(rowSums(is.na(ft)), main = "NAs per feature (after)", xlab = "count NAs")
dev.off()

# --- Save cleaned table ----------------------------------------
readr::write_csv(ft %>% tibble::rownames_to_column(var = "feature_id"), output_path)

# Summary tables 
summary_tbl <- tibble::tibble(
  metric = c("n_features", "n_samples", "NAs_total"),
  before = c(nrow(feature_table_raw), ncol(feature_table_raw), sum(is.na(feature_table_raw))),
  after  = c(nrow(ft),               ncol(ft),               sum(is.na(ft)))
)
readr::write_csv(summary_tbl, file.path(qc_dir, "cleaning_summary.csv"))

# Additional change summary (min/mean before vs after)
min_before  <- suppressWarnings(min(as.matrix(feature_table_raw), na.rm = TRUE))
min_after   <- suppressWarnings(min(as.matrix(ft),               na.rm = TRUE))
mean_before <- suppressWarnings(mean(as.matrix(feature_table_raw), na.rm = TRUE))
mean_after  <- suppressWarnings(mean(as.matrix(ft),               na.rm = TRUE))

change_tbl <- tibble::tibble(
  metric = c("min_before", "min_after", "mean_before", "mean_after"),
  value  = c(min_before,    min_after,   mean_before,    mean_after)
)
readr::write_csv(change_tbl, file.path(qc_dir, "imputation_summary.csv"))

cat("\n=== SUMMARY ===\n"); print(summary_tbl)
message("Cleaning finished → ", output_path)

#  Session info 
writeLines(capture.output(sessionInfo()),
           file.path(logs_dir, paste0("sessionInfo_cleaning_", Sys.Date(), ".txt")))











