# ================================================================
# 01_build_annotation_candidates_neg.R
#
# Goal:
#   Build a table of important NEG-mode features for annotation
#   and create an input file for MetaboAnalyst (mummichog).
# ================================================================

# --- 0. Load required packages ----------------------------------

library(dplyr)  # data manipulation
library(readr)  # read/write CSV

message("Packages loaded: dplyr, readr")


# --- Create output folder (annotation) ---------------------------

dir.create("results/annotation", recursive = TRUE, showWarnings = FALSE)

# --- 1. Define file paths ---------------------------------------

vip_path     <- "results/plsda_vip_neg_highVIP.csv"
limma_path   <- "results/limma_results_neg.csv"
featdef_path <- "data/processed/feature_definitions.csv"

out_anno_path   <- "results/annotation/annotation_candidates_neg.csv"
out_metabo_path <- "results/annotation/metaboanalyst_input_neg.csv"


# --- 2. Read input tables ---------------------------------------

vip_df <- read_csv(vip_path, show_col_types = FALSE)
message("High-VIP table loaded from: ", vip_path,
        " | n = ", nrow(vip_df), " features.")

limma_df_raw <- read_csv(limma_path, show_col_types = FALSE)
message("Limma results loaded from: ", limma_path,
        " | n = ", nrow(limma_df_raw), " rows.")

featdef_df <- read_csv(featdef_path, show_col_types = FALSE)
message("Feature definitions loaded from: ", featdef_path,
        " | n = ", nrow(featdef_df), " features.")


# --- 3. Harmonise column names ----------------------------------
# Make sure limma table has a column called "feature_id".

limma_df <- limma_df_raw

if (!"feature_id" %in% names(limma_df)) {
  names(limma_df)[1] <- "feature_id"  # assume first column is feature ID
}

if (!"feature_id" %in% names(limma_df)) {
  stop("Could not find or create a 'feature_id' column in limma results.")
}

# Check that all tables have "feature_id"
stopifnot("feature_id" %in% names(vip_df))
stopifnot("feature_id" %in% names(featdef_df))


# --- 4. Merge VIP + limma + feature definitions -----------------
# Start from VIP features and add limma + feature metadata.

anno_candidates <- vip_df %>%
  left_join(limma_df,  by = "feature_id") %>%
  left_join(featdef_df, by = "feature_id")

message("Merged annotation table created. Rows: ", nrow(anno_candidates))


# --- 5. Add helper columns and sort -----------------------------
# Add -log10(FDR) if adj.P.Val is available.

if ("adj.P.Val" %in% names(anno_candidates)) {
  anno_candidates <- anno_candidates %>%
    mutate(
      neg_log10_adjP = ifelse(
        is.na(adj.P.Val),
        NA_real_,
        -log10(adj.P.Val)
      )
    )
}

# Sort by VIP (high to low), then by FDR (low to high) if present.
if ("adj.P.Val" %in% names(anno_candidates)) {
  anno_candidates <- anno_candidates %>%
    arrange(desc(vip_mean), adj.P.Val)
} else {
  anno_candidates <- anno_candidates %>%
    arrange(desc(vip_mean))
}


# --- 6. Save full annotation candidate table --------------------

write_csv(anno_candidates, out_anno_path)

message("Annotation candidate table saved to: ", out_anno_path)
message("Columns in annotation table: ", paste(names(anno_candidates), collapse = ", "))


# --- 7. Create simplified MetaboAnalyst input -------------------
# Use:
#   mz      = mzmed
#   p.value = P.Value (if available)
#   Extra columns (rt, logFC, adj.P.Val, vip_mean) are for our reference.

if (!"mzmed" %in% names(anno_candidates)) {
  stop("Column 'mzmed' not found in annotation table. Check feature_definitions.")
}

metabo_input <- anno_candidates %>%
  filter(!is.na(mzmed)) %>%  # need valid mz
  transmute(
    feature_id = feature_id,
    mz         = mzmed,
    p.value    = if ("P.Value" %in% names(anno_candidates)) P.Value else NA_real_,
    rt         = if ("rtmed" %in% names(anno_candidates)) rtmed else NA_real_,
    logFC      = if ("logFC" %in% names(anno_candidates)) logFC else NA_real_,
    adj.P.Val  = if ("adj.P.Val" %in% names(anno_candidates)) adj.P.Val else NA_real_,
    vip_mean   = vip_mean
  )

write_csv(metabo_input, out_metabo_path)

message("MetaboAnalyst input table saved to: ", out_metabo_path)
message("Rows in MetaboAnalyst input: ", nrow(metabo_input))

message("01_build_annotation_candidates_neg.R completed successfully.")