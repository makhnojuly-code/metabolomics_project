###############################################################
# 03_integrate_mummichog_neg.R
#
# Goal:
#   Join Mummichog results (KEGG matches) to feature statistics
#   using rounded m/z and RT values (NEG mode).
#
# Inputs:
#   results/annotation/annotation_candidates_neg.csv
#   results/annotation/mummichog_matched_compound_all.csv
#
# Outputs:
#   results/annotation/features_with_mummichog_neg.csv
#   results/annotation/mummichog_hits_neg.csv
###############################################################

library(dplyr)
library(readr)

# --- Create output folder (annotation) --------------------------

dir.create("results/annotation", recursive = TRUE, showWarnings = FALSE)

# --- 1. Define file paths ---------------------------------------

anno_path <- "results/annotation/annotation_candidates_neg.csv"
mumm_path <- "results/annotation/mummicho_matched_compound_all.csv"

out_full <- "results/annotation/features_with_mummichog_neg.csv"
out_hits <- "results/annotation/mummicho_hits_neg.csv"

# --------------------------------------------------------------
# 2. Load data
# --------------------------------------------------------------

anno <- read_csv(anno_path, show_col_types = FALSE)
mumm <- read_csv(mumm_path, show_col_types = FALSE)

message("Loaded annotation candidates: ", nrow(anno))
message("Loaded Mummichog hits: ", nrow(mumm))

# Check required columns
stopifnot("mzmed" %in% names(anno))
stopifnot("rtmed" %in% names(anno))
stopifnot("Query.Mass" %in% names(mumm))
stopifnot("Retention.Time" %in% names(mumm))

# --------------------------------------------------------------
# 3. Prepare matching keys (rounded m/z and RT)
# --------------------------------------------------------------

anno_key <- anno %>%
  mutate(
    mz_round = round(mzmed, 4),
    rt_round = round(rtmed, 1)
  )

mumm_key <- mumm %>%
  mutate(
    mz_round = round(Query.Mass, 4),
    rt_round = round(Retention.Time, 1)
  )

# --------------------------------------------------------------
# 4. Join annotation candidates with Mummichog results
# --------------------------------------------------------------

feat_mumm <- anno_key %>%
  left_join(mumm_key, by = c("mz_round", "rt_round"))

message("Integrated table rows: ", nrow(feat_mumm))

write_csv(feat_mumm, out_full)
message("Saved full integrated table: ", out_full)

# --------------------------------------------------------------
# 5. Extract only rows with KEGG matches
# --------------------------------------------------------------

hits <- feat_mumm %>%
  filter(!is.na(Matched.Compound)) %>%
  select(
    feature_id,
    mzmed,
    rtmed,
    logFC,
    P.Value,
    adj.P.Val,
    vip_mean,
    Query.Mass,
    Retention.Time,
    Matched.Compound,
    Matched.Form,
    Mass.Diff,
    Empirical.Compound
  )

message("Features with KEGG matches: ", nrow(hits))

write_csv(hits, out_hits)
message("Saved Mummichog hit table: ", out_hits)

message("03_integrate_mummichog_neg.R completed successfully.")