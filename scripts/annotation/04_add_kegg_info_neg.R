###############################################################
# 04_add_kegg_info_neg.R
#
# Goal:
#   1) Get KEGG compound info (name, formula, exact mass)
#      for all matched compounds from Mummichog.
#   2) Join KEGG info back to the full feature table (NEG).
#
# Inputs:
#   results/annotation/mummichog_hits_neg.csv
#   results/annotation/features_with_mummichog_neg.csv
#
# Outputs:
#   results/annotation/kegg_compounds_neg.csv
#   results/annotation/features_with_mummichog_kegg_neg.csv
###############################################################

library(dplyr)
library(readr)
library(KEGGREST)
library(purrr)

# --- Create output folder (annotation) --------------------------

dir.create("results/annotation", recursive = TRUE, showWarnings = FALSE)

# --- 1. Define file paths ---------------------------------------

hits_path <- "results/annotation/mummichog_hits_neg.csv"
full_path <- "results/annotation/features_with_mummichog_neg.csv"

out_kegg <- "results/annotation/kegg_compounds_neg.csv"
out_full <- "results/annotation/features_with_mummichog_kegg_neg.csv"

# --------------------------------------------------------------
# 2. Load hits and collect KEGG IDs
# --------------------------------------------------------------

hits <- read_csv(hits_path, show_col_types = FALSE)
message("Loaded Mummichog hits: ", nrow(hits))

kegg_ids <- hits %>%
  filter(!is.na(Matched.Compound)) %>%
  pull(Matched.Compound) %>%
  unique()

message("Unique KEGG IDs: ", length(kegg_ids))

# --------------------------------------------------------------
# 3. Function to download KEGG compound information
# --------------------------------------------------------------

get_kegg_compound <- function(cid) {
  entry <- tryCatch(
    KEGGREST::keggGet(cid)[[1]],
    error = function(e) NULL
  )
  
  if (is.null(entry)) {
    return(tibble(
      Matched.Compound = cid,
      kegg_name    = NA_character_,
      kegg_formula = NA_character_,
      kegg_mass    = NA_real_
    ))
  }
  
  tibble(
    Matched.Compound = cid,
    kegg_name    = entry$NAME[1],
    kegg_formula = entry$FORMULA,
    kegg_mass    = as.numeric(entry$EXACT_MASS)
  )
}

# --------------------------------------------------------------
# 4. Download KEGG info for all IDs
# --------------------------------------------------------------

kegg_table <- map_dfr(kegg_ids, get_kegg_compound)

write_csv(kegg_table, out_kegg)
message("Saved KEGG compound table: ", out_kegg)

# --------------------------------------------------------------
# 5. Join KEGG info into full feature table
# --------------------------------------------------------------

feat_mumm <- read_csv(full_path, show_col_types = FALSE)

feat_mumm_kegg <- feat_mumm %>%
  left_join(kegg_table, by = "Matched.Compound")

write_csv(feat_mumm_kegg, out_full)
message("Saved final KEGG-integrated feature table: ", out_full)

message("04_add_kegg_info_neg.R completed successfully.")