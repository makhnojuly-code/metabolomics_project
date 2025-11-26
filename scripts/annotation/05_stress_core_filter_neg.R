###############################################################
# 05_stress_core_filter_neg.R
#
# Goal:
#   From NEG-mode feature table with Mummichog + KEGG:
#     1) Select "core features" (strong stress signals).
#     2) Aggregate them to "core metabolites" (KEGG level).
#     3) Add simple stress-related groups.
#
# Input:
#   results/annotation/features_with_mummichog_kegg_neg.csv
#
# Outputs:
#   results/annotation/core_features_neg.csv
#   results/annotation/core_metabolites_neg.csv
###############################################################

#########################
# 0. Load libraries
#########################

library(dplyr)
library(readr)
library(stringr)

message("Packages loaded: dplyr, readr, stringr")

# Create output folder if needed
dir.create("results/annotation", recursive = TRUE, showWarnings = FALSE)

#########################
# 1. Define thresholds
#########################

pval_cutoff  <- 0.05  # adj.P.Val <= this
logfc_cutoff <- 1.0   # |logFC|   >= this
vip_cutoff   <- 1.0   # VIP mean  >= this

# Minimal number of features per KEGG metabolite
min_features_per_metabolite <- 1

message("Using thresholds:")
message("  adj.P.Val  <= ", pval_cutoff)
message("  |logFC|    >= ", logfc_cutoff)
message("  VIP_mean   >= ", vip_cutoff)
message("  min features per metabolite = ", min_features_per_metabolite)

#########################
# 2. Load input table
#########################

in_path <- "results/annotation/features_with_mummichog_kegg_neg.csv"

feat_df <- read_csv(in_path, show_col_types = FALSE)

message("Feature table loaded from: ", in_path)
message("Rows: ", nrow(feat_df), " | Columns: ", ncol(feat_df))

# Required columns for this script
required_cols <- c(
  "feature_id",
  "logFC",
  "adj.P.Val",
  "vip_mean",
  "Matched.Compound",  # KEGG ID
  "kegg_name",
  "kegg_formula"
)

if (!all(required_cols %in% names(feat_df))) {
  stop("Some required columns are missing in features_with_mummichog_kegg_neg.csv. 
Please check the names: ", paste(required_cols, collapse = ", "))
}

#########################
# 3. Select Core Features
#########################
# Keep only:
#   - KEGG-matched features (Matched.Compound not NA)
#   - non-NA stats
#   - passing thresholds (adj.P.Val, logFC, VIP)

core_features <- feat_df %>%
  filter(!is.na(Matched.Compound)) %>%
  filter(!is.na(adj.P.Val),
         !is.na(logFC),
         !is.na(vip_mean)) %>%
  filter(
    adj.P.Val <= pval_cutoff,
    abs(logFC) >= logfc_cutoff,
    vip_mean  >= vip_cutoff
  ) %>%
  arrange(adj.P.Val)

message("Core features after filtering: ", nrow(core_features))

message("Example rows from core_features:")
print(
  core_features %>%
    select(
      feature_id, mzmed, rtmed,
      logFC, adj.P.Val, vip_mean,
      Matched.Compound, kegg_name
    ) %>%
    head(10)
)

#########################
# 4. Aggregate Core Features to Metabolites
#########################
# One KEGG metabolite can have multiple LC-MS features.
# Aggregate by Matched.Compound (KEGG ID).

core_metabolites <- core_features %>%
  group_by(Matched.Compound, kegg_name, kegg_formula) %>%
  summarise(
    n_features     = n(),
    min_adjP       = min(adj.P.Val, na.rm = TRUE),
    max_abs_logFC  = max(abs(logFC), na.rm = TRUE),
    mean_logFC     = mean(logFC, na.rm = TRUE),
    max_VIP        = max(vip_mean, na.rm = TRUE),
    .groups        = "drop"
  ) %>%
  filter(n_features >= min_features_per_metabolite) %>%
  arrange(min_adjP)

# Simple direction label based on mean_logFC
core_metabolites <- core_metabolites %>%
  mutate(
    direction = case_when(
      mean_logFC >  0.5  ~ "Up in stress",
      mean_logFC < -0.5  ~ "Down in stress",
      TRUE               ~ "Weak / mixed"
    )
  )

message("Core metabolites after aggregation: ", nrow(core_metabolites))

message("Example rows from core_metabolites:")
print(
  core_metabolites %>%
    select(
      Matched.Compound, kegg_name,
      n_features, min_adjP,
      mean_logFC, max_abs_logFC, max_VIP,
      direction
    ) %>%
    head(15)
)

#########################
# 5. Create Short Names for Metabolites
#########################
# Make shorter names for plots: first synonym, no brackets, max 25 chars.

core_metabolites <- core_metabolites %>%
  mutate(
    short_name = kegg_name %>%
      str_replace(";.*", "") %>%              # keep first synonym
      str_replace("\\s*\\(.*\\)", "") %>%     # remove (...)
      str_trunc(25, side = "right", ellipsis = "…")
  )

message("Example of short_name:")
print(
  core_metabolites %>%
    select(kegg_name, short_name) %>%
    head(10)
)

#########################
# 6. First simple "stress core" grouping
#########################
# Rough manual classification based on name patterns.

assign_group <- function(name) {
  name_low <- tolower(name)
  
  # JA / oxylipin
  if (str_detect(name_low, "hpot|opda|jasmon|eot|oxo|hydroperoxy|octadeca"))
    return("JA / oxylipin")
  
  # Aromatic AA / shikimate
  if (str_detect(name_low, "anthranil|aminobenzo|hydroxybenzo|shik"))
    return("Aromatic AA / shikimate")
  
  # Flavonoids / phenylpropanoids
  if (str_detect(name_low, "quercetin|kaempferol|flav|anthocyan|caffe|coumar|ferul|sinap"))
    return("Flavonoid / phenylprop")
  
  # Fatty acids / lipid-related
  if (str_detect(name_low, "linole|oleic|palmit|stear|hexadec"))
    return("Fatty acids / lipids")
  
  # Sugar phosphates
  if (str_detect(name_low, "glucose|fructose|mannose|galactose|phosphate"))
    return("Sugar phosphates")
  
  # Nucleotides
  if (str_detect(name_low, "ump|uridine|pseudouridine|nucleotide|purine|pyrimidine"))
    return("Nucleotides")
  
  # Plant hormones
  if (str_detect(name_low, "gibberell|auxin|cytokin|abscisic|salicy"))
    return("Hormone / signaling")
  
  # Terpenoids
  if (str_detect(name_low, "caryophyll|bisabol|humulen|barbaten|sesquiterp|monoterp|diterp"))
    return("Terpenoids")
  
  return("Other / unknown")
}

core_metabolites <- core_metabolites %>%
  rowwise() %>%
  mutate(stress_core_group = assign_group(kegg_name)) %>%
  ungroup()

#########################
# 7. Save output tables
#########################

out_core_feat_path <- "results/annotation/core_features_neg.csv"
out_core_meta_path <- "results/annotation/core_metabolites_neg.csv"

write_csv(core_features,    out_core_feat_path)
write_csv(core_metabolites, out_core_meta_path)

message("Core features saved to:    ", out_core_feat_path)
message("Core metabolites saved to: ", out_core_meta_path)

#########################
# 8. Summary messages
#########################

message("Summary:")
message("  Total features in original table: ", nrow(feat_df))
message("  Core features (strong + KEGG):    ", nrow(core_features))
message("  Core metabolites (KEGG-level):    ", nrow(core_metabolites))

message("05_stress_core_filter_neg.R completed successfully.")