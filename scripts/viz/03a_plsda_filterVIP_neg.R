# ================================================================
# 03a_plsda_highVIP_neg.R
#
# Goal:
#   Extract a focused list of important features from the PLS-DA VIP
#   table for NEG mode (features with VIP_mean > 1.0).
#
# Context:
#   This script is a small follow-up to 03_plsda_neg.R.
#   The PLS-DA script saves a full VIP table for all features:
#     results/plsda_vip_neg.csv
#
#   Here, we select only high-VIP features, which will be used as
#   candidates for downstream biological interpretation, annotation
#   and pathway analysis.
#
# Input:
#   results/plsda_vip_neg.csv
#
# Output:
#   results/plsda_vip_neg_highVIP.csv
#     - subset of features with VIP_mean > 1.0
#
# Notes:
#   - The threshold VIP_mean > 1.0 is a common rule of thumb for
#     identifying important features in PLS-DA.
# ================================================================

# --- 0. Load packages (base R is enough here) -------------------

# No additional packages are strictly required.
# If dplyr is available, it could be used, but here we keep it simple.


# --- 1. Define file paths ---------------------------------------

vip_path   <- "results/plsda_vip_neg.csv"
out_path   <- "results/plsda_vip_neg_highVIP.csv"


# --- 2. Read full VIP table -------------------------------------

vip_df <- read.csv(vip_path, check.names = FALSE)

if (!"vip_mean" %in% names(vip_df)) {
  stop("Column 'vip_mean' not found in VIP table: ", vip_path)
}

message("Full VIP table loaded: ", nrow(vip_df), " features.")


# --- 3. Filter high-VIP features (VIP_mean > 1.0) ---------------

vip_high <- subset(vip_df, vip_mean > 1)

message("Number of high-VIP features (VIP_mean > 1.0): ",
        nrow(vip_high))


# --- 4. Save filtered list --------------------------------------

write.csv(
  vip_high,
  file = out_path,
  row.names = FALSE
)

message("High-VIP feature list saved to: ", out_path)
message("03a_plsda_highVIP_neg.R completed successfully.")