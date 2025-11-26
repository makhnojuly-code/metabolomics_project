###############################################################
# 06_validate_core_neg.R
#
# Goal:
#   Validate core_features_neg and core_metabolites_neg before
#   visualization and biological interpretation (NEG mode).
#
# Output:
#   results/annotation/validation_report_neg.txt
###############################################################

library(dplyr)
library(readr)
library(stringr)

message("Running validation for NEG-mode core tables...")

# --- Create output folder (annotation) --------------------------

dir.create("results/annotation", recursive = TRUE, showWarnings = FALSE)

#########################
# 1. Load tables
#########################

core_feat_path <- "results/annotation/core_features_neg.csv"
core_meta_path <- "results/annotation/core_metabolites_neg.csv"

core_features    <- read_csv(core_feat_path, show_col_types = FALSE)
core_metabolites <- read_csv(core_meta_path, show_col_types = FALSE)

#########################
# 2. Start building report
#########################

report <- c()
report <- c(report, "VALIDATION REPORT (NEG mode)")
report <- c(report, "----------------------------------")
report <- c(report, paste("Core features:", nrow(core_features)))
report <- c(report, paste("Core metabolites:", nrow(core_metabolites)))
report <- c(report, "")

#########################
# 3. Column check
#########################
# Here we check if the key columns exist in core_metabolites.

required_meta_cols <- c(
  "Matched.Compound",
  "kegg_name",
  "n_features",
  "min_adjP",
  "mean_logFC",
  "max_abs_logFC",
  "max_VIP",
  "direction",
  "short_name",
  "stress_core_group"   # from script 05
)

missing_cols <- setdiff(required_meta_cols, names(core_metabolites))

if (length(missing_cols) > 0) {
  report <- c(report, "Missing columns in core_metabolites:")
  report <- c(report, paste(missing_cols, collapse = ", "))
} else {
  report <- c(report, "All required columns present.")
}

report <- c(report, "")

#########################
# 4. Duplicate KEGG ID check
#########################

dup_ids <- core_metabolites$Matched.Compound[
  duplicated(core_metabolites$Matched.Compound)
]

if (length(dup_ids) > 0) {
  report <- c(report, "Duplicate KEGG IDs found:")
  report <- c(report, paste(dup_ids, collapse = ", "))
} else {
  report <- c(report, "No duplicate KEGG IDs.")
}

report <- c(report, "")

#########################
# 5. NA check
#########################

na_cols <- names(core_metabolites)[colSums(is.na(core_metabolites)) > 0]

if (length(na_cols) > 0) {
  report <- c(report, "Columns containing NA values:")
  for (col in na_cols) {
    report <- c(
      report,
      paste(col, "→", sum(is.na(core_metabolites[[col]])), "missing")
    )
  }
} else {
  report <- c(report, "No NA values in core_metabolites.")
}

report <- c(report, "")

#########################
# 6. Join consistency check
#    Every KEGG ID in core_metabolites must exist in core_features
#########################

core_K <- unique(core_features$Matched.Compound)
meta_K <- unique(core_metabolites$Matched.Compound)

missing_in_features <- setdiff(meta_K, core_K)

if (length(missing_in_features) > 0) {
  report <- c(report,
              "Some KEGG IDs in core_metabolites do not exist in core_features:")
  report <- c(report, paste(missing_in_features, collapse = ", "))
} else {
  report <- c(report, "All metabolites come from core features.")
}

report <- c(report, "")

#########################
# 7. Validate direction-label consistency
#########################

incorrect_up <- core_metabolites %>%
  filter(direction == "Up in stress", mean_logFC <= 0.5)

incorrect_down <- core_metabolites %>%
  filter(direction == "Down in stress", mean_logFC >= -0.5)

if (nrow(incorrect_up) == 0 && nrow(incorrect_down) == 0) {
  report <- c(report, "Direction labels consistent with mean_logFC.")
} else {
  report <- c(report, "Direction label inconsistencies found:")
  if (nrow(incorrect_up) > 0) {
    report <- c(
      report,
      paste("Up in stress but mean_logFC <= 0.5:", nrow(incorrect_up))
    )
  }
  if (nrow(incorrect_down) > 0) {
    report <- c(
      report,
      paste("Down in stress but mean_logFC >= -0.5:", nrow(incorrect_down))
    )
  }
}

report <- c(report, "")

#########################
# 8. Pathway assignment check (optional)
#########################
# If later you add a column with pathways, you can rename this
# logic. For now we check that the column exists before using it.

if ("pathway_assignment" %in% names(core_metabolites)) {
  
  all_paths <- core_metabolites$pathway_assignment
  
  if (all(all_paths != "")) {
    report <- c(report,
                "pathway_assignment column exists and is filled.")
  } else {
    report <- c(report,
                "Empty values detected in pathway_assignment.")
  }
  
  num_paths <- sapply(strsplit(all_paths, ";"), length)
  
  report <- c(report,
              paste("Average pathways per metabolite:",
                    round(mean(num_paths), 2)))
  report <- c(report,
              paste("Max pathways per metabolite:", max(num_paths)))
  report <- c(report,
              paste("Min pathways per metabolite:", min(num_paths)))
  report <- c(report, "")
  
} else {
  report <- c(report,
              "Column 'pathway_assignment' not found (skipping pathway checks).")
  report <- c(report, "")
}

#########################
# 9. Save report
#########################

out_path <- "results/annotation/validation_report_neg.txt"
writeLines(report, out_path)

message("Validation completed. Report saved to: ", out_path)