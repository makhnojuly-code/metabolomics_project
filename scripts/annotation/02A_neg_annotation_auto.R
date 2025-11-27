###############################################################
# 02A_neg_annotation_auto.R — Automated NEG annotation (Mummichog)
#
# Purpose:
#   This script replaces the manual “Mummichog” step normally
#   performed in the MetaboAnalyst web interface.
#
#   It runs Mummichog locally via the MetaboAnalystR package
#   using exactly the same parameters and workflow that the
#   online platform uses.
#
# Input:
#   results/annotation/metabo_mummichog_NEG.txt
#      - the peak list prepared in Script 02
#      - same file previously uploaded to the MetaboAnalyst web interface
#
# Output (auto-generated):
#   results/annotation/mummicho_matched_compound_all.csv
#   results/annotation/mummicho_pathway_enrichment.csv
#
# Notes:
#   • This script is specific for NEG-mode analysis.
#   • The KEGG library used here is Arabidopsis (“ath_kegg”).
#   • Parameters (ppm, p-value cutoff, enrichment mode) match
#     those used in the MetaboAnalyst GUI export.
###############################################################

## ------------------------------------------------------------
## 0. Load libraries
## ------------------------------------------------------------
library(MetaboAnalystR)   # Mummichog / pathway analysis
library(readr)            # fast file reading
library(dplyr)            # data handling
library(RJSONIO)          # required internally by MetaboAnalystR
library(fitdistrplus)     # required for background distribution fitting

## ------------------------------------------------------------
## 1. Define paths
## ------------------------------------------------------------
input_file <- "results/annotation/metabo_mummichog_NEG.txt"

# Output directory (same as manual annotation for compatibility)
out_dir <- "results/annotation"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## ------------------------------------------------------------
## 2. Validate input file
## ------------------------------------------------------------
if (!file.exists(input_file)) {
  stop(
    "Input file not found: ", input_file, "\n",
    "Make sure Script 02 created metabo_mummichog_NEG.txt."
  )
}

message("Using input file: ", input_file)
message("Output directory: ", out_dir)

## ------------------------------------------------------------
## 3. Initialize MetaboAnalystR object (NEG mode)
## ------------------------------------------------------------
# MetaboAnalystR (current version) requires 4 arguments:
#   data type     = "mass_all"
#   analysis type = "mummichog"
#   paired data   = FALSE
#   number of permutations (initial) = 150
mSet <- InitDataObjects("mass_all", "mummichog", FALSE, 150)

## ------------------------------------------------------------
## 4. Define peak format
## ------------------------------------------------------------
# "rmp" = row format with: m/z, p-value, t-score, retention time.
# This must match the structure of metabo_mummichog_NEG.txt.
mSet <- SetPeakFormat(mSet, "rmp")

## ------------------------------------------------------------
## 5. Instrument parameters
## ------------------------------------------------------------
# These parameters must match those used in the MetaboAnalyst GUI export:
#   5.0 ppm         — mass accuracy
#   "negative"      — ionization mode
#   "no"            — no MS2 information
#   0.02            — p-value threshold used for ranking peaks
mSet <- UpdateInstrumentParameters(
  mSet,
  5.0,
  "negative",
  "no",
  0.02
)

## ------------------------------------------------------------
## 6. Load peak list
## ------------------------------------------------------------
mSet <- Read.PeakListData(mSet, input_file)

# Optional: sanity check on peak list structure and valid values
mSet <- SanityCheckMummichogData(mSet)

## ------------------------------------------------------------
## 7. Configure enrichment analysis
## ------------------------------------------------------------
# Integrated method = combines Mummichog + GSEA in a single workflow
mSet <- SetPeakEnrichMethod(mSet, "integ", "v2")

# Strict significance threshold used in original MetaboAnalyst script
mSet <- SetMummichogPval(mSet, 1e-5)

## ------------------------------------------------------------
## 8. Run Mummichog analysis (Arabidopsis KEGG)
## ------------------------------------------------------------
# Arguments:
#   "ath_kegg"  — KEGG library for Arabidopsis thaliana
#   "current"   — latest version of pathway database
#   3           — minimum number of hits per pathway
#   100         — number of permutations (can be increased)
mSet <- PerformPSEA(mSet, "ath_kegg", "current", 3, 100)

## ------------------------------------------------------------
## 9. Find Mummichog output files
## ------------------------------------------------------------
# MetaboAnalystR writes CSV files directly into the working directory.
all_mum <- list.files(pattern = "^mummi.*\\.csv$", full.names = TRUE)

message("Detected Mummichog-related CSV files:")
print(all_mum)

# Identify matched-metabolite file
comp_candidates <- all_mum[grepl("matched_compound_all", all_mum)]

# Identify pathway-enrichment file
pth_candidates  <- all_mum[grepl("pathway_enrichment", all_mum)]

# Safety check
if (length(comp_candidates) == 0L || length(pth_candidates) == 0L) {
  stop(
    "Mummichog output files were not generated.\n",
    "Files detected:\n",
    paste(all_mum, collapse = "\n")
  )
}

compounds_file <- comp_candidates[1]
pathways_file  <- pth_candidates[1]

## ------------------------------------------------------------
## 10. Save copies to results/annotation (standardized names)
## ------------------------------------------------------------
out_compounds <- file.path(out_dir, "mummicho_matched_compound_all.csv")
out_pathways  <- file.path(out_dir, "mummicho_pathway_enrichment.csv")

file.copy(compounds_file, out_compounds, overwrite = TRUE)
file.copy(pathways_file,  out_pathways,  overwrite = TRUE)

message("✔ Automatic NEG annotation completed successfully.")
message("Matched compounds saved to:  ", out_compounds)
message("Pathway enrichment saved to: ", out_pathways)

