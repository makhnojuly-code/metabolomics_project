###############################################################
# 02_export_mummichog_input_neg.R
#
# Goal:
#   Create the input file for MetaboAnalyst (Mummichog).
#   Format required: m.z | p.value | t.score | r.t
#
# Output:
#   results/annotation/metabo_mummichog_NEG.txt
###############################################################

library(dplyr)   # data manipulation
library(readr)   # read CSV

# --- Create output folder (annotation) --------------------------

dir.create("results/annotation", recursive = TRUE, showWarnings = FALSE)

# --- 1. Define file paths ---------------------------------------

input_path  <- "results/annotation/metaboanalyst_input_neg.csv"
output_path <- "results/annotation/metabo_mummichog_NEG.txt"

# --- 2. Read input table ----------------------------------------

df <- read_csv(input_path, show_col_types = FALSE)
message("Loaded MetaboAnalyst input: ", nrow(df), " rows.")

# --- 3. Build Mummichog table ----------------------------------
# Required columns:
#   m.z     = mz
#   p.value = p.value
#   t.score = logFC
#   r.t     = rt

mummichog_in <- df %>%
  transmute(
    `m.z`   = mz,
    p.value = p.value,
    t.score = logFC,
    `r.t`   = rt
  )

# --- 4. Save as TAB-delimited text file -------------------------

write.table(
  mummichog_in,
  file = output_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("Mummichog input saved to: ", output_path)
message("Upload this file to MetaboAnalyst → MS Peaks to Pathways.")