source("scripts/preprocessing/00_setup.R")

#Read assay (negative mode) 
a_neg <- readr::read_tsv(
  here::here("data","metadata","a_MTBLS4823_LC-MS_negative_reverse-phase_metabolite_profiling.txt"),
  show_col_types = FALSE
)

# Ensure target folders exist
dir.create(here::here("data","metadata"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("results","qc"),    recursive = TRUE, showWarnings = FALSE)

#  Build sample sheet 
# keep columns we actually use; derive file_name and a simple type
sample_sheet_neg <- a_neg %>%
  dplyr::mutate(
    file_name = basename(`Raw Spectral Data File`),
    sample_id = `Sample Name`,
    mode      = `Parameter Value[Scan polarity]`,
    type      = dplyr::case_when(
      stringr::str_detect(`Raw Spectral Data File`, stringr::regex("Blank", ignore_case = TRUE)) ~ "Blank",
      TRUE ~ "Sample"
    )
  ) %>%
  dplyr::select(file_name, sample_id, type, mode)

# Optional: mark QC by filename if present in dataset naming
sample_sheet_neg <- sample_sheet_neg %>%
  dplyr::mutate(
    type = dplyr::case_when(
      stringr::str_detect(file_name, stringr::regex("qc|quality", ignore_case = TRUE)) ~ "QC",
      TRUE ~ type
    )
  )

# Save sample sheet 
readr::write_csv(sample_sheet_neg, here::here("data","metadata","sample_sheet_neg.csv"))
message("sample_sheet_neg.csv has been created in data/metadata/")

#  Consistency checks with processed data (only if available) 
xdata_path <- here::here("data","intermediate","xdata_after_fill.rds")
if (file.exists(xdata_path)) {
  xdata <- readRDS(xdata_path)
  
  # normalize names for matching (lowercase, no extension)
  disk_names <- tolower(tools::file_path_sans_ext(basename(fileNames(xdata))))
  meta_names <- tolower(tools::file_path_sans_ext(sample_sheet_neg$file_name))
  
  #  all processed files must appear in metadata
  all_ok <- all(disk_names %in% meta_names)
  message("Name match (processed ↔ metadata): ", all_ok)
  
  #  quick counts
  message("N files in xdata: ", length(disk_names))
  message("N rows in sample_sheet_neg: ", nrow(sample_sheet_neg))
  
  #  list any files missing in metadata
  mism <- base::setdiff(disk_names, meta_names)
  if (length(mism)) {
    readr::write_csv(
      data.frame(missing_in_metadata = mism),
      here::here("results","qc","group_mismatch.csv")
    )
    warning("Some mzML files are missing in sample_sheet_neg.csv. See results/qc/group_mismatch.csv")
  }
} else {
  message(" Skipping xdata checks: data/intermediate/xdata_after_fill.rds not found (run preprocessing first).")
}