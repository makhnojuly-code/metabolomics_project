# scripts/preprocessing/01_preprocessing.R
source("scripts/preprocessing/00_setup.R")

# # switch for RT alignment: FALSE = serial mode (more stable)
use_parallel <- FALSE  

# Path to raw data 
raw_dir <- here::here("data", "raw", "neg_mode")

#  Find mzML files 
raw_files <- list.files(
  raw_dir,
  pattern = "\\.mzML$",
  full.names = TRUE,
  ignore.case = TRUE
)
stopifnot(length(raw_files) > 0)

# Decide groups: prefer metadata; robust match; fallback to filename patterns
meta_path <- here::here("data","metadata","sample_sheet_neg.csv")

files       <- basename(raw_files)
files_clean <- tolower(tools::file_path_sans_ext(files))

if (file.exists(meta_path)) {
  sheet <- readr::read_csv(meta_path, show_col_types = FALSE)
  
  meta_names <- basename(sheet$file_name)
  meta_clean <- tolower(tools::file_path_sans_ext(meta_names))
  
  idx <- match(files_clean, meta_clean)
  
  # take type from metadata ("Blank"/"Sample"), use lower-case
  sample_groups <- tolower(sheet$type[idx])  # may contain NA if not matched
  
  # mark QC by filename only where metadata is NA
  is_qc_name <- grepl("(^|[^a-z])qc([^a-z]|$)|quality", files, ignore.case = TRUE)
  sample_groups[is_qc_name & is.na(sample_groups)] <- "qc"
  
  # fallback by filename ONLY for those not found in metadata
  need_fb <- is.na(sample_groups)
  if (any(need_fb)) {
    is_blank_fb <- grepl("blank|blk|blanc|zero", files[need_fb], ignore.case = TRUE)
    is_qc_fb    <- grepl("(^|[^a-z])qc([^a-z]|$)|quality", files[need_fb], ignore.case = TRUE)
    sample_groups[need_fb] <- ifelse(is_blank_fb, "blank",
                                     ifelse(is_qc_fb, "qc", "sample"))
    warning("Labeled by filename (not in metadata): ",
            paste(files[need_fb], collapse = ", "))
  }
} else {
  # No metadata → detect by filename
  is_blank <- grepl("blank|blk|blanc|zero", files, ignore.case = TRUE)
  is_qc    <- grepl("(^|[^a-z])qc([^a-z]|$)|quality", files, ignore.case = TRUE)
  sample_groups <- ifelse(is_blank, "blank", ifelse(is_qc, "qc", "sample"))
}

stopifnot(length(sample_groups) == length(files))



# Import raw (onDisk saves memory)
raw_data <- readMSData(files = raw_files, mode = "onDisk")
saveRDS(raw_data, here::here("data", "intermediate", "raw_data_imported.rds"))


# Peak picking (detect chromatographic peaks)
cwp <- CentWaveParam(
  ppm = 10,
  peakwidth = c(5, 40),
  snthresh = 10,
  prefilter = c(3, 1000),
  noise = 500,
  mzdiff = -0.001,
  mzCenterFun = "wMean"
)
xdata <- findChromPeaks(raw_data, param = cwp)
saveRDS(xdata, here::here("data", "intermediate", "xdata_peaks.rds"))

# Retention time alignment
# Load checkpoint  (after peak picking)
xdata <- readRDS(here::here("data", "intermediate","xdata_peaks.rds"))

# Pick RT reference automatically: prefer a QC file; else the file with TIC closest to the median (most typical)
tic_by_file <- vapply(split(tic(xdata), fromFile(xdata)), sum, numeric(1))
is_qc_name  <- grepl("qc|quality", basename(fileNames(xdata)), ignore.case = TRUE)
cands       <- if (any(is_qc_name)) which(is_qc_name) else seq_along(tic_by_file)
center_idx  <- cands[ which.min(abs(tic_by_file[cands] - median(tic_by_file[cands]))) ]

# Temporary safety mode:
# Run RT alignment in single-thread mode because Obiwarp can be unstable in parallel.
# After alignment finishes, the original parallel backend is restored automatically.
if (!use_parallel) {
  old_bp <- BiocParallel::bpparam()
  BiocParallel::register(BiocParallel::SerialParam())
  on.exit(BiocParallel::register(old_bp), add = TRUE)
}

# RT alignment might fail on some datasets.
# tryCatch prevents the pipeline from stopping — if Obiwarp fails,
# we keep the original (unaligned) xdata and continue.

xdata <- tryCatch(
  {
    adjustRtime(xdata, param = ObiwarpParam(centerSample = center_idx))
  },
  error = function(e) {
    message("!!! RT alignment failed. Returning UNALIGNED xdata. Error: ", e$message)
    xdata
  }
)
saveRDS(xdata, here::here("data", "intermediate","xdata_rtime.rds"))
gc()



# saves a small text report, so later it is easy to check
# if alignment worked and which sample was the reference
# without re-running analysis
rt_summary <- summary(unlist(rtime(xdata))) #normally prints to the console.
capture.output(rt_summary, file = here::here("results", "qc", "qc_rt_summary.txt"))#takes this printed text and saves it into a .txt file — so we can keep this information on disk, not only on the screen.
writeLines(
  c(
    paste0("hasAdjustedRtime: ", hasAdjustedRtime(xdata)),
    paste0("RT center sample index: ", center_idx, " (", basename(fileNames(xdata))[center_idx], ")")
  ),
  here::here("results", "qc", "qc_rt_alignment_report.txt")
)

# Feature grouping (align peaks across all samples)
stopifnot(length(sample_groups) == length(fileNames(xdata)))#checks that the number of files equals the number of groups (sample / blank / qc).
#If they don’t match — stop, this means an error

xdata <- groupChromPeaks(
  xdata,
  param = PeakDensityParam(sampleGroups = factor(sample_groups))
) ## group peaks across samples into features, using sample_groups (sample / blank / qc)
saveRDS(xdata, here::here("data", "intermediate","xdata_grouped.rds"))
cat("Features after grouping:", nrow(featureDefinitions(xdata)), "\n")

# QC step checks how many peaks were found in each file. If one file has much fewer (or much more) peaks compared to the others, it can be a sign of a problem (bad injection, low signal, instrument issue). Saving this simple table helps to quickly identify suspicious samples before statistical analysis.

pk <- chromPeaks(xdata)
write.csv(as.data.frame(table(pk[,"sample"])),
          here::here("results", "qc", "qc_peaks_per_file.csv"),
          row.names = FALSE)

# Fill missing peak intensities (reduce NAs and make feature matrix more complete)
xdata <- fillChromPeaks(xdata)
saveRDS(xdata, here::here("data", "intermediate","xdata_after_fill.rds"))
xdata <- readRDS(here::here("data", "intermediate","xdata_after_fill.rds"))

# Quick NA snapshot before cleaning
ft_raw <- featureValues(xdata, value = "into", method = "medret")
na_feat_ratio <- rowMeans(is.na(ft_raw))
write.csv(
  data.frame(features = nrow(ft_raw),
             samples  = ncol(ft_raw),
             NA_total = sum(is.na(ft_raw)),
             NA_feat_median = median(na_feat_ratio, na.rm = TRUE)),
  here::here("results", "qc", "qc_fill_summary.csv"),
  row.names = FALSE
)

#  Export tables
stopifnot(hasChromPeaks(xdata), hasFeatures(xdata))

feature_table <- featureValues(xdata, value="into", method="medret")
write.csv(feature_table, here::here("data","processed","feature_table_neg.csv"), row.names = TRUE)

fdef <- as.data.frame(featureDefinitions(xdata))
fdef_simple <- fdef[, sapply(fdef, function(x) !is.list(x))] %>%
  tibble::rownames_to_column("feature_id")
write.csv(fdef_simple, here::here("data","processed","feature_definitions.csv"), row.names = FALSE)

ft_df <- as.data.frame(feature_table) %>% tibble::rownames_to_column("feature_id")
tidy_tbl <- ft_df %>%
  pivot_longer(cols = -feature_id, names_to="sample", values_to="intensity") %>%
  dplyr::left_join(fdef_simple, by="feature_id")
write.csv(tidy_tbl, here::here("data","processed","tidy_feature_table.csv"), row.names = FALSE)

# Record environment
writeLines(capture.output(sessionInfo()),
           here::here("results", "logs", paste0("sessionInfo_", Sys.Date(), ".txt"))
)

# Quick check
ft <- read.csv(here::here("data","processed","feature_table_neg.csv"), row.names = 1)
dim(ft); summary(ft[,1:5])


