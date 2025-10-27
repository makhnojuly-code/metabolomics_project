source("scripts/preprocessing/00_setup.R")


# 1. Specify the path to the raw data
raw_dir <- here::here("data", "raw", "neg_mode")

# 2. Search for mzML files in the folder

raw_files <- list.files(
  raw_dir,
  pattern = "\\.mzML$",
  full.names = TRUE,
  ignore.case = TRUE
)

# 3. Import raw data (onDisk = saves memory, loads data on demand)
raw_data <- readMSData(files = raw_files, mode = "onDisk")
saveRDS(raw_data, here::here("results", "raw_data_imported.rds"))

# 4. Peak picking (detect chromatographic peaks)
cwp <- CentWaveParam()
xdata <- findChromPeaks(raw_data, param = cwp)
saveRDS(xdata, here::here("results", "xdata_s4_peaks.rds"))   


# 5. Retention time alignment

# Load checkpoint from step 4 (after peak picking)
xdata <- readRDS(here::here("results","xdata_s4_peaks.rds"))

# Temporarily switch to single-thread mode
# (Sometimes parallel Obiwarp alignment can cause errors — this avoids it)
old_bp <- BiocParallel::bpparam()
BiocParallel::register(BiocParallel::SerialParam())

# Perform RT alignment
# 'centerSample = 130' sets sample #130 as the reference for alignment
xdata <- adjustRtime(xdata, param = ObiwarpParam(centerSample = 130))
saveRDS(xdata, here::here("results","xdata_s5_rtime.rds"))


# 5.1 Quick QC after RT alignment
xdata <- readRDS(here::here("results","xdata_s5_rtime.rds"))
# Check if retention time was adjusted
# If this returns TRUE → alignment worked
hasAdjustedRtime(xdata)

# Plot Base Peak Chromatograms (BPC)
# This shows the shape of the signal for each sample
bpis <- chromatogram(xdata, aggregationFun = "max")
cols <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(fileNames(xdata)))

plot(
  bpis,
  col = cols,
  main = "BPC after RT alignment",
  xlab = "Retention time (sec)",
  ylab = "Intensity"
)

# Plot Total Ion Current (TIC)
# This shows total signal per file
tic_values <- tic(xdata)
file_ids   <- fromFile(xdata)

boxplot(
  split(tic_values, file_ids),
  main = "TIC per file (after RT alignment)",
  xlab = "File index",
  ylab = "TIC intensity",
  col = "lightblue"
)

# Check retention time range
# This helps to see where the signal ends
rt_summary <- summary(unlist(rtime(xdata)))
rt_summary

# Restore the previous parallel backend
# (so the next steps can run in parallel again)
BiocParallel::register(old_bp)
# “Base peak chromatograms after RT alignment show consistent signal across all samples. Major peaks align well in retention time, indicating successful RT correction. No flat or anomalous chromatograms were detected. Minor variation at the end of the run (1000–1200 sec) is likely due to natural differences between runs.”

# 6. Feature grouping (align peaks across all samples)
xdata <- groupChromPeaks(xdata, param = PeakDensityParam(sampleGroups = factor(rep("all", length(fileNames(xdata))))))
saveRDS(xdata, here::here("results", "xdata_s6_grouped.rds")) 


# 6.5. Fill missing peaks by integration over RT–m/z windows
xdata <- fillChromPeaks(xdata)                                   
saveRDS(xdata, here::here("results", "xdata_after_fill.rds"))    
xdata <- readRDS(here::here("results", "xdata_after_fill.rds"))


# 7. Load the final xdata object (already processed)
# (this file already has groupChromPeaks and fillChromPeaks done)

xdata <- readRDS(here::here("results", "xdata_after_fill.rds"))

# 7a. Quick check that everything is ready

hasChromPeaks(xdata)   # should return TRUE
hasFeatures(xdata)     # should return TRUE

# 7b. Export feature intensity table (wide format: features × samples) 
# good for normalization, PCA, etc.

feature_table <- featureValues(xdata, value = "into", method = "medret")
write.csv(feature_table,
          here::here("data", "processed", "feature_table_neg.csv"),
          row.names = TRUE)

# 7c. Export feature definitions (mz, RT, npeaks, etc.) 
# remove list columns because write.csv cannot save lists

fdef <- as.data.frame(featureDefinitions(xdata))
fdef_simple <- fdef[, sapply(fdef, function(x) !is.list(x))] %>%
  tibble::rownames_to_column("feature_id")

write.csv(fdef_simple,
          here::here("data", "processed", "feature_definitions.csv"),
          row.names = FALSE)

# 7d. Export tidy table (long format: feature_id, sample, intensity, mz, RT...)
# good for statistics, clustering, and plotting

ft_df <- as.data.frame(feature_table) %>%
  tibble::rownames_to_column("feature_id")

tidy_tbl <- ft_df %>%
  pivot_longer(cols = -feature_id, names_to = "sample", values_to = "intensity") %>%
  dplyr::left_join(fdef_simple, by = "feature_id")

write.csv(tidy_tbl,
          here::here("data", "processed", "tidy_feature_table.csv"),
          row.names = FALSE)
