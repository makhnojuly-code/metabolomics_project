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

# загрузка чекпойнта 4
xdata <- readRDS(here::here("results","xdata_s4_peaks.rds"))

# временно в одиночный поток
old_bp <- BiocParallel::bpparam()
BiocParallel::register(BiocParallel::SerialParam())

# шаг 5 — RT alignment
xdata <- adjustRtime(xdata, param = ObiwarpParam(centerSample = 130))
saveRDS(xdata, here::here("results","xdata_s5_rtime.rds"))

# вернуть прежний backend (если нужен)
BiocParallel::register(old_bp)


# 6. Feature grouping (align peaks across all samples)
xdata <- groupChromPeaks(xdata, param = PeakDensityParam(sampleGroups = factor(rep("all", length(fileNames(xdata))))))
saveRDS(xdata, here::here("results", "xdata_s6_grouped.rds")) 


# 6.5. Fill missing peaks by integration over RT–m/z windows
xdata <- fillChromPeaks(xdata)                                   
saveRDS(xdata, here::here("results", "xdata_after_fill.rds"))    

# 7. Export feature table (peaks × samples)
feature_table <- featureValues(xdata, value="into", method="medret")
write.csv(feature_table, "data/processed/feature_table_neg.csv")

# Check the list of files (samples) used in the analysis
fileNames(xdata)

