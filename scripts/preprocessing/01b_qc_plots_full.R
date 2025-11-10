
source("scripts/preprocessing/00_setup.R")


# scripts/preprocessing/01b_qc_plots_full.R
# Full QC plots (fast & colorful, pastel palette)
# RT-shift + full BPC (all traces) + TIC boxplot



#  Settings 
# set path to the xdata file (this file already has aligned retention times)
xdata_rds <- here::here("data", "intermediate", "xdata_rtime.rds")  # after RT alignment
stopifnot(file.exists(xdata_rds))
# folder for saving QC plots
out_dir   <- here::here("results", "qc")
if (!dir.exists(out_dir)) dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Cross-platform parallel backend
bp <- if (.Platform$OS.type == "windows") {
  SnowParam(workers = max(1, parallel::detectCores() - 1))
} else {
  MulticoreParam(workers = max(1, parallel::detectCores() - 1))
}

# Soft pastel colors
palette_soft <- function(n, alpha = 0.8) {
  cols <- c(brewer.pal(8, "Pastel1"), brewer.pal(8, "Set2"))
  cols <- rep(cols, length.out = n)
  adjustcolor(cols, alpha.f = alpha)
}

# Bright vivid palette (higher chroma, evenly spaced hues)
palette_bright <- function(n, alpha = 0.85) {
  hues <- seq(15, 375, length.out = n + 1L)
  cols <- grDevices::hcl(h = hues[1:n], c = 85, l = 50)
  grDevices::adjustcolor(cols, alpha.f = alpha)
}

# PNG device (uses ragg if available)
open_png <- function(path, width = 2200, height = 1400, res = 160) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(path, width = width, height = height, units = "px", res = res)
  } else {
    png(path, width = width, height = height, res = res, type = "cairo")
  }
}

#  Universal Chromatogram(s) -> data.frame 
# Works for: Chromatogram, Chromatograms, or list of Chromatogram
chrom_to_df <- function(CH) {
  # Single Chromatogram
  if (inherits(CH, "Chromatogram")) {
    return(data.table::data.table(
      file = 1L,
      rt = as.numeric(MSnbase::rtime(CH)),
      intensity = as.numeric(MSnbase::intensity(CH))
    ))
  }
  # Matrix-like Chromatograms (columns ~= files)
  if (inherits(CH, "Chromatograms")) {
    out <- vector("list", ncol(CH))
    for (i in seq_len(ncol(CH))) {
      ch_i <- CH[, i, drop = FALSE]                # keep Chromatograms class
      rt_i <- MSnbase::rtime(ch_i)
      it_i <- MSnbase::intensity(ch_i)
      out[[i]] <- data.table::data.table(
        file = i,
        rt = as.numeric(unlist(rt_i, use.names = FALSE)),
        intensity = as.numeric(unlist(it_i, use.names = FALSE))
      )
    }
    return(data.table::rbindlist(out, use.names = TRUE))
  }
  # List of Chromatogram
  if (is.list(CH) && length(CH) > 0 &&
      all(vapply(CH, inherits, logical(1), what = "Chromatogram"))) {
    out <- lapply(seq_along(CH), function(i) {
      data.table::data.table(
        file = i,
        rt = as.numeric(MSnbase::rtime(CH[[i]])),
        intensity = as.numeric(MSnbase::intensity(CH[[i]]))
      )
    })
    return(data.table::rbindlist(out, use.names = TRUE))
  }
  stop("chrom_to_df(): unsupported object class: ", paste(class(CH), collapse = "/"))
}

# Downsample to ~1 point per pixel along X (per file), keeping peak shape
# fun = "max" for BPC, "sum" for TIC
downsample_1px <- function(df, px = 2200, fun = c("max","sum")) {
  stopifnot(all(c("file","rt","intensity") %in% names(df)))
  fun <- match.arg(fun)
  setDT(df)
  df <- df[is.finite(rt) & is.finite(intensity)][order(file, rt)]
  rng <- df[, range(rt, na.rm = TRUE)]
  bin <- (rng[2] - rng[1]) / px
  if (!is.finite(bin) || bin <= 0) return(df[])
  df[, bin_id := floor((rt - rng[1]) / bin)]
  if (fun == "max") {
    out <- df[, .(rt = mean(rt, na.rm = TRUE),
                  intensity = max(intensity, na.rm = TRUE)),
              by = .(file, bin_id)]
  } else {
    out <- df[, .(rt = mean(rt, na.rm = TRUE),
                  intensity = sum(intensity, na.rm = TRUE)),
              by = .(file, bin_id)]
  }
  setorder(out, file, rt)
  out[]
}

#  1) Load data 
xdata <- readRDS(xdata_rds)

#  2) RT-shift plot (fast)
open_png(file.path(out_dir, "qc_adjusted_rtime_full.png"), width = 1800, height = 1100, res = 160)
par(mar = c(4,4,1,1))
plotAdjustedRtime(xdata, lwd = 2, col = "#444444")
mtext("Adjusted retention time (all files)", side = 3, line = -1.2, cex = 1.0)
dev.off()

#  3) BPC (all traces) — FAST & COLOR
# Extract BPC per file (MS1 only), in parallel
bpc_ch <- chromatogram(filterMsLevel(xdata, 1L),
                       aggregationFun = "max",
                       BPPARAM = bp)

# Chromatogram(s) -> data.frame (file, rt, intensity)
bpc_df <- chrom_to_df(bpc_ch)
bpc_df$file      <- as.character(bpc_df$file)
bpc_df$rt        <- as.numeric(bpc_df$rt)
bpc_df$intensity <- as.numeric(bpc_df$intensity)
bpc_df           <- bpc_df[is.finite(bpc_df$intensity), c("file","rt","intensity")]

# Downsample to 1 point / pixel
bpc_ds <- downsample_1px(bpc_df, px = 2200, fun = "max")

files_bpc <- unique(bpc_ds$file)
# <<< яркая палитра для BPC
cols_bpc  <- palette_bright(length(files_bpc), alpha = 0.8); names(cols_bpc) <- files_bpc

open_png(file.path(out_dir, "qc_bpc_full.png"), width = 2200, height = 1400, res = 160)
par(mar = c(4,4,1,1))
plot(NA,
     xlim = range(bpc_ds$rt, na.rm = TRUE),
     ylim = range(bpc_ds$intensity, na.rm = TRUE),
     xlab = "RT (sec)", ylab = "Intensity",
     main = sprintf("BPC after RT alignment (all %d traces)", length(files_bpc)))
invisible(lapply(files_bpc, function(f) {
  d <- bpc_ds[bpc_ds$file == f]
  lines(d$rt, d$intensity, col = cols_bpc[f], lwd = 1.6)   # thicker line
}))
legend("topright", legend = files_bpc, col = cols_bpc, lwd = 1, cex = 0.7, bty = "n")
dev.off()

# 3a) OPTIONAL: BPC in blocks of 12 files (easier to see, less lines per plot) 
# this will create multiple images: qc_bpc_block_XX.png (if you don't need this -> comment out)
files_bpc <- unique(bpc_ds$file)
chunks <- split(files_bpc, ceiling(seq_along(files_bpc)/12))
for (i in seq_along(chunks)) {
  fs <- chunks[[i]]
  open_png(file.path(out_dir, sprintf("qc_bpc_block_%02d.png", i)),
           width = 2200, height = 1400, res = 160)
  par(mar = c(4,4,1,1))
  plot(NA,
       xlim = range(bpc_ds$rt, na.rm = TRUE),
       ylim = range(bpc_ds$intensity, na.rm = TRUE),
       xlab = "RT (sec)", ylab = "Intensity",
       main = sprintf("BPC after RT alignment — files %s",
                      paste(fs, collapse = ", ")))
  invisible(lapply(fs, function(f) {
    d <- bpc_ds[bpc_ds$file == f]
    lines(d$rt, d$intensity, col = cols_bpc[f], lwd = 1.6)
  }))
  legend("topright", legend = fs, col = cols_bpc[fs], lwd = 1, cex = 0.7, bty = "n")
  dev.off()
}

#  4) TIC boxplot (FAST, per-file distribution) 
tic_ch <- chromatogram(filterMsLevel(xdata, 1L),
                       aggregationFun = "sum",
                       BPPARAM = bp)

# convert chromatograms to data.frame (file, rt, intensity)
tic_df <- chrom_to_df(tic_ch)
tic_df$file      <- as.character(tic_df$file)
tic_df$rt        <- as.numeric(tic_df$rt)
tic_df$intensity <- as.numeric(tic_df$intensity)
tic_df           <- tic_df[is.finite(tic_df$intensity), c("file","rt","intensity")]

# choose color per file (here pastel palette)
files_tic <- sort(unique(tic_df$file))
cols_tic  <- palette_soft(length(files_tic), alpha = 0.8)  
names(cols_tic) <- files_tic
bx_cols <- cols_tic[files_tic]

# list of vectors for boxplot (one vector per file)
tic_split <- split(tic_df$intensity, tic_df$file)

open_png(file.path(out_dir, "qc_tic_full.png"), width = 1600, height = 1000, res = 160)
par(mar = c(6,5,2,1))
boxplot(tic_split,
        col     = bx_cols,
        border  = "grey50",   # neutral mid-grey border
        outline = TRUE,
        las     = 2,
        xlab    = "",
        ylab    = "TIC (sum of intensities per scan)",
        main    = "TIC per file (after RT)")
mtext("File index", side = 1, line = 4)
dev.off()

cat("FAST FULL QC plots saved to:", normalizePath(out_dir), "\n")




