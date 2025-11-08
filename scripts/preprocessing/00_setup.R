# ======================================================
# 00_packages.R
# Load and configure all required packages for LC–MS workflow
# ======================================================

# 1. Install pacman if not installed
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)

# 2. Load all required packages
pacman::p_load(
  conflicted,     # handle function name conflicts
  tidyverse,      # data handling and wrangling
  here,           # project paths
  MSnbase,        # LC–MS data import
  xcms,           # LC–MS preprocessing
  pheatmap,       # visualization (heatmaps)
  RColorBrewer,   # color palettes
  pander,         # reporting
  MsExperiment,   # data container for MS
  BiocManager,    # Bioconductor manager
  BiocParallel    # parallel processing backend
)



# 3. Register parallel backend (auto-select based on OS)
# Detect the operating system and choose the appropriate backend for parallel processing
# On Windows: use SnowParam (MulticoreParam is not supported)
# On Linux / macOS / WSL: use MulticoreParam (faster)
if (.Platform$OS.type == "windows") {
  # Use one less than the total number of CPU cores to keep the system responsive
  register(SnowParam(workers = parallel::detectCores() - 1))
} else {
  # Same idea for Linux/macOS — use multicore processing for better performance
  register(MulticoreParam(workers = max(1, parallel::detectCores() - 1)))
}

# 4. Resolve function conflicts to avoid ambiguity
conflicts_prefer(
  dplyr::filter,
  dplyr::rename,
  xcms::groups,
  xcms::collect,
  MSnbase::combine
)

# 5. Confirmation message
message("✅ All basic packages loaded and parallel backend successfully registered.")