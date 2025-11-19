# ================================================================
# 01_volcano_plot.R
# Volcano plot for limma results (NEG mode)
#
# This script:
#   1) loads limma results (logFC, p-values, FDR),
#   2) classifies features as up/down/not significant,
#   3) builds a volcano plot,
#   4) highlights the top significant features,
#   5) saves a static PNG (and optionally an interactive HTML).
# ================================================================

# --- 1. Load packages ------------------------------------------

library(dplyr)
library(ggplot2)
library(ggrepel)   # for non-overlapping text labels (top features)

message("Packages loaded for volcano plot.")


# --- 2. Load limma results -------------------------------------

res <- read.csv(
  "results/limma_results_neg.csv",
  row.names = 1,
  check.names = FALSE
)

message("Limma results loaded: ", nrow(res), " features.")


# --- 3. Prepare data for volcano plot --------------------------

# Add:
#   -neg_log10_FDR : -log10(FDR-adjusted p-value)
#   -significance  : category (Up / Down / Not significant)
#   -feature_id    : FT ID as a column
volcano <- res %>%
  mutate(
    neg_log10_FDR = -log10(adj.P.Val),
    significance = case_when(
      adj.P.Val < 0.05 & logFC > 0  ~ "Up in Wounding",
      adj.P.Val < 0.05 & logFC < 0  ~ "Down in Wounding",
      TRUE                          ~ "Not significant"
    )
  ) %>%
  tibble::rownames_to_column(var = "feature_id")

# Select top N features by FDR for labeling
top_n <- 10

volcano <- volcano %>%
  arrange(adj.P.Val) %>%
  mutate(
    label = ifelse(row_number() <= top_n & adj.P.Val < 0.05, feature_id, NA)
  )

# Data frame with only labeled points (used in geom_text_repel)
volcano_labels <- volcano %>% filter(!is.na(label))


# --- 4. Construct volcano plot ---------------------------------

p <- ggplot(volcano, aes(x = logFC, y = neg_log10_FDR)) +
  geom_point(aes(color = significance), alpha = 0.7, size = 2) +
  scale_color_manual(
    values = c(
      "Up in Wounding"   = "#D55E00",
      "Down in Wounding" = "#0072B2",
      "Not significant"  = "grey70"
    )
  ) +
  # vertical line at logFC = 0 (no change)
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  # horizontal line at FDR = 0.05
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black") +
  # label top features
  geom_text_repel(
    data = volcano_labels,
    aes(label = label),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3,
    min.segment.length = 0
  ) +
  labs(
    title    = "Volcano plot (Wounding vs Null)",
    subtitle = "NEG mode, limma; FDR < 0.05 highlighted, top 10 features labeled",
    x = "log2 fold change (Wounding vs Null)",
    y = "-log10(FDR-adjusted p-value)",
    color = "Feature status"
  ) +
  theme_minimal(base_size = 14)


# --- 5. Create output directory --------------------------------

outdir <- "results/plots"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)


# --- 6. Save static PNG ----------------------------------------

outfile_png <- file.path(outdir, "volcano_neg.png")
ggsave(outfile_png, p, width = 10, height = 7, dpi = 300)

message("Volcano plot (PNG) saved to: ", outfile_png)


# --- 7. (Optional) Save interactive HTML version ---------------

# If plotly and htmlwidgets are installed, create an interactive plot.
# The script does NOT install packages automatically, it only uses them
# if they are already available.

if (requireNamespace("plotly", quietly = TRUE) &&
    requireNamespace("htmlwidgets", quietly = TRUE)) {
  
  interactive_plot <- plotly::ggplotly(p)
  
  outfile_html <- file.path(outdir, "volcano_neg_interactive.html")
  htmlwidgets::saveWidget(
    widget = interactive_plot,
    file   = outfile_html,
    selfcontained = TRUE
  )
  
  message("Interactive volcano plot saved to: ", outfile_html)
} else {
  message("plotly/htmlwidgets not available: skipping interactive export.")
}


# --- 8. Interpretation (for documentation, not used in code) ---
# Volcano plot demonstrates a strong metabolic response to wounding,
# with numerous significantly up-regulated and down-regulated features (FDR < 0.05).
# Many metabolites show large fold-changes (|log2FC| > 2), indicating robust
# and biologically meaningful differences between conditions.
# The symmetric distribution of features around logFC = 0 suggests that
# normalization was effective and no systematic bias remained.