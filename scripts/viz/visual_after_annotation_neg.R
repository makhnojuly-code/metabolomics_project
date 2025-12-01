###############################################################
# visual_after_annotation_neg.R
# NEG mode – visualisation after annotation
###############################################################
# Run this script from the project root so that relative paths work:
#  - data/processed/
#  - results/annotation/
#  - results/stress_core_neg/

# Packages ---------------------------------------------------
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(viridis)
library(pheatmap)
library(conflicted)


library(ggrepel)    # pretty labels on volcano
library(ggridges)   # ridgeplots
library(circlize)   # chord diagram
library(ggpubr)
library(reshape2)
library(gridExtra)

conflicts_prefer(dplyr::mutate)
conflicts_prefer(dplyr::summarise)
conflicts_prefer(dplyr::arrange)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::group_by)

# Input files ------------------------------------------------
core_path   <- "results/stress_core_neg/core_metabolites_neg.csv"
enrich_path <- "results/annotation/mummicho_pathway_enrichment.csv"

# Output folders ---------------------------------------------
viz_base <- "results/viz_neg"
vol_dir  <- file.path(viz_base, "volcano")
path_dir <- file.path(viz_base, "pathways")
hm_dir   <- file.path(viz_base, "heatmaps")
box_dir  <- file.path(viz_base, "boxplots")
chord_dir <- file.path(viz_base, "chord")

dir.create(vol_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(path_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(hm_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(box_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(chord_dir, recursive = TRUE, showWarnings = FALSE)

# Load data --------------------------------------------------
core       <- read_csv(core_path,   show_col_types = FALSE)
enrich_raw <- read_csv(enrich_path, show_col_types = FALSE)
names(enrich_raw)[1] <- "pathway"

enrich <- enrich_raw %>%
  rename(
    total_size  = Total_Size,
    hits        = Hits,
    sig_hits    = Sig_Hits,
    mummichog_p = Mummichog_Pvals,
    gsea_p      = GSEA_Pvals,
    combined_p  = Combined_Pvals
  ) %>%
  mutate(
    rich_factor   = sig_hits / total_size,
    neg_log10_p   = -log10(combined_p),
    pathway_clean = str_trunc(pathway, 60)
  )

message("Loaded core metabolites: ", nrow(core))
message("Loaded KEGG enrichment rows: ", nrow(enrich))

###############################################################
# 1) Volcano plot
###############################################################

# Volcano plot (NEG mode)
# This figure shows how core metabolites change under stress.
# Each point is one metabolite.
# X-axis: log2 fold change (stress vs control).
# Y-axis: -log10 adjusted p-value.
# Red points = metabolites increased under stress.
# Blue points = metabolites decreased under stress.
# Only significant metabolites are shown in this core dataset.
# Input data: core_metabolites_neg.csv


p_cut  <- 0.05     # adjusted p-value cutoff
fc_cut <- 1.0      # |log2FC| >= 1

vol_df <- core %>%
  dplyr::mutate(
    neg_log10_p = -log10(min_adjP),
    regulation = dplyr::case_when(
      min_adjP <= p_cut & mean_logFC >=  fc_cut ~ "Up in stress",
      min_adjP <= p_cut & mean_logFC <= -fc_cut ~ "Down in stress",
      TRUE                                      ~ "Not significant"
    ),
    regulation = factor(
      regulation,
      levels = c("Down in stress", "Up in stress", "Not significant")
    )
  )

vol_cols <- c(
  "Down in stress"   = "#377EB8",
  "Up in stress"     = "#E41A1C",
  "Not significant"  = "grey80"
)

p_volcano <- ggplot(
  vol_df,
  aes(x = mean_logFC, y = neg_log10_p, colour = regulation)
) +
  geom_point(size = 2.8, alpha = 0.9) +
  scale_colour_manual(values = vol_cols, name = "Regulation") +
  geom_vline(xintercept = c(-fc_cut, fc_cut),
             linetype = "dashed", colour = "grey60", linewidth = 0.5) +
  geom_hline(yintercept = -log10(p_cut),
             linetype = "dashed", colour = "grey60", linewidth = 0.5) +
  labs(
    title = "Volcano plot — NEG core metabolites",
    x = "log2 fold change (stress / control)",
    y = "-log10(min adjusted p-value)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

message("Showing updated volcano plot...")
print(p_volcano)

ggsave(
  file.path(vol_dir, "neg_core_volcano_scientific.png"),
  p_volcano, width = 7, height = 5.5, dpi = 320
)
ggsave(
  file.path(vol_dir, "neg_core_volcano_scientific.pdf"),
  p_volcano, width = 7, height = 5.5
)

# Figure: Volcano plot — NEG core metabolites.
# This volcano plot summarizes how strongly each core metabolite responds to wounding stress.
# Each point represents one metabolite. 
# The x-axis shows the log2 fold change (stress vs control), so metabolites on the right are
# higher in stress, and metabolites on the left are lower in stress.
# The y-axis shows the statistical significance (−log10 of the adjusted p-value):
# points higher on the plot are more significant.
# Red points are metabolites that are up-regulated in stress,
# blue points are metabolites that are down-regulated.
# Together, the plot highlights a clear core stress-response pattern with several strongly
# induced metabolites and a smaller group that is suppressed.

###############################################################
# 2) KEGG bubble plot 
###############################################################

# KEGG pathway enrichment (NEG mode)
# This plot shows the top enriched KEGG pathways.
# Rich factor = (number of significant metabolites) / (total metabolites in pathway).
# Bubble size = number of significant metabolites.
# Bubble colour = -log10(combined p-value).
# Brighter colour means higher statistical significance.
# Input data: mummicho_pathway_enrichment.csv

enrich_plot <- enrich %>%
  dplyr::filter(
    !is.na(rich_factor),
    !is.na(neg_log10_p),
    sig_hits > 0
  ) %>%
  dplyr::arrange(neg_log10_p) %>%
  dplyr::slice_tail(n = 20) %>%
  dplyr::mutate(
    pathway_clean = factor(pathway_clean, levels = pathway_clean)
  )

p_rich <- ggplot(
  enrich_plot,
  aes(
    x      = rich_factor,
    y      = pathway_clean,
    size   = sig_hits,
    colour = neg_log10_p
  )
) +
  geom_point(alpha = 0.9) +
  scale_colour_viridis(
    option = "A",
    name = expression(-log[10](combined~p))
  ) +
  scale_size(range = c(3, 11), name = "Significant hits") +
  labs(
    title = "KEGG pathway enrichment (NEG)",
    x = "Rich factor (Sig_Hits / Total_Size)",
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.text.y      = element_text(size = 11),
    panel.grid.minor = element_blank(),
    legend.position  = "right"
  )

message("Showing KEGG enrichment (rich factor) bubble...")
print(p_rich)

ggsave(
  file.path(path_dir, "neg_kegg_richfactor.png"),
  p_rich, width = 7, height = 6, dpi = 320
)
ggsave(
  file.path(path_dir, "neg_kegg_richfactor.pdf"),
  p_rich, width = 7, height = 6
)

# Figure: KEGG pathway enrichment (NEG mode).
# Bubble plot showing enriched KEGG pathways.
# Rich factor (x-axis) indicates the proportion of significant metabolites in each pathway.
# Bubble size corresponds to the number of significant hits, and color indicates –log10(combined p-value), where brighter colors represent higher significance.

# The enriched pathways show that the plant activates its main stress-response systems.
# The strongest changes appear in flavonoid and anthocyanin biosynthesis, amino acid metabolism, and lipid metabolism.
# These pathways are typically involved in protection against oxidative stress, membrane stabilization, and the production of defensive metabolites.
# Overall, the profile suggests a clear biochemical defense response to stress.

###############################################################
# 3) Heatmap: top significant annotated metabolites (NEG mode)
###############################################################

# Heatmap helps to quickly see patterns in metabolite behavior across samples.
# It shows which metabolites increase or decrease in each group and whether
# samples cluster together based on similar metabolic profiles.

# Input:
#   - feature_table_neg_processed_for_stats_filtered.csv
# Processed intensity matrix for NEG mode (features × samples)
#   - sample_sheet_neg_filtered.csv
# Metadata for samples (sample names and treatment groups)
# - core_features_neg.csv
# Table of significant metabolites with statistics
#(feature_id, adjusted p-values, KEGG names, annotations)


expr_path          <- "data/processed/feature_table_neg_processed_for_stats_filtered.csv"
meta_path          <- "data/processed/sample_sheet_neg_filtered.csv"
core_features_path <- "results/annotation/core_features_neg.csv"

top_n <- 30
message("Heatmap: top ", top_n, " annotated metabolites.")

expr_df       <- read.csv(expr_path,          check.names = FALSE)
meta          <- read.csv(meta_path,          check.names = FALSE)
core_features <- read.csv(core_features_path, check.names = FALSE)

feature_ids <- expr_df[[1]]
expr_mat    <- as.matrix(expr_df[, -1])
rownames(expr_mat) <- feature_ids

stopifnot(all(meta$file_name %in% colnames(expr_mat)))
expr_mat <- expr_mat[, meta$file_name]

sig_core <- core_features %>%
  dplyr::arrange(adj.P.Val) %>%
  dplyr::slice_head(n = top_n)

sel_ids  <- sig_core$feature_id
heat_mat <- expr_mat[sel_ids, , drop = FALSE]

message("Heatmap will show ", nrow(heat_mat), " metabolites.")

labels_raw <- ifelse(
  !is.na(sig_core$kegg_name) & sig_core$kegg_name != "",
  sig_core$kegg_name,
  sig_core$Matched.Compound
)

labels_clean <- labels_raw %>%
  sub(";.*$", "", .) %>%
  sub("^\\([^\\)]+\\)", "", .) %>%
  trimws() %>%
  stringr::str_to_sentence() %>%
  stringr::str_trunc(width = 45)

rownames(heat_mat) <- labels_clean

annotation_col <- data.frame(
  Group = factor(meta$group)
)
rownames(annotation_col) <- meta$file_name

ann_colors <- list(
  Group = c(
    "Null"     = "#3CB371",
    "Wounding" = "#FFFF00"
  )
)

my_colors <- colorRampPalette(c("#0000CD", "#FFFFFF", "#CD0000"))(200)

dir.create(hm_dir, showWarnings = FALSE, recursive = TRUE)

out_png <- file.path(hm_dir, paste0("neg_nature_heatmap_top", top_n, ".png"))
out_pdf <- file.path(hm_dir, paste0("neg_nature_heatmap_top", top_n, ".pdf"))

message("Showing heatmap preview in console...")

heat_obj <- pheatmap(
  mat                      = heat_mat,
  color                    = my_colors,
  scale                    = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "ward.D2",
  show_rownames            = TRUE,
  show_colnames            = FALSE,
  fontsize_row             = 7,
  annotation_col           = annotation_col,
  annotation_colors        = ann_colors,
  main                     = paste0("Top ", top_n, " NEG features (Z-scored)")
)
print(heat_obj)

message("Saving PNG and PDF versions...")

pheatmap(
  mat                      = heat_mat,
  color                    = my_colors,
  scale                    = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "ward.D2",
  show_rownames            = TRUE,
  show_colnames            = FALSE,
  fontsize_row             = 7,
  annotation_col           = annotation_col,
  annotation_colors        = ann_colors,
  main                     = NA,
  filename                 = out_png,
  width                    = 6,
  height                   = 7
)

pheatmap(
  mat                      = heat_mat,
  color                    = my_colors,
  scale                    = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "ward.D2",
  show_rownames            = TRUE,
  show_colnames            = FALSE,
  fontsize_row             = 7,
  annotation_col           = annotation_col,
  annotation_colors        = ann_colors,
  main                     = NA,
  filename                 = out_pdf,
  width                    = 6,
  height                   = 7
)

message("Saved heatmap PNG to: ", out_png)
message("Saved heatmap PDF to: ", out_pdf)

# Figure: Heatmap of top significant metabolites (NEG mode).
# This heatmap displays the metabolites that change the most between the 
# Null and Wounding groups. Each row represents one metabolite 
# and each column represents one sample. Colors show Z-scored intensity values,
# where red means higher levels, blue means lower levels, and white indicates average 
# abundance. Because each metabolite is scaled separately, the heatmap clearly 
# shows the patterns of up- and down-regulation across the two groups.

###############################################################
# 4) Boxplots for top 16 significant metabolites (NEG mode)
###############################################################

# Boxplots of the most significant NEG metabolites.
# Each panel (M1–M16) shows one metabolite.
# Input:
#   - expr_mat: processed intensity matrix (NEG mode),
#   - meta: sample information with group (Null / Wounding),
#   - core_features: statistics table (logFC, adj.P.Val, annotation).
# Only metabolites with good data in both groups and with annotation are used.
# Intensities are transformed to log10(intensity + pseudo)
# so very small values are stable and can be compared between groups.

n_total_box <- 16
message("Selecting mixed UP + DOWN metabolites (max ", n_total_box,
        ") with valid data in BOTH groups...")

if ("logFC" %in% names(core_features)) {
  core_features$fc_val <- core_features$logFC
} else if ("mean_logFC" %in% names(core_features)) {
  core_features$fc_val <- core_features$mean_logFC
} else {
  stop("Could not find a fold-change column (logFC / mean_logFC) in core_features.")
}

sig_core <- core_features %>% arrange(adj.P.Val)

valid_features <- sig_core$feature_id[
  sapply(sig_core$feature_id, function(fid) {
    vals <- expr_mat[fid, ]
    g1 <- vals[meta$group == "Null"]
    g2 <- vals[meta$group == "Wounding"]
    sum(g1 > 0, na.rm = TRUE) > 1 && sum(g2 > 0, na.rm = TRUE) > 1
  })
]

sig_valid <- sig_core %>%
  dplyr::filter(feature_id %in% valid_features) %>%
  dplyr::filter(
    (!is.na(kegg_name)        & kegg_name        != "") |
      (!is.na(Matched.Compound) & Matched.Compound != "")
  )

message("Valid annotated features (both groups): ", nrow(sig_valid))

sig_up   <- sig_valid %>% filter(fc_val > 0)
sig_down <- sig_valid %>% filter(fc_val < 0)

message("UP-regulated features:   ", nrow(sig_up))
message("DOWN-regulated features: ", nrow(sig_down))

target_down <- min(ceiling(n_total_box / 2), nrow(sig_down))
target_up   <- min(n_total_box - target_down, nrow(sig_up))

if (target_up + target_down < n_total_box) {
  remaining <- n_total_box - (target_up + target_down)
  extra_up   <- max(0, min(remaining, nrow(sig_up)   - target_up))
  extra_down <- max(0, min(remaining - extra_up,
                           nrow(sig_down) - target_down))
  target_up   <- target_up   + extra_up
  target_down <- target_down + extra_down
}

message("Final target: ", target_up, " UP and ", target_down, " DOWN.")

sig_mix <- bind_rows(
  sig_up   %>% slice_head(n = target_up)   %>% mutate(direction = "Up in stress"),
  sig_down %>% slice_head(n = target_down) %>% mutate(direction = "Down in stress")
) %>%
  arrange(adj.P.Val)

message("Total selected for boxplots: ", nrow(sig_mix))

sig_mix$label_plot <- paste0("M", seq_len(nrow(sig_mix)))

box_df <- as.data.frame(expr_mat[sig_mix$feature_id, , drop = FALSE])
box_df$feature_id <- rownames(box_df)

box_long <- reshape2::melt(
  box_df,
  id.vars       = "feature_id",
  variable.name = "file_name",
  value.name    = "intensity"
) %>%
  left_join(meta[, c("file_name", "group")], by = "file_name") %>%
  left_join(sig_mix[, c("feature_id", "label_plot", "direction", "fc_val")],
            by = "feature_id") %>%
  dplyr::filter(intensity > 0) %>%
  dplyr::filter(!is.na(label_plot))

box_long$group <- factor(box_long$group, levels = c("Null", "Wounding"))

pseudo <- min(box_long$intensity, na.rm = TRUE) / 2
box_long$log10_intensity <- log10(box_long$intensity + pseudo)

p_table <- box_long %>%
  dplyr::group_by(label_plot) %>%
  dplyr::summarise(
    p_value = tryCatch(
      wilcox.test(log10_intensity ~ group)$p.value,
      error = function(e) NA_real_
    ),
    y_pos = max(log10_intensity, na.rm = TRUE) + 0.15,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    stars = dplyr::case_when(
      is.na(p_value)       ~ "",
      p_value < 1e-4       ~ "****",
      p_value < 1e-3       ~ "***",
      p_value < 1e-2       ~ "**",
      p_value < 0.05       ~ "*",
      TRUE                 ~ ""
    )
  )

box_long$label_plot <- factor(
  box_long$label_plot,
  levels = sig_mix$label_plot
)
p_table$label_plot <- factor(
  p_table$label_plot,
  levels = sig_mix$label_plot
)

group_colors <- c("Null" = "#4DA3FF", "Wounding" = "#FFAA00")

p_mixed <- ggplot(box_long, aes(group, log10_intensity, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.6) +
  facet_wrap(~ label_plot, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = group_colors) +
  labs(
    x = NULL,
    y = "log10(intensity + pseudo)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text      = element_text(size = 10, face = "bold"),
    legend.position = "none",
    plot.title      = element_blank()
  ) +
  geom_segment(
    data = p_table,
    aes(x = 1, xend = 2, y = y_pos, yend = y_pos),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = p_table,
    aes(x = 1.5, y = y_pos + 0.05, label = stars),
    inherit.aes = FALSE,
    size = 3
  )

print(p_mixed)

# Save boxplots ----------------------------------------------
box_dir <- file.path(viz_base, "boxplots")
dir.create(box_dir, recursive = TRUE, showWarnings = FALSE)

png_file <- file.path(
  box_dir,
  "neg_boxplots_top16_updown_log10_pseudo.png"
)
pdf_file <- file.path(
  box_dir,
  "neg_boxplots_top16_updown_log10_pseudo.pdf"
)

ggsave(
  filename = png_file,
  plot     = p_mixed,
  width    = 10,
  height   = 7,
  dpi      = 300
)
ggsave(
  filename = pdf_file,
  plot     = p_mixed,
  width    = 10,
  height   = 7
)

message("Boxplots (UP + DOWN, n = 16) saved to:")
message("PNG: ", png_file)
message("PDF: ", pdf_file)

# Figure: Boxplots of 16 key NEG metabolites (M1–M16).
# For each metabolite, blue boxes show Null samples and orange boxes show Wounding samples.
# The y-axis gives log10(intensity + pseudo); points are individual samples.
# Panels where the Wounding box is higher than Null indicate metabolites increased in stress,
# while higher boxes in Null indicate decreased levels in stress.
# Horizontal brackets with stars mark Wilcoxon test significance
# (* p < 0.05, ** p < 0.01, *** p < 0.001, **** p < 0.0001).

###############################################################
# Mapping table for boxplot panels (M → metabolite)
###############################################################

# This table links each boxplot panel (M1–M16) to the real metabolite.
# For every M label we store: feature ID, annotation, log2 fold change
# and adjusted p-value.
# We save a full version (for methods / supplement) and a short,
# nicely formatted version for figures and main text.

mapping_table <- sig_mix %>%
  dplyr::mutate(
    raw_name = dplyr::if_else(
      !is.na(kegg_name) & kegg_name != "",
      kegg_name,
      Matched.Compound
    ),
    clean_name = raw_name %>%
      sub(";.*$", "", .) %>%
      sub("^\\([^\\)]+\\)", "", .) %>%
      trimws() %>%
      stringr::str_to_sentence() %>%
      stringr::str_trunc(width = 80)
  ) %>%
  dplyr::select(
    M_label      = label_plot,
    feature_id,
    logFC        = fc_val,
    adj.P.Val,
    kegg_name,
    Matched.Compound,
    clean_name
  )

message("Boxplot mapping table (M → metabolite):")
print(mapping_table)

box_dir <- file.path(viz_base, "boxplots")
dir.create(box_dir, showWarnings = FALSE, recursive = TRUE)

full_csv <- file.path(box_dir, "neg_boxplot_mapping_M_labels_full.csv")
write.csv(mapping_table, full_csv, row.names = FALSE)

mapping_pretty <- mapping_table %>%
  dplyr::mutate(
    Direction = dplyr::case_when(
      logFC > 0  ~ "Up in stress",
      logFC < 0  ~ "Down in stress",
      TRUE       ~ "No change"
    ),
    logFC_round    = round(logFC, 2),
    adj.P.Val_fmt  = formatC(adj.P.Val, format = "e", digits = 2)
  ) %>%
  dplyr::select(
    M_label,
    Metabolite = clean_name,
    Direction,
    logFC_round,
    adj.P.Val_fmt
  )

pretty_csv <- file.path(box_dir, "neg_boxplot_mapping_M_labels_pretty.csv")
write.csv(mapping_pretty, pretty_csv, row.names = FALSE)

pdf_file <- file.path(box_dir, "neg_boxplot_mapping_M_labels_pretty.pdf")
pdf(pdf_file, height = 8, width = 10)
gridExtra::grid.table(mapping_pretty)
dev.off()

message("Mapping tables saved:")
message("  Full CSV : ", full_csv)
message("  Pretty CSV: ", pretty_csv)
message("  Pretty PDF: ", pdf_file)

# Table: Mapping of M labels (M1–M16) to annotated metabolites,
# showing the direction of change in stress (Up / Down),
# rounded log2 fold change and formatted adjusted p-values.

###############################################################
# 5) Chord diagram: ALL KEGG pathways vs Up / Down
###############################################################

# A chord diagram is used to visualise connections between two groups of data.
# Here it shows how metabolites (Up or Down in stress) are linked to KEGG pathways.
# This gives an overview of how many metabolites belong to each pathway
# and how strongly each pathway participates in the stress response.
# Input for chord diagrams:
#   - mummichog_pathway_enrichment_integ.csv  (pathway ↔ metabolite hits)
#   - mummicho_hits_neg.csv                   (logFC and regulation per metabolite)



message("Preparing data for chord diagram (ALL pathways vs Up/Down)...")

path_enrich <- read.csv("mummichog_pathway_enrichment_integ.csv")
hits        <- read.csv("results/annotation/mummicho_hits_neg.csv")
# --- 1) Pathway ↔ Empirical.Compound mapping -------------------------
path_map <- path_enrich %>%
  select(Pathway = X, cpd.hits) %>%
  filter(!is.na(cpd.hits), cpd.hits != "") %>%
  separate_rows(cpd.hits, sep = ";") %>%     
  mutate(Empirical.Compound = trimws(cpd.hits)) %>%
  select(Pathway, Empirical.Compound)

# --- 2) LogFC + Up/Down info -----------------------------------------
hits_reg <- hits %>%
  mutate(
    regulation = case_when(
      logFC > 0 ~ "Up in stress",
      logFC < 0 ~ "Down in stress",
      TRUE      ~ "No change"
    )
  ) %>%
  filter(regulation != "No change")

# --- 3) Join pathways with FC info -----------------------------------
path_with_fc <- path_map %>%
  left_join(
    hits_reg %>% select(Empirical.Compound, logFC, regulation),
    by = "Empirical.Compound"
  ) %>%
  filter(!is.na(regulation))

# --- 4) Summary table for chord diagram -------------------------------
chord_df <- path_with_fc %>%
  count(Pathway, regulation, name = "n")

message("Chord: summary table created with ", nrow(chord_df), " rows.")

# --- 5) Sector order --------------------------------------------------
path_order <- chord_df %>%
  distinct(Pathway) %>%
  arrange(Pathway) %>%
  pull()

reg_order    <- c("Up in stress", "Down in stress")
sector_order <- c(path_order, reg_order)

# --- 6) Colors --------------------------------------------------------
grid.col <- c(
  setNames(rep("grey70", length(path_order)), path_order),
  "Up in stress"   = "#1B7837",
  "Down in stress" = "#B2182B"
)

# --- 7) Function to draw the chord -----------------------------------
plot_chord_all <- function() {
  circos.clear()
  par(mar = c(1, 1, 1, 1))
  
  circos.par(
    canvas.xlim = c(-1.15, 1.15),
    canvas.ylim = c(-1.15, 1.15),
    gap.after    = c(rep(2, length(path_order) - 1), 8, 2, 2),
    track.margin = c(0.01, 0.01)
  )
  
  chordDiagram(
    x               = chord_df,
    order           = sector_order,
    grid.col        = grid.col,
    transparency    = 0.25,
    annotationTrack = "grid"
  )
  
  circos.trackPlotRegion(
    track.index = 1,
    bg.border   = NA,
    panel.fun   = function(x, y) {
      lab  <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      circos.text(
        x          = mean(xlim),
        y          = 0,
        labels     = lab,
        cex        = 0.45,
        facing     = "clockwise",
        niceFacing = TRUE,
        adj        = c(0, 0.5)
      )
    }
  )
}

# --- 8) Show in console -----------------------------------------------
message("Rendering chord diagram in console...")
plot_chord_all()

# --- 9) Save to PNG + PDF ---------------------------------------------
png(
  filename = file.path(chord_dir, "neg_chord_ALL_pathways_updown.png"),
  width    = 11,
  height   = 11,
  units    = "in",
  res      = 300
)
plot_chord_all()
dev.off()

pdf(
  file   = file.path(chord_dir, "neg_chord_ALL_pathways_updown.pdf"),
  width  = 10,
  height = 10
)
plot_chord_all()
dev.off()

message("Chord diagram saved:")
message("  PNG: ", file.path(chord_dir, 'neg_chord_ALL_pathways_updown.png'))
message("  PDF: ", file.path(chord_dir, 'neg_chord_ALL_pathways_updown.pdf'))

# Figure: Chord diagram of all enriched KEGG pathways (NEG mode).
# Grey sectors show individual pathways; green and red sectors show
# metabolites that are up- or down-regulated in stress.
# The density of links indicates how strongly each pathway is involved
# in the stress response and helps to choose the main pathways
# for the focused chord plot.

###############################################################
# 6) Chord diagram: main 12 pathways (P1–P12, Up/Down)
###############################################################

# Chord diagram of the 12 most important KEGG pathways.
# Shows how metabolites from each pathway are Up or Down in stress.
# P1–P12 = pathways, green/red = direction of regulation.
# Each ribbon connects a pathway to Up/Down.
# Input files: pathway_enrichment_integ + hits_neg.


# 1) Load input data ---------------------------------------
path_enrich <- read.csv("mummichog_pathway_enrichment_integ.csv")
hits        <- read.csv("results/annotation/mummicho_hits_neg.csv")

# 2) Pathway ↔ Empirical.Compound mapping -----------------
path_map <- path_enrich %>%
  select(Pathway = X, cpd.hits) %>%
  filter(!is.na(cpd.hits), cpd.hits != "") %>%
  separate_rows(cpd.hits, sep = ";") %>%
  mutate(Empirical.Compound = trimws(cpd.hits)) %>%
  select(Pathway, Empirical.Compound)

# 3) Add logFC & regulation to hits ------------------------
hits_reg <- hits %>%
  mutate(
    regulation = case_when(
      logFC > 0 ~ "Up in stress",
      logFC < 0 ~ "Down in stress",
      TRUE      ~ "No change"
    )
  ) %>%
  filter(regulation != "No change")

# 4) Combine pathways + regulation per metabolite ----------
path_with_fc <- path_map %>%
  left_join(hits_reg %>% select(Empirical.Compound, logFC, regulation),
            by = "Empirical.Compound") %>%
  filter(!is.na(regulation))

# 5) Select top pathways (manually defined) -------------
top_paths <- c(
  "alpha-Linolenic acid metabolism",
  "Galactose metabolism",
  "Starch and sucrose metabolism",
  "Valine, leucine and isoleucine biosynthesis",
  "Phenylpropanoid biosynthesis",
  "Inositol phosphate metabolism",
  "Pentose phosphate pathway",
  "Amino sugar and nucleotide sugar metabolism",
  "Carbon fixation by Calvin cycle",
  "Phenylalanine, tyrosine and tryptophan biosynthesis",
  "Flavone and flavonol biosynthesis",
  "Anthocyanin biosynthesis"
)

path_meta <- data.frame(
  Pathway = top_paths,
  ID = paste0("P", seq_along(top_paths)),
  Category = c(
    "Lipids",         # P1
    "Carbohydrates",  # P2
    "Carbohydrates",  # P3
    "Amino acids",    # P4
    "Secondary",      # P5
    "Signaling",      # P6
    "Carbohydrates",  # P7
    "Carbohydrates",  # P8
    "Energy",         # P9
    "Amino acids",    # P10
    "Secondary",      # P11 (flavonoids)
    "Secondary"       # P12 (anthocyanins)
  ),
  stringsAsFactors = FALSE
)

# 6) Metabolite-level mapping for selected pathways --------
metab_map_top <- path_with_fc %>%
  filter(Pathway %in% path_meta$Pathway) %>%
  left_join(path_meta, by = "Pathway") %>%
  select(ID, regulation, Empirical.Compound, Category) %>%
  distinct() %>%
  mutate(weight = 1)

# 7) Colors ------------------------------------------------
path_cols <- c(
  P1  = "#E69F00",
  P2  = "#56B4E9",
  P3  = "#009E73",
  P4  = "#F0E442",
  P5  = "#0072B2",
  P6  = "#D55E00",
  P7  = "#CC79A7",
  P8  = "#999999",
  P9  = "#8DA0CB",
  P10 = "#66C2A5",
  P11 = "#CD3278",
  P12 = "#2F4F4F"
)

reg_cols <- c(
  "Up in stress"   = "#1B7837",
  "Down in stress" = "#B2182B"
)

grid.col  <- c(path_cols, reg_cols)
link_cols <- grid.col[metab_map_top$ID]

path_order   <- path_meta$ID
reg_order    <- c("Up in stress", "Down in stress")
sector_order <- c(path_order, reg_order)

plot_chord_main <- function() {
  circos.clear()
  par(mar = c(5, 5, 5, 5))
  
  circos.par(
    gap.after    = c(rep(8, length(path_order) - 1), 14, 6, 6),
    track.margin = c(0.01, 0.01)
  )
  
  chordDiagram(
    x               = metab_map_top[, c("ID", "regulation", "weight")],
    order           = sector_order,
    grid.col        = grid.col,
    col             = link_cols,
    transparency    = 0.40,
    annotationTrack = "grid",
    link.sort       = TRUE,
    link.lwd        = 1,
    link.border     = NA
  )
  
  circos.trackPlotRegion(
    track.index = 1,
    bg.border   = NA,
    panel.fun   = function(x, y) {
      lab  <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
      
      circos.text(
        x          = mean(xlim),
        y          = ylim[2] + 0.3,
        labels     = lab,
        cex        = 0.7,
        facing     = "outside",
        niceFacing = TRUE
      )
    }
  )
  
  par(xpd = TRUE)
  legend(
    "bottomleft",
    inset  = c(-0.04, 0.10),
    legend = paste0(
      path_meta$ID, " – ", path_meta$Pathway,
      " (", path_meta$Category, ")"
    ),
    pch    = 15,
    pt.cex = 1.2,
    col    = grid.col[path_meta$ID],
    cex    = 0.80,
    bty    = "n",
    title  = "Pathways"
  )
}



# Save PNG
png(
  filename = file.path(chord_dir, "chord_main_pathways.png"),
  width    = 14,
  height   = 8,
  units    = "in",
  res      = 300
)
plot_chord_main()
dev.off()

# Save PDF
pdf(
  file   = file.path(chord_dir, "chord_main_pathways.pdf"),
  width  = 14,
  height = 8
)
plot_chord_main()
dev.off()

# This diagram shows how metabolites from 12 major KEGG pathways are regulated under stress.
# Most pathways send many links to the Up-in-stress sector, showing a broad activation of lipid, 
# carbohydrate and amino-acid metabolism.
# The flavonoid and anthocyanin pathways also contribute several up-regulated
# metabolites, indicating enhanced production of protective phenolic compounds during stress.
# 

message("Done. Figures saved under: ", viz_base)




