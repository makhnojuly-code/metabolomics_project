############################################################
# Metabolomics visualisation (NEG mode) - Shiny app
############################################################

############################
# 0. LIBRARIES
############################
library(shiny)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(viridis)
library(pheatmap)
library(conflicted)
library(ggrepel)
library(ggridges)
library(circlize)
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

############################
# 1. INPUT PATHS & FOLDERS
############################

core_path   <- "results/stress_core_neg/core_metabolites_neg.csv"
enrich_path <- "results/annotation/mummicho_pathway_enrichment.csv"

viz_base  <- "results/viz_neg"
vol_dir   <- file.path(viz_base, "volcano")
path_dir  <- file.path(viz_base, "pathways")
chord_dir <- file.path(viz_base, "chord")

dir.create(vol_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(path_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(chord_dir, recursive = TRUE, showWarnings = FALSE)

############################
# 2. LOAD DATA FOR VOLCANO + KEGG
############################

core        <- read_csv(core_path,   show_col_types = FALSE)
enrich_raw  <- read_csv(enrich_path, show_col_types = FALSE)
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

############################
# 3. BASE VOLCANO OBJECTS
############################

p_cut_default  <- 0.05
fc_cut_default <- 1.0

vol_df <- core %>%
  mutate(
    neg_log10_p = -log10(min_adjP),
    regulation = case_when(
      min_adjP <= p_cut_default & mean_logFC >=  fc_cut_default ~ "Up in stress",
      min_adjP <= p_cut_default & mean_logFC <= -fc_cut_default ~ "Down in stress",
      TRUE                                                      ~ "Not significant"
    ),
    regulation = factor(
      regulation,
      levels = c("Down in stress", "Up in stress", "Not significant")
    )
  )

vol_cols <- c(
  "Down in stress"  = "#377EB8",
  "Up in stress"    = "#E41A1C",
  "Not significant" = "grey80"
)

############################
# 4. PREPARE DATA FOR CHORD
#    (top 10 pathways P1–P10, Up/Down)
############################

# Adjust working directory or paths if needed
path_enrich <- read.csv("mummichog_pathway_enrichment_integ.csv")
hits        <- read.csv("results/annotation/mummicho_hits_neg.csv")

# 4.1 Pathway ↔ Empirical.Compound mapping
path_map <- path_enrich %>%
  select(Pathway = X, cpd.hits) %>%
  filter(!is.na(cpd.hits), cpd.hits != "") %>%
  separate_rows(cpd.hits, sep = ";") %>%
  mutate(Empirical.Compound = trimws(cpd.hits)) %>%
  select(Pathway, Empirical.Compound)

# 4.2 Add logFC and regulation information
hits_reg <- hits %>%
  mutate(
    regulation = case_when(
      logFC > 0 ~ "Up in stress",
      logFC < 0 ~ "Down in stress",
      TRUE      ~ "No change"
    )
  ) %>%
  filter(regulation != "No change")

# 4.3 Join pathways with regulation per metabolite
path_with_fc <- path_map %>%
  left_join(
    hits_reg %>% select(Empirical.Compound, logFC, regulation),
    by = "Empirical.Compound"
  ) %>%
  filter(!is.na(regulation))

# 4.4 Define top 10 pathways and metadata
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
  "Phenylalanine, tyrosine and tryptophan biosynthesis"
)

path_meta <- data.frame(
  Pathway  = top_paths,
  ID       = paste0("P", seq_along(top_paths)),
  Category = c(
    "Lipids",
    "Carbohydrates",
    "Carbohydrates",
    "Amino acids",
    "Secondary",
    "Signaling",
    "Carbohydrates",
    "Carbohydrates",
    "Energy",
    "Amino acids"
  ),
  stringsAsFactors = FALSE
)

# 4.5 Metabolite-level mapping for selected pathways
metab_map_top <- path_with_fc %>%
  filter(Pathway %in% path_meta$Pathway) %>%
  left_join(path_meta, by = "Pathway") %>%
  select(ID, regulation, Empirical.Compound, Category) %>%
  distinct() %>%
  mutate(weight = 1)

# 4.6 Color palettes and sector order -----------------------

# Bright palette: one distinct color per pathway (P1–P10)
path_palette_bright <- c(
  P1  = "#E69F00",
  P2  = "#56B4E9",
  P3  = "#009E73",
  P4  = "#F0E442",
  P5  = "#0072B2",
  P6  = "#D55E00",
  P7  = "#CC79A7",
  P8  = "#999999",
  P9  = "#8DA0CB",
  P10 = "#66C2A5"
)

# Category-based palette: same color for the same pathway category
category_palette <- c(
  "Lipids"         = "#E69F00",
  "Carbohydrates"  = "#56B4E9",
  "Amino acids"    = "#009E73",
  "Secondary"      = "#CC79A7",
  "Signaling"      = "#F0E442",
  "Energy"         = "#D55E00"
)

# Fixed colors for regulation sectors
reg_cols_fixed <- c(
  "Up in stress"   = "#1B7837",
  "Down in stress" = "#B2182B"
)

# Order of sectors on the circle
path_order   <- path_meta$ID
reg_order    <- c("Up in stress", "Down in stress")
sector_order <- c(path_order, reg_order)

# 4.7 Function to draw the chord diagram (extended) ---------
# Arguments:
#   cex_legend   - legend text size
#   palette_mode - color palette for pathways
#   legend_order - order of legend entries

plot_chord_main <- function(
    cex_legend   = 0.80,
    palette_mode = c("Bright (per pathway)", "By category"),
    legend_order = c("By ID", "By pathway name", "By category")
) {
  palette_mode <- match.arg(palette_mode)
  legend_order <- match.arg(legend_order)
  
  circos.clear()
  par(mar = c(5, 5, 5, 5))
  
  # --- choose colors for pathway sectors -------------------
  if (palette_mode == "Bright (per pathway)") {
    path_cols <- path_palette_bright
  } else { # "By category"
    # map each pathway ID to a category color
    path_cols <- setNames(
      category_palette[path_meta$Category],
      path_meta$ID
    )
  }
  
  # full grid.col for all sectors (pathways + regulation)
  grid.col <- c(path_cols, reg_cols_fixed)
  
  # link colors follow pathway colors (one color per pathway ID)
  link_cols <- grid.col[metab_map_top$ID]
  
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
  
  # --- legend: build data frame and re-order ----------------
  legend_df <- path_meta
  
  if (legend_order == "By ID") {
    legend_df <- dplyr::arrange(legend_df, ID)
  } else if (legend_order == "By pathway name") {
    legend_df <- dplyr::arrange(legend_df, Pathway)
  } else if (legend_order == "By category") {
    legend_df <- dplyr::arrange(legend_df, Category, ID)
  }
  
  legend_labels <- paste0(
    legend_df$ID, " – ", legend_df$Pathway,
    " (", legend_df$Category, ")"
  )
  legend_cols <- grid.col[legend_df$ID]
  
  par(xpd = TRUE)
  legend(
    "bottomleft",
    inset  = c(-0.04, 0.10),
    legend = legend_labels,
    pch    = 15,
    pt.cex = 1.2,
    col    = legend_cols,
    cex    = cex_legend,
    bty    = "n",
    title  = "Pathways"
  )
}

############################
# 5. UI
############################

ui <- navbarPage(
  title = "Metabolomics visualisation (NEG mode)",
  
  # 5.1 Volcano tab ----------------------------------------------------
  tabPanel(
    "Volcano",
    sidebarLayout(
      sidebarPanel(
        h4("Volcano settings"),
        sliderInput(
          "vol_p_cut",
          "Adjusted p-value cutoff:",
          min   = 0.001,
          max   = 0.10,
          step  = 0.001,
          value = p_cut_default
        ),
        sliderInput(
          "vol_fc_cut",
          "|log2FC| cutoff:",
          min   = 0.2,
          max   = 2.0,
          step  = 0.1,
          value = fc_cut_default
        ),
        sliderInput(
          "vol_alpha",
          "Point alpha:",
          min   = 0.2,
          max   = 1.0,
          step  = 0.1,
          value = 0.9
        ),
        downloadButton("download_volcano", "Download PNG")
      ),
      mainPanel(
        plotOutput("volcano_plot", height = "600px"),
        br(),
        helpText("Each point is a core metabolite. Red = up in stress, blue = down in stress.")
      )
    )
  ),
  
  # 5.2 KEGG bubble tab ----------------------------------------------
  tabPanel(
    "KEGG bubble",
    sidebarLayout(
      sidebarPanel(
        h4("KEGG settings"),
        sliderInput(
          "kegg_top_n",
          "Number of top pathways:",
          min   = 5,
          max   = 40,
          step  = 5,
          value = 20
        ),
        downloadButton("download_kegg", "Download PNG")
      ),
      mainPanel(
        plotOutput("kegg_plot", height = "650px"),
        br(),
        helpText("Bubble size = significant hits, colour = -log10(combined p-value).")
      )
    )
  ),
  
  # 5.3 Chord diagram tab --------------------------------------------
  tabPanel(
    "Chord diagram",
    sidebarLayout(
      sidebarPanel(
        h4("Chord settings"),
        
        sliderInput(
          "chord_cex_legend",
          "Legend text size:",
          min   = 0.5,
          max   = 1.4,
          step  = 0.1,
          value = 0.8
        ),
        
        selectInput(
          "chord_palette_mode",
          "Color palette:",
          choices  = c("Bright (per pathway)", "By category"),
          selected = "Bright (per pathway)"
        ),
        
        selectInput(
          "chord_legend_order",
          "Legend order:",
          choices  = c("By ID", "By pathway name", "By category"),
          selected = "By ID"
        ),
        
        downloadButton("download_chord", "Download PNG")
      ),
      mainPanel(
        plotOutput("chord_plot", height = "650px"),
        br(),
        helpText("Chord diagram showing the main 10 KEGG pathways (P1–P10) and regulation (Up / Down in stress).")
      )
    )
  )
)

############################
# 6. SERVER
############################

server <- function(input, output, session) {
  
  # 6.1 Volcano -------------------------------------------------------
  volcano_filtered <- reactive({
    p_cut  <- input$vol_p_cut
    fc_cut <- input$vol_fc_cut
    
    vol_df %>%
      mutate(
        regulation = case_when(
          min_adjP <= p_cut & mean_logFC >=  fc_cut ~ "Up in stress",
          min_adjP <= p_cut & mean_logFC <= -fc_cut ~ "Down in stress",
          TRUE                                      ~ "Not significant"
        ),
        regulation = factor(
          regulation,
          levels = c("Down in stress", "Up in stress", "Not significant")
        ),
        neg_log10_p = -log10(min_adjP)
      )
  })
  
  output$volcano_plot <- renderPlot({
    df <- volcano_filtered()
    
    ggplot(
      df,
      aes(x = mean_logFC, y = neg_log10_p, colour = regulation)
    ) +
      geom_point(size = 2.8, alpha = input$vol_alpha) +
      scale_colour_manual(values = vol_cols, name = "Regulation") +
      geom_vline(
        xintercept = c(-input$vol_fc_cut, input$vol_fc_cut),
        linetype   = "dashed",
        colour     = "grey60",
        linewidth  = 0.5
      ) +
      geom_hline(
        yintercept = -log10(input$vol_p_cut),
        linetype   = "dashed",
        colour     = "grey60",
        linewidth  = 0.5
      ) +
      labs(
        title = "Volcano plot — NEG core metabolites",
        x     = "log2 fold change (stress / control)",
        y     = "-log10(min adjusted p-value)"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title       = element_text(face = "bold", hjust = 0.5),
        legend.position  = "right",
        legend.title     = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
  })
  
  output$download_volcano <- downloadHandler(
    filename = function() {
      paste0(
        "neg_core_volcano_p", input$vol_p_cut,
        "_fc", input$vol_fc_cut, ".png"
      )
    },
    content = function(file) {
      df <- volcano_filtered()
      png(file, width = 2000, height = 1600, res = 300)
      print(
        ggplot(
          df,
          aes(x = mean_logFC, y = neg_log10_p, colour = regulation)
        ) +
          geom_point(size = 2.8, alpha = 0.9) +
          scale_colour_manual(values = vol_cols, name = "Regulation") +
          geom_vline(
            xintercept = c(-input$vol_fc_cut, input$vol_fc_cut),
            linetype   = "dashed",
            colour     = "grey60",
            linewidth  = 0.5
          ) +
          geom_hline(
            yintercept = -log10(input$vol_p_cut),
            linetype   = "dashed",
            colour     = "grey60",
            linewidth  = 0.5
          ) +
          labs(
            title = "Volcano plot — NEG core metabolites",
            x     = "log2 fold change (stress / control)",
            y     = "-log10(min adjusted p-value)"
          ) +
          theme_minimal(base_size = 14) +
          theme(
            plot.title       = element_text(face = "bold", hjust = 0.5),
            legend.position  = "right",
            legend.title     = element_text(face = "bold"),
            panel.grid.minor = element_blank()
          )
      )
      dev.off()
    }
  )
  
  # 6.2 KEGG bubble ---------------------------------------------------
  kegg_filtered <- reactive({
    n_top <- input$kegg_top_n
    
    enrich %>%
      filter(
        !is.na(rich_factor),
        !is.na(neg_log10_p),
        sig_hits > 0
      ) %>%
      arrange(neg_log10_p) %>%
      slice_tail(n = n_top) %>%
      mutate(
        pathway_clean = factor(pathway_clean, levels = pathway_clean)
      )
  })
  
  output$kegg_plot <- renderPlot({
    df <- kegg_filtered()
    
    ggplot(
      df,
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
        name   = expression(-log[10](combined~p))
      ) +
      scale_size(range = c(3, 11), name = "Significant hits") +
      labs(
        title = "KEGG pathway enrichment (NEG)",
        x     = "Rich factor (Sig_Hits / Total_Size)",
        y     = NULL
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title       = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text.y      = element_text(size = 11),
        panel.grid.minor = element_blank(),
        legend.position  = "right"
      )
  })
  
  output$download_kegg <- downloadHandler(
    filename = function() {
      paste0("neg_kegg_richfactor_top", input$kegg_top_n, ".png")
    },
    content = function(file) {
      df <- kegg_filtered()
      png(file, width = 2000, height = 1800, res = 300)
      print(
        ggplot(
          df,
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
            name   = expression(-log[10](combined~p))
          ) +
          scale_size(range = c(3, 11), name = "Significant hits") +
          labs(
            x = "Rich factor (Sig_Hits / Total_Size)",
            y = NULL
          ) +
          theme_minimal(base_size = 13) +
          theme(
            plot.title       = element_blank(),
            axis.text.y      = element_text(size = 11),
            panel.grid.minor = element_blank(),
            legend.position  = "right"
          )
      )
      dev.off()
    }
  )
  
  # 6.3 Chord diagram -----------------------------------------------
  output$chord_plot <- renderPlot({
    plot_chord_main(
      cex_legend   = input$chord_cex_legend,
      palette_mode = input$chord_palette_mode,
      legend_order = input$chord_legend_order
    )
  })
  
  output$download_chord <- downloadHandler(
    filename = function() {
      paste0(
        "chord_main_pathways_",
        gsub(" ", "_", input$chord_palette_mode), "_",
        gsub(" ", "_", input$chord_legend_order),
        "_legend-", input$chord_cex_legend, ".png"
      )
    },
    content = function(file) {
      png(file, width = 2800, height = 1600, res = 300)
      plot_chord_main(
        cex_legend   = input$chord_cex_legend,
        palette_mode = input$chord_palette_mode,
        legend_order = input$chord_legend_order
      )
      dev.off()
    }
  )
}

############################
# 7. RUN APP
############################

shinyApp(ui, server)
