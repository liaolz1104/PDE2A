# Shared visual and export helpers for the PDE2A manuscript redraw.
# This file changes presentation only; it does not recalculate statistics.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(grid)
  library(ragg)
})

COL_CORE <- "#3C5488"
COL_EDGE <- "#E64B35"
COL_PARA <- "#00A087"
COL_NEUTRAL <- "#B8B8B8"
COL_NEUTRAL_DARK <- "#5F6368"
COL_PDE2A <- "#7C3AED"
COL_LIGHT_BLUE <- "#DDE7F2"
COL_DARK_BLUE <- "#1F4E79"
COL_LIGHT_RED <- "#F7D9D4"

# Submission typography: all figure text is rendered as 7-point Arial.
# ggplot2 text geoms use millimetres, while theme text uses points.
FONT_PT <- 7
FONT_MM <- FONT_PT / ggplot2::.pt

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing frozen input: ", path)
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
             quote = "", comment.char = "", na.strings = c("NA", ""))
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  # Keep TSV records structurally valid when wrapped plot labels contain line breaks.
  factor_cols <- vapply(x, is.factor, logical(1))
  x[factor_cols] <- lapply(x[factor_cols], as.character)
  character_cols <- vapply(x, is.character, logical(1))
  x[character_cols] <- lapply(x[character_cols], function(z) {
    gsub("[\r\n]+", "\\\\n", z)
  })
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

wrap_text <- function(x, width = 32L) {
  vapply(x, function(z) paste(strwrap(as.character(z), width = width), collapse = "\n"), character(1))
}

clean_term <- function(x) {
  y <- gsub("\\s*\\([A-Za-z_]+:[0-9]+\\)\\s*$", "", x)
  y <- gsub("\\s*R-HSA-[0-9]+\\s*$", "", y)
  y <- gsub("_", " ", y, fixed = TRUE)
  trimws(y)
}

# Standardized display terminology. Original database terms remain unchanged in
# the source-data `pathway` column for provenance.
standardize_pathway_name <- function(x, collection = NA_character_) {
  y <- clean_term(x)
  replacements <- c(
    "E2F Targets" = "E2F Targets",
    "G2-M Checkpoint" = "G2/M Checkpoint",
    "Interferon Alpha Response" = "Interferon-alpha Response",
    "Interferon Gamma Response" = "Interferon-gamma Response",
    "mTORC1 Signaling" = "mTORC1 Signaling",
    "Myc Targets V1" = "MYC Targets V1",
    "Estrogen Response Early" = "Estrogen Response Early",
    "KRAS Signaling Dn" = "KRAS Signaling DN",
    "Cytoplasmic Translation" = "Cytoplasmic translation",
    "Myofibril Assembly" = "Myofibril assembly",
    "Muscle Contraction" = "Muscle contraction",
    "Peptide Biosynthetic Process" = "Peptide biosynthetic process",
    "Sarcomere Organization" = "Sarcomere organization",
    "Metal Ion Transport" = "Metal ion transport",
    "DNA Metabolic Process" = "DNA metabolic process",
    "DNA Repair" = "DNA repair",
    "DNA-templated DNA Replication" = "DNA-templated DNA replication",
    "DNA Replication" = "DNA replication",
    "Double-Strand Break Repair Via Homologous Recombination" = "Double-strand break repair via homologous recombination",
    "Double-Strand Break Repair" = "Double-strand break repair",
    "Eukaryotic Translation Elongation" = "Eukaryotic translation elongation",
    "Formation Of A Pool Of Free 40S Subunits" = "Formation of a pool of free 40S subunits",
    "Eukaryotic Translation Termination" = "Eukaryotic translation termination",
    "Peptide Chain Elongation" = "Peptide chain elongation",
    "Viral mRNA Translation" = "Viral mRNA translation",
    "L13a-mediated Translational Silencing Of Ceruloplasmin Expression" = "L13a-mediated translational silencing of ceruloplasmin expression",
    "Cell Cycle, Mitotic" = "Mitotic cell cycle",
    "Cell Cycle Checkpoints" = "Cell cycle checkpoints",
    "M Phase" = "M phase",
    "Mitotic Prometaphase" = "Mitotic prometaphase",
    "Mitotic Metaphase And Anaphase" = "Mitotic metaphase and anaphase",
    "Mitotic Anaphase" = "Mitotic anaphase",
    "Regulation Of Angiogenesis" = "Regulation of angiogenesis",
    "Blood Vessel Endothelial Cell Migration" = "Blood vessel endothelial cell migration",
    "Epithelium Development" = "Epithelium development",
    "Cellular Response To Transforming Growth Factor Beta Stimulus" = "Cellular response to transforming growth factor beta stimulus",
    "Cellular Response To Growth Factor Stimulus" = "Cellular response to growth factor stimulus",
    "Positive Regulation Of Angiogenesis" = "Positive regulation of angiogenesis",
    "External Encapsulating Structure Organization" = "External encapsulating structure organization",
    "Extracellular Structure Organization" = "Extracellular structure organization",
    "Extracellular Matrix Organization" = "Extracellular matrix organization",
    "Activation Of Matrix Metalloproteinases" = "Activation of matrix metalloproteinases",
    "Collagen Degradation" = "Collagen degradation",
    "Extracellular Matrix Disassembly" = "Extracellular matrix disassembly",
    "Cellular Component Disassembly" = "Cellular component disassembly",
    "Degradation Of Extracellular Matrix" = "Degradation of extracellular matrix",
    "CaM Pathway" = "CaM pathway"
  )
  hit <- match(y, names(replacements))
  y[!is.na(hit)] <- unname(replacements[hit[!is.na(hit)]])
  y
}

standardize_functional_module <- function(x) {
  replacements <- c(
    "Other / unresolved" = "Other or unresolved",
    "Vascular / endothelial" = "Vascular and endothelial programs",
    "Signaling / regulatory" = "Signaling and regulatory processes",
    "ECM / stromal remodeling" = "Extracellular matrix and stromal remodeling",
    "GPCR / adenylate cyclase signaling" = "GPCR and adenylate cyclase signaling",
    "Epithelial differentiation" = "Epithelial differentiation",
    "Immune / interferon" = "Immune and interferon responses",
    "Metabolism / endocrine" = "Metabolic and endocrine processes",
    "NO / cGMP signaling" = "Nitric oxide and cGMP signaling",
    "Growth-factor signaling" = "Growth factor signaling",
    "Cyclic nucleotide signaling" = "Cyclic nucleotide signaling",
    "Hemostasis / platelet" = "Hemostasis and platelet function",
    "Translation / ribosome" = "Translation and ribosome biogenesis",
    "Cell cycle / proliferation" = "Cell cycle and proliferation",
    "DNA replication / repair" = "DNA replication and repair",
    "Stress / hypoxia" = "Stress and hypoxia responses"
  )
  y <- as.character(x)
  hit <- match(y, names(replacements))
  y[!is.na(hit)] <- unname(replacements[hit[!is.na(hit)]])
  y
}

theme_manuscript <- function(base_size = FONT_PT) {
  base_size <- FONT_PT
  theme_classic(base_family = "Arial", base_size = base_size) +
    theme(
      text = element_text(family = "Arial", size = FONT_PT),
      plot.title = element_text(face = "bold", size = FONT_PT, hjust = 0),
      plot.subtitle = element_text(size = FONT_PT, colour = COL_NEUTRAL_DARK, hjust = 0,
                                   margin = margin(b = 3)),
      axis.title = element_text(size = FONT_PT),
      axis.text = element_text(size = FONT_PT, colour = "#222222"),
      axis.line = element_line(linewidth = 0.35, colour = "#222222"),
      axis.ticks = element_line(linewidth = 0.3, colour = "#222222"),
      axis.ticks.length = unit(1.5, "mm"),
      legend.title = element_text(size = FONT_PT, face = "bold"),
      legend.text = element_text(size = FONT_PT),
      legend.key.height = unit(3.5, "mm"),
      legend.key.width = unit(4, "mm"),
      plot.margin = margin(5, 6, 5, 5),
      plot.tag = element_text(face = "bold", size = FONT_PT),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = FONT_PT)
    )
}

theme_heatmap <- function(base_size = FONT_PT) {
  base_size <- FONT_PT
  theme_minimal(base_family = "Arial", base_size = base_size) +
    theme(
      text = element_text(family = "Arial", size = FONT_PT),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(colour = "#222222", size = FONT_PT),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold", size = FONT_PT, hjust = 0),
      plot.subtitle = element_text(size = FONT_PT, colour = COL_NEUTRAL_DARK, hjust = 0,
                                   margin = margin(b = 3)),
      legend.title = element_text(size = FONT_PT, face = "bold"),
      legend.text = element_text(size = FONT_PT),
      plot.margin = margin(5, 6, 5, 5),
      plot.tag = element_text(face = "bold", size = FONT_PT)
    )
}

save_png <- function(plot, path, width_mm, height_mm) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(filename = path, width = width_mm / 25.4, height = height_mm / 25.4,
                 units = "in", res = 300, type = "cairo", bg = "white")
  print(plot)
  invisible(dev.off())
  if (!file.exists(path) || file.info(path)$size <= 0) stop("Figure export failed: ", path)
}

save_bundle <- function(output_root, figure_name, plot, width_mm, height_mm,
                        sources, plotting_script) {
  if (identical(Sys.getenv("PDE2A_SKIP_BUNDLE_EXPORT", unset = "0"), "1")) {
    return(invisible(NULL))
  }
  fig_dir <- file.path(output_root, figure_name)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  if (length(sources) == 1L) {
    write_tsv(sources[[1]], file.path(fig_dir, "source_data.tsv"))
  } else {
    source_index <- data.frame(
      panel = names(sources),
      source_file = paste0("source_data_", names(sources), ".tsv"),
      rows = vapply(sources, nrow, integer(1)),
      stringsAsFactors = FALSE
    )
    for (nm in names(sources)) {
      write_tsv(sources[[nm]], file.path(fig_dir, paste0("source_data_", nm, ".tsv")))
    }
    write_tsv(source_index, file.path(fig_dir, "source_data.tsv"))
  }
  save_png(plot, file.path(fig_dir, paste0(figure_name, ".png")), width_mm, height_mm)
  file.copy(plotting_script, file.path(fig_dir, "plotting_script.R"), overwrite = TRUE)
  invisible(fig_dir)
}

region_label <- function(x) {
  out <- ifelse(tolower(x) == "c", "Core",
                ifelse(tolower(x) == "ti", "Edge", "Paracancerous"))
  factor(out, levels = c("Core", "Edge", "Paracancerous"))
}

direction_factor <- function(x) {
  factor(x, levels = c("Edge-down", "NS", "Edge-up"))
}

direction_scale <- scale_colour_manual(
  values = c("Edge-down" = COL_CORE, "NS" = COL_NEUTRAL, "Edge-up" = COL_EDGE),
  drop = FALSE
)

fill_direction_scale <- scale_fill_manual(
  values = c("Edge-down" = COL_CORE, "NS" = COL_NEUTRAL, "Edge-up" = COL_EDGE),
  drop = FALSE
)

region_colour_scale <- scale_colour_manual(
  values = c("Core" = COL_CORE, "Edge" = COL_EDGE, "Paracancerous" = COL_PARA),
  drop = FALSE
)

region_fill_scale <- scale_fill_manual(
  values = c("Core" = COL_CORE, "Edge" = COL_EDGE, "Paracancerous" = COL_PARA),
  drop = FALSE
)
