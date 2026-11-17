#!/usr/bin/env Rscript

# Recreate one Checkpoint 6 development figure from source_data.tsv.
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
fig_dir <- if (length(args) >= 1) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
fig_name <- basename(fig_dir)
source_file <- file.path(fig_dir, "source_data.tsv")
if (!file.exists(source_file) || file.info(source_file)$size == 0) stop("Missing source data: ", source_file)
d <- read.delim(source_file, check.names = FALSE, stringsAsFactors = FALSE)

plot_gsea <- function(d, title, output) {
  d$pathway_label <- factor(d$pathway_label, levels = rev(unique(d$pathway_label)))
  d$direction <- factor(d$direction, levels = c("Core-enriched", "Edge-enriched"))
  p <- ggplot(d, aes(NES, pathway_label, size = geneSetSize, color = minus_log10_padj)) +
    geom_vline(xintercept = 0, linewidth = 0.35, color = "#777777") +
    geom_point(alpha = 0.85) +
    scale_color_gradient(low = "#56B4E9", high = "#D55E00") +
    labs(x = "Normalized enrichment score (Edge - Core)", y = NULL,
         size = "Gene-set size", color = expression(-log[10](FDR)), title = title) +
    theme_classic(base_size = 10.5) + theme(legend.position = "right")
  ggsave(file.path(fig_dir, output), p, width = 8.0, height = max(4.5, 0.28 * nrow(d) + 1.6), dpi = 300, bg = "white")
}

plot_ora <- function(d, title, output) {
  d$pathway_label <- factor(d$pathway_label, levels = rev(unique(d$pathway_label)))
  p <- ggplot(d, aes(foldEnrichment, pathway_label, size = overlapSize, color = minus_log10_padj)) +
    geom_point(alpha = 0.85) +
    facet_grid(collection ~ ., scales = "free_y", space = "free_y") +
    scale_color_gradient(low = "#56B4E9", high = "#D55E00") +
    labs(x = "Fold enrichment", y = NULL, size = "Overlap genes",
         color = expression(-log[10](FDR)), title = title) +
    theme_classic(base_size = 10.5) + theme(strip.background = element_blank(), strip.text.y = element_text(angle = 0))
  ggsave(file.path(fig_dir, output), p, width = 8.2, height = max(5.0, 0.27 * nrow(d) + 1.8), dpi = 300, bg = "white")
}

if (fig_name == "Figure_Enrichment_01_GSEA_Hallmark") {
  plot_gsea(d, "Hallmark preranked GSEA", "Figure_Enrichment_01_GSEA_Hallmark.png")
} else if (fig_name == "Figure_Enrichment_02_GSEA_GO_Reactome") {
  plot_gsea(d, "GO BP and Reactome preranked GSEA", "Figure_Enrichment_02_GSEA_GO_Reactome.png")
} else if (fig_name == "Figure_Enrichment_03_ORA_Edge_down") {
  plot_ora(d, "Edge-down DEG over-representation analysis", "Figure_Enrichment_03_ORA_Edge_down.png")
} else if (fig_name == "Figure_Enrichment_04_ORA_Edge_up") {
  plot_ora(d, "Edge-up DEG over-representation analysis", "Figure_Enrichment_04_ORA_Edge_up.png")
} else if (fig_name %in% c("Figure_Candidate_Context_global", "Figure_Candidate_Context_PDE2A")) {
  d$direction <- factor(d$direction, levels = c("Edge-down", "Edge-up"))
  p <- ggplot(d, aes(shrunken_LFC, minus_log10_padj, color = direction,
                     size = log10_baseMean_plus1, alpha = patient_direction_consistency)) +
    geom_point() +
    geom_vline(xintercept = 0, linewidth = 0.35, color = "#777777") +
    scale_color_manual(values = c("Edge-down" = "#0072B2", "Edge-up" = "#D55E00")) +
    scale_alpha_continuous(range = c(0.35, 0.95), breaks = 3:5) +
    scale_x_continuous(expand = expansion(mult = c(0.10, 0.10))) +
    scale_y_continuous(expand = expansion(mult = c(0.12, 0.08))) +
    labs(x = "Shrunken log2 fold change (Edge - Core)", y = expression(-log[10](adjusted~italic(P))),
         color = NULL, size = "log10(baseMean + 1)", alpha = "Patient consistency",
         title = "Functional candidate context among 84 FDR-significant genes") +
    theme_classic(base_size = 10.5) + theme(legend.position = "right")
  if (fig_name == "Figure_Candidate_Context_PDE2A") {
    p <- p + geom_text_repel(data = d[d$gene == "PDE2A", , drop = FALSE], aes(label = gene),
                             size = 3.2, nudge_y = 0.12, direction = "y",
                             min.segment.length = 0, show.legend = FALSE)
  }
  output <- paste0(fig_name, ".png")
  ggsave(file.path(fig_dir, output), p, width = 7.0, height = 5.2, dpi = 300, bg = "white")
} else {
  stop("Unrecognized figure directory: ", fig_name)
}
