#!/usr/bin/env Rscript

# Plot Checkpoint 7 development figures from their frozen source_data.tsv.
# Outputs PNG at 300 dpi only. The two scatter versions share all aesthetics;
# the PDE2A-context version adds only a text label.

args <- commandArgs(trailingOnly = TRUE)
fig_dir <- if (length(args) >= 1) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
fig_name <- basename(fig_dir)
source_file <- file.path(fig_dir, "source_data.tsv")
if (!file.exists(source_file) || file.info(source_file)$size == 0) stop("Missing source_data.tsv: ", source_file)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})
d <- read.delim(source_file, check.names = FALSE, stringsAsFactors = FALSE)
edge_colors <- c("Edge-down" = "#0072B2", "Edge-up" = "#D55E00")
theme_set(theme_classic(base_size = 11))

if (fig_name == "Figure_Convergence_01_GeneModuleHeatmap") {
  d$gene_name <- factor(d$gene_name, levels = rev(unique(d$gene_name[order(d$selection_rank)])))
  module_totals <- aggregate(pathway_count ~ functional_module, d, sum)
  module_order <- module_totals$functional_module[order(-module_totals$pathway_count, module_totals$functional_module)]
  d$functional_module <- factor(d$functional_module, levels = module_order)
  p <- ggplot(d, aes(functional_module, gene_name, fill = pathway_count)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient(low = "#F7FBFF", high = "#08519C", name = "Significant\npathway count") +
    labs(title = "Gene–functional-module leading-edge convergence",
         subtitle = unique(d$selection_rule), x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1),
          plot.subtitle = element_text(size = 8.5, color = "#555555"),
          legend.position = "right")
  ggsave(file.path(fig_dir, paste0(fig_name, ".png")), p, width = 10.5, height = 10.2, dpi = 300, bg = "white")
} else if (fig_name %in% c("Figure_Convergence_02_Global", "Figure_Convergence_02_PDE2A_context")) {
  d$neg_log10_padj <- -log10(pmax(d$padj, .Machine$double.xmin))
  d$point_size <- d$total_significant_leading_edge_count
  p <- ggplot(d, aes(shrunken_log2FC, neg_log10_padj)) +
    geom_vline(xintercept = 0, color = "#777777", linewidth = 0.45) +
    geom_point(aes(size = point_size, color = direction, alpha = patient_direction_consistency), stroke = 0) +
    scale_color_manual(values = edge_colors, name = NULL) +
    scale_size_continuous(range = c(2.0, 10.0), name = "Significant leading-edge\npathway count") +
    scale_alpha_continuous(range = c(0.45, 0.95), breaks = 3:5, name = "Patient consistency") +
    labs(title = "DEG–pathway convergence among 84 FDR-significant genes",
         x = "Shrunken log2 fold change (Edge - Core)", y = expression(-log[10](adjusted~P))) +
    coord_cartesian(clip = "off") +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.08))) +
    scale_y_continuous(expand = expansion(mult = c(0.18, 0.12))) +
    theme(legend.position = "right", plot.margin = margin(8, 16, 8, 8))
  if (fig_name == "Figure_Convergence_02_PDE2A_context") {
    label <- d[d$gene_name == "PDE2A", , drop = FALSE]
    if (nrow(label) != 1L) stop("PDE2A is not unique in scatter source data")
    p <- p + geom_text_repel(data = label, aes(label = gene_name), color = edge_colors[["Edge-down"]],
                             seed = 20260829, box.padding = 0.35, point.padding = 0.3,
                             nudge_y = 0.16, direction = "y", min.segment.length = 0,
                             max.overlaps = Inf, show.legend = FALSE)
  }
  ggsave(file.path(fig_dir, paste0(fig_name, ".png")), p, width = 7.6, height = 5.6, dpi = 300, bg = "white")
} else {
  stop("Unknown Checkpoint 7 figure directory: ", fig_name)
}
