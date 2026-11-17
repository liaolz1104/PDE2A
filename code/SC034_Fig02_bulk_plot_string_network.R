#!/usr/bin/env Rscript

# Plot frozen Checkpoint 8 source data. Development output is PNG 300 dpi only.

args <- commandArgs(trailingOnly = TRUE)
fig_dir <- if (length(args) >= 1) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
fig_name <- basename(fig_dir)
source_file <- file.path(fig_dir, "source_data.tsv")
if (!file.exists(source_file) || file.info(source_file)$size == 0) stop("Missing source_data.tsv")
suppressPackageStartupMessages({library(ggplot2); library(ggrepel)})
d <- read.delim(source_file, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
colors <- c("Edge-down" = "#0072B2", "Edge-up" = "#D55E00")
theme_set(theme_classic(base_size = 11))

if (grepl("^Figure_PPI_", fig_name)) {
  nodes <- d[d$record_type == "node", , drop = FALSE]
  edges <- d[d$record_type == "edge", , drop = FALSE]
  nodes$degree <- as.numeric(nodes$degree)
  nodes$weighted_degree_strength <- as.numeric(nodes$weighted_degree_strength)
  top <- nodes[order(-nodes$degree, -nodes$weighted_degree_strength, nodes$gene), , drop = FALSE]
  top <- head(top, 6)
  p <- ggplot() +
    geom_segment(data = edges, aes(x, y, xend = xend, yend = yend, alpha = combined_score),
                 color = "#A7ADB4", linewidth = 0.45, show.legend = FALSE) +
    geom_point(data = nodes, aes(x, y, size = degree, color = direction), alpha = 0.90) +
    scale_color_manual(values = colors, name = NULL) +
    scale_size_continuous(range = c(2.2, 10), name = "Degree") +
    coord_equal(clip = "off") +
    theme_void(base_size = 11) +
    theme(legend.position = "right", plot.margin = margin(15, 20, 15, 20))
  if (fig_name == "Figure_PPI_01_84DEG_Global") {
    p <- p + labs(title = "STRING network among 84 FDR-significant DEGs",
                  subtitle = "Input-only functional network; combined score >= 0.400") +
      geom_text_repel(data = top, aes(x, y, label = gene), seed = 20260829,
                      size = 3.2, force = 5, force_pull = 0.15, box.padding = 0.55,
                      point.padding = 0.45, max.time = 5, max.iter = 20000,
                      max.overlaps = Inf, min.segment.length = 0)
  } else if (fig_name == "Figure_PPI_02_PDE2A_Context") {
    pde <- nodes[nodes$gene == "PDE2A", , drop = FALSE]
    p <- p + labs(title = "PDE2A in the frozen 84-DEG STRING network",
                  subtitle = "Identical layout and aesthetics to the global network") +
      geom_text_repel(data = top, aes(x, y, label = gene), seed = 20260829,
                      size = 3.2, force = 5, force_pull = 0.15, box.padding = 0.55,
                      point.padding = 0.45, max.time = 5, max.iter = 20000,
                      max.overlaps = Inf, min.segment.length = 0) +
      geom_point(data = pde, aes(x, y, size = degree), shape = 21, fill = NA,
                 color = "#7C3AED", stroke = 1.2, show.legend = FALSE) +
      geom_text_repel(data = pde, aes(x, y, label = gene), seed = 20260829,
                      color = "#7C3AED", fontface = "bold", size = 3.5,
                      nudge_y = 0.08, force = 5, box.padding = 0.55,
                      max.overlaps = Inf, min.segment.length = 0)
  } else if (fig_name == "Figure_PPI_03_PDE2A_LocalNeighborhood") {
    p <- p + labs(title = "PDE2A first-degree neighborhood within the 84-DEG network",
                  subtitle = "No external STRING interactor added") +
      geom_text_repel(data = nodes, aes(x, y, label = paste0(gene, "\nLFC=", sprintf("%.2f", shrunken_log2FC),
                                                             "; LE=", leading_edge_count)),
                      seed = 20260829, size = 3.0, max.overlaps = Inf, min.segment.length = 0)
  }
  ggsave(file.path(fig_dir, paste0(fig_name, ".png")), p, width = 9.6, height = 7.6, dpi = 300, bg = "white")
} else if (fig_name == "Figure_Prioritization_01_EvidenceMatrix") {
  gene_order <- unique(d$gene[order(d$display_order)])
  field_order <- c("shrunken_log2FC", "patient_consistency", "leading_edge_count",
                   "functional_module_count", "PPI_degree", "enzyme_yes_no")
  labels <- c(shrunken_log2FC = "Shrunken LFC", patient_consistency = "Patient consistency",
              leading_edge_count = "Leading-edge count", functional_module_count = "Module count",
              PPI_degree = "PPI degree", enzyme_yes_no = "Enzyme")
  d$gene <- factor(d$gene, levels = rev(gene_order))
  d$evidence_field <- factor(d$evidence_field, levels = field_order, labels = labels[field_order])
  d$value_label <- ifelse(as.character(d$evidence_field) == "Shrunken LFC", sprintf("%.2f", d$raw_value), sprintf("%.0f", d$raw_value))
  p <- ggplot(d, aes(evidence_field, gene)) +
    geom_tile(aes(fill = within_field_display_scale), color = "white", linewidth = 0.6) +
    geom_text(aes(label = value_label), size = 2.7) +
    scale_fill_gradient(low = "#F7FBFF", high = "#2171B5", guide = "none") +
    labs(title = "Evidence matrix for Edge-down signaling/regulatory candidates",
         subtitle = unique(d$sorting_rule), x = NULL, y = NULL,
         caption = "Color is scaled within each column for display only; no values were summed or weighted.") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          plot.subtitle = element_text(size = 9, color = "#555555"),
          plot.caption = element_text(size = 8, hjust = 0))
  ggsave(file.path(fig_dir, paste0(fig_name, ".png")), p, width = 8.5,
         height = max(7.5, 0.24 * length(gene_order) + 2.2), dpi = 300, bg = "white")
} else {
  stop("Unknown Checkpoint 8 figure directory: ", fig_name)
}
