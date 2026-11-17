#!/usr/bin/env Rscript

# Reproducible development plotting for one Checkpoint 11 figure directory.
suppressPackageStartupMessages({
  library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 26_plot_checkpoint11_nomination.R <figure_directory>")
d <- normalizePath(args[[1]], mustWork = TRUE)
name <- basename(d)
src <- file.path(d, "source_data.tsv")
if (!file.exists(src) || file.info(src)$size == 0) stop("Missing source_data.tsv")
x <- read.delim(src, check.names = FALSE, stringsAsFactors = FALSE)
out <- file.path(d, paste0(name, ".png"))
theme_pub <- theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold"))

if (name == "Figure_Nomination_01_Funnel") {
  x$stage <- factor(x$stage, levels = rev(x$stage))
  p <- ggplot(x, aes(n, stage)) +
    geom_segment(aes(x = 0, xend = n, yend = stage), linewidth = 7, color = "#7C3AED", lineend = "round") +
    geom_text(aes(label = format(n, big.mark = ",")), hjust = -0.15, size = 3.5) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = "Number of genes", y = NULL, title = "Transparent candidate nomination workflow") + theme_pub
  ggsave(out, p, width = 8.2, height = 4.8, dpi = 300)
} else if (name == "Figure_Nomination_02_LOPO_Stability") {
  lfc_cols <- grep("^LOPO_P(08|09|34|36|38)_LFC$", names(x), value = TRUE)
  ord <- x$gene[order(-x$LOPO_negative_count, x$LOPO_median_LFC, x$gene)]
  show <- head(ord, 20)
  z <- do.call(rbind, lapply(lfc_cols, function(v) data.frame(
    gene = x$gene, omitted_patient = gsub("^LOPO_|_LFC$", "", v), LFC = x[[v]], stringsAsFactors = FALSE)))
  z <- z[z$gene %in% show, ]
  z$gene <- factor(z$gene, levels = rev(show))
  z$omitted_patient <- factor(z$omitted_patient, levels = c("P08", "P09", "P34", "P36", "P38"))
  lim <- max(abs(z$LFC), na.rm = TRUE)
  p <- ggplot(z, aes(omitted_patient, gene, fill = LFC)) + geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-lim, lim), name = "Shrunken\nlog2FC") +
    labs(x = "Omitted patient", y = NULL, title = "LOPO stability: top 20 stable Edge-down candidates") + theme_pub
  ggsave(out, p, width = 7.2, height = 8.0, dpi = 300)
} else if (name == "Figure_Nomination_03_LeadingEdge_CoreMap") {
  annotation <- c("gene", "pareto_front", "full_shrunken_LFC", "LOPO_negative_count", "patient_consistency")
  module_cols <- setdiff(names(x), annotation)
  z <- do.call(rbind, lapply(module_cols, function(v) data.frame(gene = x$gene, module = v, pathway_count = x[[v]])))
  z <- z[z$pathway_count > 0, , drop = FALSE]
  gene_order <- x$gene[order(x$pareto_front, x$full_shrunken_LFC)]
  z$gene <- factor(z$gene, levels = rev(gene_order))
  p <- ggplot(z, aes(module, gene, fill = pathway_count)) + geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient(low = "#F3E8FF", high = "#5B21B6", name = "Leading-edge\npathway count") +
    labs(x = "Functional module", y = NULL, title = "NES<0 leading-edge core map") +
    theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(out, p, width = 11.5, height = max(5.5, 0.28 * length(gene_order) + 2), dpi = 300)
} else if (name %in% c("Figure_Nomination_04_Pareto_Global", "Figure_Nomination_04_PDE2A_Context")) {
  x$front_group <- ifelse(x$pareto_front == 1, "Pareto Front 1", "Other fronts")
  p <- ggplot(x, aes(full_shrunken_LFC, edge_down_leading_edge_count,
                     size = functional_module_count, color = front_group)) +
    geom_point(alpha = 0.82) +
    scale_color_manual(values = c("Pareto Front 1" = "#7C3AED", "Other fronts" = "#6B7280")) +
    scale_size_continuous(range = c(2.5, 9), breaks = scales::pretty_breaks(4)) +
    labs(x = "Shrunken log2FC (Edge - Core)", y = "NES<0 significant leading-edge memberships",
         size = "Functional\nmodules", color = NULL, title = "Non-dominated candidate context") + theme_pub
  if (name == "Figure_Nomination_04_PDE2A_Context") {
    pde <- x[x$gene == "PDE2A", ]
    p <- p + geom_text(data = pde, aes(label = gene), inherit.aes = TRUE,
                       nudge_x = 0.08, nudge_y = 3, hjust = 0, size = 3.8, fontface = "bold", show.legend = FALSE)
  }
  ggsave(out, p, width = 8.2, height = 5.8, dpi = 300)
} else stop("Unknown Checkpoint 11 figure directory: ", name)
