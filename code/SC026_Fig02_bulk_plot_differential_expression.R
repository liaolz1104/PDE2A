#!/usr/bin/env Rscript

# Recreate one Checkpoint 5 development figure from its source_data.tsv.
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)
fig_dir <- if (length(args) >= 1) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
fig_name <- basename(fig_dir)
source_file <- file.path(fig_dir, "source_data.tsv")
if (!file.exists(source_file) || file.info(source_file)$size == 0) stop("Missing source data: ", source_file)
d <- read.delim(source_file, check.names = FALSE, stringsAsFactors = FALSE)

if (fig_name == "Figure_DEG_01_MA") {
  d$significance <- factor(d$significance, levels = c("NS", "padj < 0.05"))
  p <- ggplot(d, aes(baseMean, log2FoldChange_raw, color = significance)) +
    geom_point(size = 0.75, alpha = 0.55, na.rm = TRUE) +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "#555555") +
    scale_x_log10() +
    scale_color_manual(values = c("NS" = "#B7B7B7", "padj < 0.05" = "#D55E00")) +
    labs(x = "Mean normalized count", y = "log2 fold change (Edge - Core)", color = NULL,
         title = "Paired differential expression: MA plot") +
    theme_classic(base_size = 11) + theme(legend.position = "top")
  ggsave(file.path(fig_dir, "Figure_DEG_01_MA.png"), p, width = 6.2, height = 4.8, dpi = 300, bg = "white")
} else if (fig_name == "Figure_DEG_02_Volcano_global") {
  d$class <- factor(d$class, levels = c("Edge-down", "NS", "Edge-up"))
  p <- ggplot(d, aes(log2FoldChange_shrunk, minus_log10_padj, color = class)) +
    geom_point(size = 0.8, alpha = 0.6, na.rm = TRUE) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.35, color = "#777777") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.35, color = "#777777") +
    geom_text_repel(data = d[d$auto_label, , drop = FALSE], aes(label = display_label),
                    size = 2.7, box.padding = 0.25, point.padding = 0.15,
                    max.overlaps = Inf, min.segment.length = 0, show.legend = FALSE) +
    scale_color_manual(values = c("Edge-down" = "#0072B2", "NS" = "#B7B7B7", "Edge-up" = "#D55E00")) +
    labs(x = "Shrunken log2 fold change (Edge - Core)", y = expression(-log[10](adjusted~italic(P))),
         color = NULL, title = "Paired differential expression: global volcano") +
    theme_classic(base_size = 11) + theme(legend.position = "top")
  ggsave(file.path(fig_dir, "Figure_DEG_02_Volcano_global.png"), p, width = 6.6, height = 5.1, dpi = 300, bg = "white")
} else if (fig_name == "Figure_DEG_03_TopDEG_heatmap") {
  required <- c("gene_label", "sample_id", "row_z", "patient_id", "region_label")
  if (!all(required %in% names(d))) stop("Heatmap source data lacks required columns")
  gene_order <- unique(d$gene_label)
  sample_order <- unique(d$sample_id)
  mat <- matrix(NA_real_, nrow = length(gene_order), ncol = length(sample_order), dimnames = list(gene_order, sample_order))
  mat[cbind(match(d$gene_label, gene_order), match(d$sample_id, sample_order))] <- d$row_z
  ann <- unique(d[, c("sample_id", "patient_id", "region_label")])
  ann <- ann[match(sample_order, ann$sample_id), , drop = FALSE]
  rownames(ann) <- ann$sample_id
  ann$sample_id <- NULL
  names(ann) <- c("Patient", "Region")
  ann$Patient <- factor(ann$Patient)
  ann$Region <- factor(ann$Region, levels = c("Core", "Edge"))
  ann_colors <- list(
    Region = c(Core = "#0072B2", Edge = "#D55E00"),
    Patient = c(P08 = "#4E79A7", P09 = "#59A14F", P34 = "#F28E2B", P36 = "#B07AA1", P38 = "#E15759")
  )
  png(file.path(fig_dir, "Figure_DEG_03_TopDEG_heatmap.png"), width = 2400, height = 2700, res = 300)
  pheatmap(mat, scale = "none", cluster_rows = TRUE, cluster_cols = TRUE,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
           border_color = NA, annotation_col = ann, annotation_colors = ann_colors,
           fontsize = 8, fontsize_row = 7.2, fontsize_col = 8.5,
           main = unique(d$selection_rule), angle_col = 45)
  dev.off()
} else {
  stop("Unrecognized figure directory: ", fig_name)
}
