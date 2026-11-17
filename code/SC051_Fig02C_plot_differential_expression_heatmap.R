#!/usr/bin/env Rscript

# Presentation-only redraw from frozen heatmap values and frozen dendrogram
# segments. No VST transformation, clustering, or DESeq2 model is recomputed.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
root <- sub("/+$", "", args[[1]])
source(file.path(root, "scripts", "30_manuscript_redraw_helpers.R"))
suppressPackageStartupMessages(library(patchwork))

input_dir <- file.path(root, "results", "07_figures", "manuscript_redraw_v3_single_panel",
                       "Figure_06_CoreEdge_SplitDendrogramHeatmap")
out_root <- file.path(root, "results", "07_figures", "manuscript_redraw_v4_Arial7pt")
figure_name <- "Figure_06_CoreEdge_SplitDendrogramHeatmap"
heat <- read_tsv(file.path(input_dir, "source_data.tsv"))
branches <- read_tsv(file.path(input_dir, "dendrogram_segments.tsv"))

heat$cluster_order <- as.integer(heat$cluster_order)
gene_key <- unique(heat[, c("cluster_order", "gene_label", "direction_block")])
gene_key <- gene_key[order(gene_key$cluster_order), ]
sample_key <- unique(heat[, c("sample_id", "sample_label", "patient_id", "region")])
sample_order <- c("P08_C", "P09_C", "P34_C", "P36_C", "P38_C",
                  "P08_TI", "P09_TI", "P34_TI", "P36_TI", "P38_TI")
sample_key <- sample_key[match(sample_order, sample_key$sample_id), ]
stopifnot(!anyNA(sample_key$sample_id), nrow(gene_key) == 63L)
heat$x <- match(heat$sample_id, sample_order)

gene_labels <- setNames(gene_key$gene_label, gene_key$cluster_order)

p_tree <- ggplot(branches) +
  geom_segment(aes(x = x0, xend = x1, y = y0 + 0.5, yend = y1 + 0.5),
               linewidth = 0.28, colour = "#444444", lineend = "square") +
  scale_x_reverse(expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_reverse(limits = c(nrow(gene_key) + 0.5, 0.5), expand = c(0, 0)) +
  theme_void(base_family = "Arial", base_size = FONT_PT) +
  theme(text = element_text(family = "Arial", size = FONT_PT),
        plot.margin = margin(5, 0, 16, 5))

p_heat <- ggplot(heat, aes(x = x, y = cluster_order, fill = row_z)) +
  geom_tile(width = 0.98, height = 0.98) +
  geom_vline(xintercept = 5.5, linewidth = 0.55, colour = "white") +
  geom_hline(yintercept = 12.5, linewidth = 0.55, colour = "white") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-max(abs(heat$row_z)), max(abs(heat$row_z))),
                       name = "Row z-score") +
  scale_x_continuous(breaks = seq_along(sample_order), labels = sample_key$sample_label,
                     expand = c(0, 0), position = "bottom") +
  scale_y_reverse(breaks = gene_key$cluster_order, labels = gene_labels,
                  limits = c(nrow(gene_key) + 0.5, 0.5), expand = c(0, 0),
                  position = "right") +
  labs(title = "Core-Edge differential expression", x = NULL, y = NULL) +
  theme_minimal(base_family = "Arial", base_size = FONT_PT) +
  theme(text = element_text(family = "Arial", size = FONT_PT),
        plot.title = element_text(family = "Arial", face = "bold", size = FONT_PT,
                                  hjust = 0, margin = margin(b = 4)),
        axis.text.x = element_text(family = "Arial", size = FONT_PT, angle = 45,
                                   hjust = 1, vjust = 1),
        axis.text.y = element_text(family = "Arial", size = FONT_PT,
                                   colour = "#222222"),
        panel.grid = element_blank(), axis.ticks = element_blank(),
        legend.title = element_text(family = "Arial", size = FONT_PT, face = "bold"),
        legend.text = element_text(family = "Arial", size = FONT_PT),
        legend.position = "right", plot.margin = margin(5, 7, 5, 0))

plot_obj <- p_tree + p_heat + plot_layout(widths = c(1.1, 8.9))
save_bundle(out_root, figure_name, plot_obj, 180, 250,
            list(data = heat), file.path(root, "scripts", "45_redraw_split_heatmap_arial7_from_frozen.R"))
write_tsv(branches, file.path(out_root, figure_name, "dendrogram_segments.tsv"))
cat("Frozen split heatmap redrawn in Arial 7 pt.\n")
