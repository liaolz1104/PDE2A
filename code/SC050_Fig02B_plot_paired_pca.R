#!/usr/bin/env Rscript
# Redraw frozen paired PCA and expression; do not refit or reposition data.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
root <- sub("/+$", "", args[[1]])
source(file.path(root, "scripts/30_manuscript_redraw_helpers.R"))
font_dir <- if (.Platform$OS.type == "windows") "C:/Windows/Fonts" else
  file.path(root, "assets/fonts/Arial")
fonts <- file.path(font_dir, c("arial.ttf", "arialbd.ttf", "ariali.ttf", "arialbi.ttf"))
stopifnot(all(file.exists(fonts)))
if (!"Arial" %in% unique(systemfonts::system_fonts()$family)) {
  systemfonts::register_font("Arial", plain = fonts[1], bold = fonts[2],
                            italic = fonts[3], bolditalic = fonts[4])
}
out <- file.path(root, "results/07_figures/manuscript_redraw_v4_Arial7pt")
pca <- read_tsv(file.path(out, "Figure_02_PairedPCA_CoreEdge/source_data.tsv"))
expr <- read_tsv(file.path(out, "Figure_03_PDE2A_PairedExpression/source_data.tsv"))
stopifnot(nrow(pca) == 10L, nrow(expr) == 10L,
          !anyDuplicated(pca$sample_id), !anyDuplicated(expr$sample_id))
for (dat in list(pca, expr)) {
  stopifnot(all(table(dat$patient_id, dat$region) == 1L))
}
style <- theme_manuscript(FONT_PT) + theme(
  text = element_text(family = "Arial"),
  plot.title = element_text(size = FONT_PT, face = "bold", margin = margin(b = 5)),
  plot.subtitle = element_text(size = FONT_PT, colour = "#555555", margin = margin(b = 10)),
  axis.text = element_text(size = FONT_PT), axis.title = element_text(size = FONT_PT),
  plot.margin = margin(9, 10, 9, 9))
colours <- c(Core = COL_CORE, Edge = COL_EDGE)
pca$region_label <- factor(pca$region_label, levels = c("Core", "Edge"))
pca$display_label <- paste(pca$patient_id, pca$region_label)
p <- ggplot(pca, aes(PC1, PC2)) +
  geom_line(aes(group = patient_id), colour = "#B2B2B2", linewidth = 0.45) +
  geom_point(aes(fill = region_label), shape = 21, colour = "white", stroke = 0.35, size = 3.1) +
  geom_text_repel(aes(label = display_label), family = "Arial", size = FONT_MM,
    colour = "#333333", seed = 20260909, box.padding = 0.6, point.padding = 0.5,
    min.segment.length = 0, segment.colour = "#888888", segment.size = 0.25,
    max.overlaps = Inf) +
  scale_fill_manual(values = colours, name = NULL) +
  scale_x_continuous(breaks = seq(-100, 100, 50), expand = expansion(mult = 0.14)) +
  scale_y_continuous(breaks = seq(-100, 150, 50), expand = expansion(mult = 0.12)) +
  labs(title = "Paired Core-Edge PCA", subtitle = "5 patients; 10 paired samples",
    x = sprintf("PC1 (%.1f%%)", pca$PC1_variance_percent[1]),
    y = sprintf("PC2 (%.1f%%)", pca$PC2_variance_percent[1])) + style +
  theme(legend.position = "top", legend.justification = "left",
        legend.margin = margin(0, 0, 3, 0))

# Endpoint labels have separate display positions; all measured points stay fixed.
expr$region_label <- factor(expr$region_label, levels = c("Core", "Edge"))
ends <- expr[expr$region == "ti", ]
ends <- ends[order(ends$log2_normalized_count_plus1), ]
ends$label_y <- ends$log2_normalized_count_plus1
for (i in 2:nrow(ends)) ends$label_y[i] <- max(ends$label_y[i], ends$label_y[i-1] + 0.27)
expr$label_y <- ends$label_y[match(expr$patient_id, ends$patient_id)]
expr$label_y[expr$region == "c"] <- NA_real_
q <- ggplot(expr, aes(x, log2_normalized_count_plus1)) +
  geom_line(aes(group = patient_id), colour = "#999999", linewidth = 0.55) +
  geom_point(aes(fill = region_label), shape = 21, colour = "white", stroke = 0.4, size = 3.3) +
  geom_segment(data = ends, aes(x = 2.035, xend = 2.16,
      y = log2_normalized_count_plus1, yend = label_y),
    colour = "#888888", linewidth = 0.28) +
  geom_text(data = ends, aes(x = 2.19, y = label_y, label = patient_id),
    hjust = 0, family = "Arial", size = FONT_MM, colour = "#333333") +
  scale_fill_manual(values = colours, guide = "none") +
  scale_x_continuous(breaks = c(1, 2), labels = c("Core", "Edge")) +
  scale_y_continuous(breaks = seq(5.5, 9.5, 1)) +
  coord_cartesian(xlim = c(0.8, 2.45), ylim = c(5.3, 9.65), expand = FALSE) +
  labs(title = "PDE2A expression",
    subtitle = sprintf("Paired DESeq2; shrunken log2FC = %.2f; FDR = %.3f\nEdge < Core in %d/5 patients",
      expr$DESeq2_shrunken_log2FC[1], expr$DESeq2_padj[1], expr$patient_direction_consistency[1]),
    x = NULL, y = "log2(normalized count + 1)") + style
script <- file.path(root, "scripts/42_polish_paired_figures_arial.R")
save_bundle(out, "Figure_02_PairedPCA_Arial_Refined", p, 160, 135, list(data = pca), script)
save_bundle(out, "Figure_03_PDE2A_Paired_Arial_Refined", q, 155, 140, list(data = expr), script)
for (name in c("Figure_02_PairedPCA_Arial_Refined", "Figure_03_PDE2A_Paired_Arial_Refined")) {
  write_tsv(data.frame(family = "Arial", file = fonts, md5 = unname(tools::md5sum(fonts))),
            file.path(out, name, "FONT_AUDIT.tsv"))
}
cat("Completed: two figures; 10 original points each; Arial TrueType.\n")
