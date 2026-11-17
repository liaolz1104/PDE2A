#!/usr/bin/env Rscript
# Presentation-only export: every finite point keeps its original coordinates.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
root <- sub("/+$", "", args[[1]])
source(file.path(root, "scripts", "30_manuscript_redraw_helpers.R"))
font_dir <- if (.Platform$OS.type == "windows") "C:/Windows/Fonts" else
  file.path(root, "assets/fonts/Arial")
fonts <- file.path(font_dir, c("arial.ttf", "arialbd.ttf", "ariali.ttf", "arialbi.ttf"))
stopifnot(all(file.exists(fonts)))
if (!"Arial" %in% unique(systemfonts::system_fonts()$family)) {
  systemfonts::register_font("Arial", plain = fonts[1], bold = fonts[2],
                            italic = fonts[3], bolditalic = fonts[4])
}
input_root <- file.path(root, "results/07_figures/manuscript_redraw_v3_single_panel")
out_root <- file.path(root, "results/07_figures/manuscript_redraw_v4_Arial7pt")
d <- read_tsv(file.path(input_root, "Figure_05_Volcano_PDE2A_Polished/source_data.tsv"))
d$plotted <- is.finite(d$log2FoldChange_shrunk) & is.finite(d$minus_log10_padj)
d$point_size_mm <- ifelse(!is.na(d$padj) & d$padj < 0.05, 0.95, 0.25)
d$point_alpha <- ifelse(!is.na(d$padj) & d$padj < 0.05, 0.95, 0.16)
x <- d[d$plotted, ]
x <- x[order(x$minus_log10_padj), ]
pde <- x[!is.na(x$gene_name) & x$gene_name == "PDE2A", ]
stopifnot(nrow(pde) == 1L)
xlim <- range(x$log2FoldChange_shrunk) + c(-0.28, 0.28)
ymax <- max(x$minus_log10_padj)
p <- ggplot(x, aes(log2FoldChange_shrunk, minus_log10_padj)) +
  geom_vline(xintercept = c(-1, 1), colour = "#B8B8B8", linewidth = 0.3, linetype = "22") +
  geom_hline(yintercept = -log10(0.05), colour = "#999999", linewidth = 0.3, linetype = "22") +
  geom_point(aes(colour = minus_log10_padj, size = point_size_mm, alpha = point_alpha),
             shape = 16, stroke = 0) +
  scale_size_identity() + scale_alpha_identity() +
  scale_colour_gradientn(colours = c("#F04E23", "#FDE78D", "#86BBD8", "#2C7FB8"),
    values = scales::rescale(c(0, -log10(0.05), ymax * 0.75, ymax)),
    limits = c(0, ymax), name = "-log10(FDR)",
    guide = guide_colourbar(barheight = unit(26, "mm"), barwidth = unit(3, "mm"))) +
  annotate("segment", x = pde$log2FoldChange_shrunk - 0.05,
    xend = pde$log2FoldChange_shrunk - 0.26,
    y = pde$minus_log10_padj + 0.02, yend = pde$minus_log10_padj + 0.11,
    colour = COL_PDE2A, linewidth = 0.4) +
  annotate("text", x = pde$log2FoldChange_shrunk - 0.30,
    y = pde$minus_log10_padj + 0.12, label = "PDE2A", hjust = 1,
    family = "Arial", fontface = "bold", size = FONT_MM, colour = COL_PDE2A) +
  scale_x_continuous(breaks = seq(-3, 4, 1)) +
  scale_y_continuous(breaks = seq(0, ceiling(ymax), 1)) +
  coord_cartesian(xlim = xlim, ylim = c(-0.035, ymax + 0.18), expand = FALSE) +
  labs(title = "Core-Edge differential expression",
    x = "Shrunken log2 fold change (Edge - Core)", y = "-log10(adjusted P)") +
  theme_manuscript(FONT_PT) +
  theme(text = element_text(family = "Arial"), legend.position = "right",
    plot.title = element_text(size = FONT_PT, face = "bold", margin = margin(b = 10)),
    axis.title = element_text(size = FONT_PT), axis.text = element_text(size = FONT_PT),
    legend.title = element_text(size = FONT_PT), legend.text = element_text(size = FONT_PT),
    plot.margin = margin(8, 8, 8, 8))
name <- "Figure_05_Volcano_Arial_Refined"
save_bundle(out_root, name, p, 180, 130, list(data = d),
            file.path(root, "scripts/41_polish_volcano_arial.R"))
write_tsv(data.frame(family = "Arial", file = fonts, md5 = unname(tools::md5sum(fonts))),
          file.path(out_root, name, "FONT_AUDIT.tsv"))
cat("Rows:", nrow(d), "; plotted:", nrow(x), "; PDE2A:",
    pde$log2FoldChange_shrunk, pde$padj, "\n")
