#!/usr/bin/env Rscript

# Presentation-only redraw from the frozen Checkpoint 11 Venn source tables.
# No DESeq2 object or statistical result is recomputed.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
root <- sub("/+$", "", args[[1]])
source(file.path(root, "scripts", "30_manuscript_redraw_helpers.R"))
suppressPackageStartupMessages(library(ggforce))

input_dir <- file.path(root, "results", "07_figures", "manuscript_redraw_v3_single_panel",
                       "Figure_06b_PatientDirection_Venn")
out_root <- file.path(root, "results", "07_figures", "manuscript_redraw_v4_Arial7pt")
figure_name <- "Figure_06b_PatientDirection_Venn"

counts <- read_tsv(file.path(input_dir, "source_data.tsv"))
membership <- read_tsv(file.path(input_dir, "gene_set_membership.tsv"))
pde <- read_tsv(file.path(input_dir, "PDE2A_location.tsv"))
patients <- c("P08", "P09", "P34", "P36", "P38")

circles <- data.frame(
  set = patients,
  x0 = c(-0.82, 0.82, -0.98, 0.00, 0.98),
  y0 = c(0.66, 0.66, -0.30, -0.73, -0.30),
  r = rep(1.38, 5), stringsAsFactors = FALSE
)
circle_colours <- c(P08 = "#E69F00", P09 = "#56B4E9", P34 = "#009E73",
                    P36 = "#CC79A7", P38 = "#D55E00")
set_labels <- data.frame(
  set = patients,
  x = c(-1.90, 1.90, -1.95, 0.00, 1.95),
  y = c(1.92, 1.92, -1.92, -2.05, -1.92),
  stringsAsFactors = FALSE
)
pde_marker <- data.frame(x = pde$marker_x[[1]], y = pde$marker_y[[1]])

plot_obj <- ggplot() +
  geom_circle(data = circles, aes(x0 = x0, y0 = y0, r = r, fill = set, colour = set),
              alpha = 0.20, linewidth = 0.65, show.legend = FALSE) +
  geom_text(data = counts[counts$shared_gene_count > 0, , drop = FALSE],
            aes(x = x, y = y, label = plot_label), family = "Arial",
            fontface = "bold", size = FONT_MM, colour = "#222222") +
  geom_point(data = pde_marker, aes(x = x, y = y), inherit.aes = FALSE,
             shape = 21, fill = NA, colour = COL_PDE2A, size = 4.0,
             stroke = 1.15, show.legend = FALSE) +
  ggrepel::geom_label_repel(
    data = pde_marker, aes(x = x, y = y, label = "PDE2A"), inherit.aes = FALSE,
    colour = COL_PDE2A, fill = "white", family = "Arial", fontface = "bold",
    label.size = 0.25, size = FONT_MM, box.padding = 0.35, point.padding = 0.20,
    min.segment.length = 0, segment.colour = COL_PDE2A, segment.size = 0.45,
    force = 1.2, show.legend = FALSE) +
  geom_text(data = set_labels, aes(x = x, y = y, label = set), family = "Arial",
            fontface = "bold", size = FONT_MM, colour = "#222222") +
  scale_fill_manual(values = circle_colours) +
  scale_colour_manual(values = circle_colours) +
  coord_equal(xlim = c(-2.28, 2.28), ylim = c(-2.20, 2.12), expand = FALSE, clip = "off") +
  labs(title = "Patient-wise direction consistency of Edge-down DEGs",
       subtitle = "66 FDR-significant Edge-down genes; membership = normalized Edge < Core",
       caption = "Direction-consistency subsets of the frozen global Edge-down DEG list; not independent patient-level DE tests.") +
  theme_void(base_family = "Arial", base_size = FONT_PT) +
  theme(text = element_text(family = "Arial", size = FONT_PT),
        plot.title = element_text(face = "bold", size = FONT_PT, hjust = 0),
        plot.subtitle = element_text(size = FONT_PT, colour = COL_NEUTRAL_DARK, hjust = 0,
                                     margin = margin(b = 4)),
        plot.caption = element_text(size = FONT_PT, colour = COL_NEUTRAL_DARK, hjust = 0,
                                    margin = margin(t = 5)),
        plot.margin = margin(8, 8, 8, 8))

save_bundle(out_root, figure_name, plot_obj, 150, 128,
            list(data = counts), file.path(root, "scripts", "44_redraw_venn_arial7_from_frozen.R"))
write_tsv(membership, file.path(out_root, figure_name, "gene_set_membership.tsv"))
write_tsv(pde, file.path(out_root, figure_name, "PDE2A_location.tsv"))
cat("Frozen Venn figure redrawn in Arial 7 pt.\n")
