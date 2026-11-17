#!/usr/bin/env Rscript

# Rebuild the primary self-cohort manuscript figures from frozen source tables.
# No statistical model, threshold, candidate set, or pathway result is recomputed.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 31_plot_manuscript_main_figures.R <project_root>")
root <- sub("/+$", "", args[[1]])
script_path <- file.path(root, "scripts", "31_plot_manuscript_main_figures.R")
source(file.path(root, "scripts", "30_manuscript_redraw_helpers.R"))

output_root <- file.path(root, "results", "07_figures", "manuscript_redraw_v1")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
candidate_root <- file.path(root, "results", "07_figures", "manuscript_candidates")

fig_source <- function(section, figure) {
  file.path(candidate_root, section, figure, "source_data.tsv")
}

# -----------------------------------------------------------------------------
# Figure 1: paired cohort structure and PDE2A regional expression
# -----------------------------------------------------------------------------
pca <- read_tsv(file.path(root, "results", "04_featurecounts", "PCA_main10_paired_coordinates.tsv"))
pca$region_label <- region_label(pca$region)
pca$sample_id <- factor(pca$sample_id, levels = pca$sample_id)

pde <- read_tsv(file.path(root, "results", "05_pde2a_audit", "PDE2A_normalized_counts_preanalysis.tsv"))
pde$region_label <- region_label(pde$region)
pde$x <- ifelse(pde$region_label == "Core", 1, 2)
pde_stats <- read_tsv(file.path(root, "results", "05_deseq2", "PDE2A_GENOMEWIDE_POSITION.tsv"))
stopifnot(nrow(pde_stats) == 1L, pde_stats$gene_name[[1]] == "PDE2A")

patients <- c("P08", "P09", "P34", "P36", "P38")
design <- data.frame(patient_id = factor(patients, levels = rev(patients)), x1 = 1, x2 = 2)
design_points <- rbind(
  data.frame(patient_id = design$patient_id, x = 1, region_label = "Core"),
  data.frame(patient_id = design$patient_id, x = 2, region_label = "Edge")
)
design_points$region_label <- factor(design_points$region_label, levels = c("Core", "Edge"))

p_design <- ggplot(design, aes(y = patient_id)) +
  geom_segment(aes(x = x1, xend = x2, yend = patient_id), linewidth = 0.45, colour = "#B0B0B0") +
  geom_point(data = design_points, aes(x = x, colour = region_label), size = 2.6) +
  region_colour_scale +
  scale_x_continuous(breaks = c(1, 2), labels = c("Core", "Edge"), limits = c(0.8, 2.2)) +
  labs(title = "Paired sampling", x = NULL, y = NULL, colour = "Region") +
  theme_manuscript(8) +
  theme(legend.position = "none", axis.line.y = element_blank(), axis.ticks.y = element_blank())

p_pca <- ggplot(pca, aes(PC1, PC2)) +
  geom_path(aes(group = patient_id), colour = "#B0B0B0", linewidth = 0.45) +
  geom_point(aes(colour = region_label), size = 2.6) +
  geom_text_repel(aes(label = sample_id), size = FONT_MM, family = "Arial", seed = 1101,
                  box.padding = 0.25, point.padding = 0.12, min.segment.length = 0,
                  max.overlaps = Inf, colour = "#222222") +
  region_colour_scale +
  labs(
    title = "Global transcriptomic structure",
    x = sprintf("PC1 (%.1f%%)", unique(pca$PC1_variance_percent)[[1]]),
    y = sprintf("PC2 (%.1f%%)", unique(pca$PC2_variance_percent)[[1]]),
    colour = "Region"
  ) +
  theme_manuscript(8) +
  theme(legend.position = "top")

stat_text <- sprintf(
  "shrunken log2FC = %.2f\nFDR = %.3f\nEdge < Core: %d/5",
  pde_stats$log2FoldChange_shrunk[[1]], pde_stats$padj[[1]],
  pde_stats$direction_consistency_edge_lower[[1]]
)
p_pde <- ggplot(pde, aes(x, log2_normalized_count_plus1, group = patient_id)) +
  geom_line(colour = "#8C8C8C", linewidth = 0.55) +
  geom_point(aes(colour = region_label), size = 2.8) +
  geom_text(data = subset(pde, region_label == "Edge"), aes(label = patient_id),
            x = 2.09, hjust = 0, size = FONT_MM, family = "Arial", colour = "#333333") +
  region_colour_scale +
  scale_x_continuous(breaks = c(1, 2), labels = c("Core", "Edge"), limits = c(0.8, 2.42)) +
  labs(title = "PDE2A expression", subtitle = stat_text,
       x = NULL, y = "log2(normalized count + 1)", colour = "Region") +
  theme_manuscript(8) +
  theme(legend.position = "none")

figure1 <- p_design + p_pca + p_pde +
  plot_layout(widths = c(0.72, 1.28, 1.05), guides = "collect") +
  plot_annotation(
    title = "Paired Core-Edge transcriptomes and PDE2A expression",
    tag_levels = "a",
    theme = theme(plot.title = element_text(family = "Arial", face = "bold", size = FONT_PT),
                  plot.tag = element_text(family = "Arial", face = "bold", size = FONT_PT))
  ) & theme(legend.position = "top")

save_bundle(output_root, "Figure1_Cohort_and_PDE2A", figure1, 183, 82,
            list(cohort_design = design_points, paired_PCA = pca,
                 PDE2A_expression = pde, PDE2A_statistics = pde_stats), script_path)

# -----------------------------------------------------------------------------
# Figure 2: transcriptome-wide differential-expression landscape
# -----------------------------------------------------------------------------
ma <- read_tsv(fig_source("DEG", "Figure_DEG_01_MA"))
ma$direction <- "NS"
ma$direction[!is.na(ma$padj) & ma$padj < 0.05 & ma$log2FoldChange_raw < 0] <- "Edge-down"
ma$direction[!is.na(ma$padj) & ma$padj < 0.05 & ma$log2FoldChange_raw > 0] <- "Edge-up"
ma$direction <- direction_factor(ma$direction)

p_ma <- ggplot(ma, aes(baseMean, log2FoldChange_raw)) +
  geom_point(data = subset(ma, direction == "NS"), aes(colour = direction), alpha = 0.32,
             size = 0.35) +
  geom_point(data = subset(ma, direction != "NS"), aes(colour = direction), alpha = 0.88, size = 0.8) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "#555555") +
  scale_x_log10(labels = label_number()) +
  direction_scale +
  labs(title = "Mean-difference view", x = "Mean normalized count", y = "log2 fold change (Edge - Core)",
       colour = NULL) +
  theme_manuscript(8) +
  theme(legend.position = "top")

volcano <- read_tsv(fig_source("DEG", "Figure_DEG_02_Volcano_global"))
volcano$class <- direction_factor(volcano$class)
volcano$auto_label <- as.logical(volcano$auto_label)
deg_full_for_volcano <- read_tsv(file.path(
  root, "results", "05_deseq2", "DESeq2_Edge_vs_Core_all_genes.tsv"
))
raw_lfc_index <- match(volcano$gene_id, deg_full_for_volcano$gene_id)
volcano$log2FoldChange_raw <- deg_full_for_volcano$log2FoldChange_raw[raw_lfc_index]
if (anyNA(raw_lfc_index)) stop("Failed to map all volcano genes to the frozen DESeq2 table")
volcano$template_class <- ifelse(
  !is.na(volcano$padj) & volcano$padj < 0.05 & abs(volcano$log2FoldChange_shrunk) >= 1,
  "FDR < 0.05 and |LFC| >= 1",
  ifelse(!is.na(volcano$padj) & volcano$padj < 0.05,
         "FDR < 0.05 and |LFC| < 1", "Not significant")
)
volcano$plot_direction <- ifelse(
  is.na(volcano$stat), NA_character_,
  ifelse(volcano$stat < 0, "Edge-down", "Edge-up")
)
volcano_plot <- subset(volcano, is.finite(log2FoldChange_raw) & is.finite(stat))
volcano_plot$template_class <- factor(
  volcano_plot$template_class,
  levels = c("Not significant", "FDR < 0.05 and |LFC| < 1", "FDR < 0.05 and |LFC| >= 1")
)
pde2a_volcano <- subset(volcano_plot, display_label == "PDE2A")
if (nrow(pde2a_volcano) != 1L) stop("Expected exactly one PDE2A row in volcano source data")
significant_for_labels <- subset(
  volcano_plot,
  !is.na(padj) & padj < 0.05 & display_label != "PDE2A" &
    !grepl("^ENSG[0-9]+", display_label)
)
label_parts <- lapply(split(significant_for_labels, significant_for_labels$plot_direction), function(z) {
  head(z[order(z$padj, -abs(z$stat)), , drop = FALSE], 8L)
})
volcano_labels <- do.call(rbind, label_parts)
volcano$template_auto_label <- volcano$gene_id %in% volcano_labels$gene_id
volcano$template_pde2a_label <- volcano$display_label == "PDE2A"

p_volcano <- ggplot(volcano_plot, aes(log2FoldChange_raw, stat)) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "#A9A9A9") +
  geom_vline(xintercept = 0, linetype = "22", linewidth = 0.35, colour = "#A9A9A9") +
  geom_point(aes(colour = plot_direction), alpha = 0.20, size = 0.34) +
  geom_point(data = subset(volcano_plot, !is.na(padj) & padj < 0.05),
             aes(colour = plot_direction), alpha = 0.90, size = 0.82) +
  geom_label_repel(
    data = volcano_labels,
    aes(label = display_label, colour = plot_direction),
    family = "Arial", size = FONT_MM, fill = "white", label.size = 0.18,
    box.padding = 0.30, point.padding = 0.10, min.segment.length = 0,
    segment.size = 0.28, seed = 1102, max.overlaps = Inf,
    force = 2.5, force_pull = 0.35, max.time = 5, max.iter = 20000,
    show.legend = FALSE
  ) +
  geom_label_repel(
    data = pde2a_volcano,
    aes(label = "PDE2A"),
    family = "Arial", fontface = "bold", size = FONT_MM,
    colour = "#7C3AED", fill = "white", label.size = 0.25,
    box.padding = 0.24, point.padding = 0.12, min.segment.length = 0,
    segment.color = "#7C3AED", segment.size = 0.40, seed = 2202,
    nudge_x = -1.05, nudge_y = 0.65, direction = "both", max.overlaps = Inf,
    force = 4, force_pull = 0.2, max.time = 5, max.iter = 20000,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = c("Edge-down" = "#2F5D8A", "Edge-up" = "#C00000"),
                      guide = "none") +
  coord_cartesian(xlim = c(-8.8, 10.6), ylim = c(-6.0, 5.7), expand = FALSE) +
  labs(x = "Raw log2 fold change (Edge - Core)",
       y = "Wald statistic") +
  theme_manuscript(7) +
  theme(
    plot.title = element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "#333333", fill = NA, linewidth = 0.45),
    axis.line = element_blank(),
    plot.margin = margin(5.5, 7.5, 5.5, 5.5)
  )

heat <- read_tsv(fig_source("DEG", "Figure_DEG_03_TopDEG_heatmap"))
sample_order <- unique(heat$sample_id)
gene_order <- unique(heat$gene_label)
heat$sample_display <- gsub("_TI$", "\nEdge", gsub("_C$", "\nCore", heat$sample_id))
sample_display_order <- unique(heat$sample_display)
heat$sample_display <- factor(heat$sample_display, levels = sample_display_order)
heat$gene_label <- factor(heat$gene_label, levels = rev(gene_order))

p_heat <- ggplot(heat, aes(sample_display, gene_label, fill = row_z)) +
  geom_tile(linewidth = 0.12, colour = "white") +
  scale_fill_gradient2(low = COL_CORE, mid = "white", high = "#B40426", midpoint = 0,
                       limits = c(-max(abs(heat$row_z), na.rm = TRUE), max(abs(heat$row_z), na.rm = TRUE)),
                       name = "Row z-score") +
  labs(title = "Top 30 genes by adjusted P value") +
  theme_heatmap(7) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1),
        legend.position = "right")

figure2 <- p_ma / (p_volcano | p_heat) +
  plot_layout(heights = c(0.72, 1.1), widths = c(1.05, 0.95), guides = "collect") +
  plot_annotation(
    title = "Paired transcriptome-wide differential expression",
    tag_levels = "a",
    theme = theme(plot.title = element_text(family = "Arial", face = "bold", size = FONT_PT),
                  plot.tag = element_text(family = "Arial", face = "bold", size = FONT_PT))
  ) & theme(legend.position = "top")

save_bundle(output_root, "Figure2_DifferentialExpression", figure2, 183, 148,
            list(MA = ma, volcano = volcano, heatmap = heat), script_path)

# -----------------------------------------------------------------------------
# Figure 3: functional programs
# -----------------------------------------------------------------------------
select_balanced <- function(x, n_each = 6L) {
  x <- x[!is.na(x$padj) & x$padj < 0.05, , drop = FALSE]
  parts <- lapply(split(x, x$direction), function(z) head(z[order(z$padj, -abs(z$NES)), , drop = FALSE], n_each))
  do.call(rbind, parts)
}

gsea_h_full <- read_tsv(fig_source("Enrichment", "Figure_Enrichment_01_GSEA_Hallmark"))
gsea_h <- select_balanced(gsea_h_full, 6L)
gsea_h$pathway_display_name <- standardize_pathway_name(gsea_h$pathway_label, "Hallmark")
gsea_h$label <- wrap_text(gsea_h$pathway_display_name, 27)
gsea_h$label <- factor(gsea_h$label, levels = gsea_h$label[order(gsea_h$NES)])

hallmark_axis_labels <- function(x) {
  as.expression(lapply(x, function(z) {
    if (identical(z, "Interferon-alpha Response")) {
      return(bquote(paste("Interferon-", alpha, " Response")))
    }
    if (identical(z, "Interferon-gamma Response")) {
      return(bquote(paste("Interferon-", gamma, " Response")))
    }
    z
  }))
}

p_hallmark <- ggplot(gsea_h, aes(NES, label)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "#777777") +
  geom_segment(aes(x = 0, xend = NES, yend = label, colour = direction), linewidth = 0.55) +
  geom_point(aes(colour = direction, size = minus_log10_padj), alpha = 0.95) +
  scale_colour_manual(values = c("Core-enriched" = COL_CORE, "Edge-enriched" = COL_EDGE)) +
  scale_y_discrete(labels = hallmark_axis_labels) +
  scale_size_continuous(range = c(1.8, 4.8), name = expression(-log[10]("FDR q-value"))) +
  labs(title = "MSigDB Hallmark Gene Set Enrichment Analysis",
       x = "Normalized Enrichment Score (NES; Edge vs. Core)",
       y = NULL, colour = "Enrichment Direction") +
  theme_manuscript(7) + theme(legend.position = "right")

gsea_gr_full <- read_tsv(fig_source("Enrichment", "Figure_Enrichment_02_GSEA_GO_Reactome"))
gsea_gr <- select_balanced(gsea_gr_full, 6L)
gsea_gr$label <- wrap_text(clean_term(gsea_gr$pathway_label), 32)
gsea_gr$label <- factor(gsea_gr$label, levels = gsea_gr$label[order(gsea_gr$NES)])

p_gr <- ggplot(gsea_gr, aes(NES, label)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "#777777") +
  geom_segment(aes(x = 0, xend = NES, yend = label, colour = direction), linewidth = 0.5) +
  geom_point(aes(colour = direction, size = minus_log10_padj), alpha = 0.95) +
  scale_colour_manual(values = c("Core-enriched" = COL_CORE, "Edge-enriched" = COL_EDGE)) +
  scale_size_continuous(range = c(1.6, 4.4), name = expression(-log[10](FDR))) +
  facet_grid(collection ~ ., scales = "free_y", space = "free_y") +
  labs(title = "GO BP and Reactome GSEA", x = "Normalized enrichment score (Edge - Core)", y = NULL,
       colour = NULL) +
  theme_manuscript(7.2) + theme(legend.position = "right")

ora_down_full <- read_tsv(fig_source("Enrichment", "Figure_Enrichment_03_ORA_Edge_down"))
ora_down <- head(ora_down_full[order(ora_down_full$padj, -ora_down_full$foldEnrichment), ], 12L)
ora_down$pathway_display_name <- standardize_pathway_name(ora_down$pathway_label, ora_down$collection)
ora_down$label <- wrap_text(ora_down$pathway_display_name, 31)
ora_down$label <- factor(ora_down$label, levels = rev(ora_down$label[order(ora_down$padj)]))

p_ora_down <- ggplot(ora_down, aes(foldEnrichment, label)) +
  geom_segment(aes(x = 0, xend = foldEnrichment, yend = label), linewidth = 0.5, colour = COL_LIGHT_BLUE) +
  geom_point(aes(size = overlapSize, colour = minus_log10_padj), alpha = 0.95) +
  scale_colour_gradient(low = "#8FB3D1", high = COL_DARK_BLUE, name = "-log10(BH-adjusted P value)") +
  scale_size_continuous(range = c(1.7, 4.6), name = "Overlapping genes") +
  labs(title = "Over-representation analysis: Edge-downregulated DEGs",
       x = "Fold enrichment", y = NULL) +
  theme_manuscript(7.2) + theme(legend.position = "right")

figure3 <- (p_hallmark | p_gr | p_ora_down) +
  plot_layout(widths = c(0.92, 1.15, 1.05), guides = "keep") +
  plot_annotation(
    title = "Functional programs distinguishing tumour Edge from Core",
    subtitle = "Negative NES denotes Core-enriched / Edge-decreased transcriptional programs",
    tag_levels = "a",
    theme = theme(plot.title = element_text(family = "Arial", face = "bold", size = FONT_PT),
                  plot.subtitle = element_text(family = "Arial", size = FONT_PT, colour = COL_NEUTRAL_DARK),
                  plot.tag = element_text(family = "Arial", face = "bold", size = FONT_PT))
  )

save_bundle(output_root, "Figure3_FunctionalPrograms", figure3, 183, 112,
            list(Hallmark_display = gsea_h, Hallmark_full = gsea_h_full,
                 GO_Reactome_display = gsea_gr, GO_Reactome_full = gsea_gr_full,
                 Edge_down_ORA_display = ora_down, Edge_down_ORA_full = ora_down_full), script_path)

# -----------------------------------------------------------------------------
# Figure 4: gene–pathway convergence, global and PDE2A-context versions
# -----------------------------------------------------------------------------
module_heat <- read_tsv(fig_source("Candidate_convergence", "Figure_Convergence_01_GeneModuleHeatmap"))
module_heat$functional_module_display <- standardize_functional_module(module_heat$functional_module)
gene_rank <- aggregate(selection_rank ~ gene_name, module_heat, min)
gene_levels <- gene_rank$gene_name[order(gene_rank$selection_rank, gene_rank$gene_name)]
module_totals <- aggregate(pathway_count ~ functional_module, module_heat, sum)
module_levels <- module_totals$functional_module[order(-module_totals$pathway_count, module_totals$functional_module)]
module_heat$gene_name <- factor(module_heat$gene_name, levels = rev(gene_levels))
module_heat$module_label <- wrap_text(module_heat$functional_module_display, 24)
module_label_levels <- wrap_text(standardize_functional_module(module_levels), 24)
module_heat$module_label <- factor(module_heat$module_label, levels = module_label_levels)

p_module <- ggplot(module_heat, aes(module_label, gene_name, fill = pathway_count)) +
  geom_tile(linewidth = 0.12, colour = "white") +
  scale_fill_gradient(low = "#F7FAFC", high = COL_DARK_BLUE, name = "Leading-edge\npathway count") +
  labs(title = "Leading-edge pathway membership across functional modules") +
  theme_heatmap(6.8) +
  theme(axis.text.x = element_text(angle = 48, hjust = 1, vjust = 1), legend.position = "right")

conv <- read_tsv(fig_source("Candidate_convergence", "Figure_Convergence_02_Global"))
conv$direction <- factor(conv$direction, levels = c("Edge-down", "Edge-up"))
conv$minus_log10_padj <- -log10(pmax(conv$padj, .Machine$double.xmin))

make_convergence_scatter <- function(label_pde2a = FALSE) {
  p <- ggplot(conv, aes(shrunken_log2FC, minus_log10_padj)) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "#777777") +
    geom_point(aes(colour = direction, size = total_significant_leading_edge_count), alpha = 0.86) +
    scale_colour_manual(values = c("Edge-down" = COL_CORE, "Edge-up" = COL_EDGE)) +
    scale_size_continuous(range = c(1.4, 5.6), name = "Significant leading-edge\npathway count") +
    labs(title = "Candidate context among 84 FDR-significant genes",
         x = "Shrunken log2 fold change (Edge - Core)",
         y = expression(-log[10](adjusted~italic(P))), colour = NULL) +
    theme_manuscript(8) + theme(legend.position = "right")
  if (label_pde2a) {
    p <- p +
      geom_point(data = subset(conv, gene_name == "PDE2A"), shape = 21, fill = NA,
                 colour = COL_PDE2A, stroke = 1.0,
                 aes(size = total_significant_leading_edge_count), show.legend = FALSE) +
      geom_text(data = subset(conv, gene_name == "PDE2A"), aes(label = gene_name),
                colour = COL_PDE2A, family = "Arial", fontface = "bold", size = FONT_MM,
                nudge_x = -0.18, nudge_y = 0.10, hjust = 1, vjust = 0,
                show.legend = FALSE)
  }
  p
}

p_conv_global <- make_convergence_scatter(FALSE)
p_conv_pde <- make_convergence_scatter(TRUE)

figure4_global <- p_module / p_conv_global +
  plot_layout(heights = c(1.25, 0.82)) +
  plot_annotation(title = "Gene–pathway convergence among significant regional DEGs", tag_levels = "a",
                  theme = theme(plot.title = element_text(family = "Arial", face = "bold", size = FONT_PT),
                                plot.tag = element_text(family = "Arial", face = "bold", size = FONT_PT)))
figure4_pde <- p_module / p_conv_pde +
  plot_layout(heights = c(1.25, 0.82)) +
  plot_annotation(title = "PDE2A within the frozen gene–pathway convergence landscape", tag_levels = "a",
                  theme = theme(plot.title = element_text(family = "Arial", face = "bold", size = FONT_PT),
                                plot.tag = element_text(family = "Arial", face = "bold", size = FONT_PT)))

conv_sources <- list(gene_module_heatmap = module_heat, candidate_context = conv)
save_bundle(output_root, "Figure4_GenePathwayConvergence_Global", figure4_global, 183, 166,
            conv_sources, script_path)
save_bundle(output_root, "Figure4_GenePathwayConvergence_PDE2A", figure4_pde, 183, 166,
            conv_sources, script_path)

# -----------------------------------------------------------------------------
# Figure 6: transparent candidate nomination, global and PDE2A-context versions
# -----------------------------------------------------------------------------
funnel <- read_tsv(fig_source("Candidate_nomination", "Figure_Nomination_01_Funnel"))
funnel$y <- max(funnel$stage_order) + 1 - funnel$stage_order
funnel$fill_class <- ifelse(grepl("PDE2A", funnel$stage), "PDE2A", "Evidence stage")

p_funnel <- ggplot(funnel, aes(y = y)) +
  geom_segment(data = funnel[-nrow(funnel), ], aes(x = 0.5, xend = 0.5, yend = y - 0.72),
               arrow = arrow(length = unit(1.8, "mm"), type = "closed"),
               linewidth = 0.45, colour = "#808080") +
  geom_label(aes(x = 0.5, label = paste0(stage, "\n", format(n, big.mark = ",")), fill = fill_class),
             family = "Arial", fontface = "bold", size = FONT_MM, linewidth = 0.25,
             label.padding = unit(2.1, "mm"), colour = "#222222") +
  scale_fill_manual(values = c("Evidence stage" = "#E8EEF5", "PDE2A" = "#EDE3FF")) +
  coord_cartesian(xlim = c(0, 1), clip = "off") +
  labs(title = "Transparent nomination workflow") +
  theme_void(base_family = "Arial", base_size = FONT_PT) +
  theme(text = element_text(family = "Arial", size = FONT_PT),
        plot.title = element_text(face = "bold", hjust = 0, size = FONT_PT), legend.position = "none",
        plot.margin = margin(4, 4, 4, 4))

lopo <- read_tsv(fig_source("Candidate_nomination", "Figure_Nomination_02_LOPO_Stability"))
lopo <- lopo[order(-lopo$LOPO_negative_count, lopo$LOPO_median_LFC, lopo$gene), ]
lopo_display <- head(lopo, 20L)
lfc_cols <- grep("^LOPO_P(08|09|34|36|38)_LFC$", names(lopo_display), value = TRUE)
lopo_long <- do.call(rbind, lapply(lfc_cols, function(v) data.frame(
  gene = lopo_display$gene,
  omitted_patient = gsub("^LOPO_|_LFC$", "", v),
  shrunken_LFC = lopo_display[[v]],
  stringsAsFactors = FALSE
)))
lopo_long$gene <- factor(lopo_long$gene, levels = rev(lopo_display$gene))
lopo_long$omitted_patient <- factor(lopo_long$omitted_patient, levels = c("P08", "P09", "P34", "P36", "P38"))

p_lopo <- ggplot(lopo_long, aes(omitted_patient, gene, fill = shrunken_LFC)) +
  geom_tile(linewidth = 0.15, colour = "white") +
  scale_fill_gradient2(low = COL_CORE, mid = "white", high = COL_EDGE, midpoint = 0,
                       name = "Shrunken\nlog2FC") +
  labs(title = "LOPO effect stability", subtitle = "Top 20 by median LOPO effect; all 35 retained in source data") +
  theme_heatmap(6.8) + theme(axis.text.x = element_text(angle = 0), legend.position = "right")

core_map <- read_tsv(fig_source("Candidate_nomination", "Figure_Nomination_03_LeadingEdge_CoreMap"))
annotation_cols <- c("gene", "pareto_front", "full_shrunken_LFC", "LOPO_negative_count", "patient_consistency")
module_cols <- setdiff(names(core_map), annotation_cols)
core_long <- data.frame(
  gene = rep(core_map$gene, times = length(module_cols)),
  functional_module = rep(module_cols, each = nrow(core_map)),
  pathway_count = as.numeric(unlist(core_map[module_cols], use.names = FALSE)),
  stringsAsFactors = FALSE
)
core_long <- merge(core_long, core_map[annotation_cols], by = "gene", all.x = TRUE, sort = FALSE)
core_long$functional_module_display <- standardize_functional_module(core_long$functional_module)
gene_order_core <- core_map$gene[order(core_map$pareto_front, core_map$full_shrunken_LFC, core_map$gene)]
module_total_core <- colSums(core_map[module_cols], na.rm = TRUE)
module_order_core <- names(sort(module_total_core, decreasing = TRUE))
core_long$gene <- factor(core_long$gene, levels = rev(gene_order_core))
core_long$module_label <- factor(wrap_text(core_long$functional_module_display, 22),
                                 levels = wrap_text(standardize_functional_module(module_order_core), 22))

p_core <- ggplot(core_long, aes(module_label, gene, fill = pathway_count)) +
  geom_tile(linewidth = 0.15, colour = "white") +
  scale_fill_gradient(low = "#F8F5FF", high = COL_PDE2A,
                      name = "Core-enriched leading-edge\npathway count (NES < 0)") +
  labs(title = "Core-enriched leading-edge pathway map") +
  theme_heatmap(6.6) +
  theme(axis.text.x = element_text(angle = 48, hjust = 1, vjust = 1), legend.position = "right")

pareto <- read_tsv(fig_source("Candidate_nomination", "Figure_Nomination_04_Pareto_Global"))
pareto$front_class <- ifelse(pareto$pareto_front == 1, "Pareto Front 1", "Other fronts")

make_pareto <- function(label_pde2a = FALSE) {
  p <- ggplot(pareto, aes(full_shrunken_LFC, edge_down_leading_edge_count)) +
    geom_point(aes(size = functional_module_count, fill = front_class), shape = 21,
               colour = "#FFFFFF", stroke = 0.35, alpha = 0.92) +
    scale_fill_manual(values = c("Pareto Front 1" = COL_PDE2A, "Other fronts" = "#7D8794")) +
    scale_size_continuous(range = c(1.8, 6.2), name = "Functional modules") +
    labs(title = "Non-dominated candidate context", x = "Shrunken log2 fold change (Edge - Core)",
         y = "NES < 0 significant leading-edge pathways", fill = NULL) +
    theme_manuscript(8) + theme(legend.position = "right")
  if (label_pde2a) {
    p <- p +
      geom_point(data = subset(pareto, gene == "PDE2A"), aes(size = functional_module_count),
                 shape = 21, fill = NA, colour = COL_PDE2A, stroke = 1.1, show.legend = FALSE) +
      geom_text_repel(data = subset(pareto, gene == "PDE2A"), aes(label = gene),
                      colour = COL_PDE2A, family = "Arial", fontface = "bold", size = FONT_MM,
                      seed = 1106, min.segment.length = 0, show.legend = FALSE)
  }
  p
}

p_pareto_global <- make_pareto(FALSE)
p_pareto_pde <- make_pareto(TRUE)

figure6_global <- (p_funnel | p_pareto_global) / (p_lopo | p_core) +
  plot_layout(widths = c(0.78, 1.22), heights = c(0.88, 1.12)) +
  plot_annotation(title = "Transparent robustness and pathway-core candidate nomination", tag_levels = "a",
                  theme = theme(plot.title = element_text(family = "Arial", face = "bold", size = FONT_PT),
                                plot.tag = element_text(family = "Arial", face = "bold", size = FONT_PT)))
figure6_pde <- (p_funnel | p_pareto_pde) / (p_lopo | p_core) +
  plot_layout(widths = c(0.78, 1.22), heights = c(0.88, 1.12)) +
  plot_annotation(title = "PDE2A within the transparent candidate-nomination landscape", tag_levels = "a",
                  theme = theme(plot.title = element_text(family = "Arial", face = "bold", size = FONT_PT),
                                plot.tag = element_text(family = "Arial", face = "bold", size = FONT_PT)))

nom_sources <- list(funnel = funnel, LOPO_all35 = lopo, LOPO_display = lopo_long,
                    leading_edge_core_map = core_long, Pareto_candidates = pareto)
save_bundle(output_root, "Figure6_CandidateNomination_Global", figure6_global, 183, 164,
            nom_sources, script_path)
save_bundle(output_root, "Figure6_CandidateNomination_PDE2A", figure6_pde, 183, 164,
            nom_sources, script_path)

write_tsv(data.frame(
  figure = c("Figure1_Cohort_and_PDE2A", "Figure2_DifferentialExpression",
             "Figure3_FunctionalPrograms", "Figure4_GenePathwayConvergence_Global",
             "Figure4_GenePathwayConvergence_PDE2A", "Figure6_CandidateNomination_Global",
             "Figure6_CandidateNomination_PDE2A"),
  scope = c("main", "main", "main", "unbiased_global", "main_context",
            "unbiased_global", "main_context"),
  batch_visible = FALSE,
  stringsAsFactors = FALSE
), file.path(output_root, "MAIN_FIGURE_INDEX.tsv"))

cat("Primary manuscript redraw completed in:", output_root, "\n")
