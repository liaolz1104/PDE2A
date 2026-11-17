#!/usr/bin/env Rscript

# Rebuild supplementary and internal-evidence figures from frozen source tables.
# Batch annotations are deliberately excluded from manuscript-facing outputs.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 32_plot_manuscript_supplementary_figures.R <project_root>")
root <- sub("/+$", "", args[[1]])
script_path <- file.path(root, "scripts", "32_plot_manuscript_supplementary_figures.R")
source(file.path(root, "scripts", "30_manuscript_redraw_helpers.R"))

output_root <- file.path(root, "results", "07_figures", "manuscript_redraw_v1")
candidate_root <- file.path(root, "results", "07_figures", "manuscript_candidates")
internal_root <- file.path(root, "results", "07_figures", "internal_QC")

fig_source <- function(section, figure) {
  file.path(candidate_root, section, figure, "source_data.tsv")
}

# -----------------------------------------------------------------------------
# Supplementary Figure S1: submission-facing sequencing/count QC and structure
# -----------------------------------------------------------------------------
star <- read_tsv(file.path(internal_root, "STAR", "Figure_QC_STAR_01_UniquelyMappedRate", "source_data.tsv"))
fc <- read_tsv(file.path(internal_root, "featureCounts", "Figure_QC_featureCounts_01_AssignmentRate", "source_data.tsv"))
sample_order <- c("P08_C", "P08_TI", "P08_P", "P09_C", "P09_TI", "P09_P",
                  "P34_C", "P34_TI", "P36_C", "P36_TI", "P38_C", "P38_TI")
star$sample_id <- factor(star$sample_id, levels = sample_order)
fc$sample_id <- factor(fc$sample_id, levels = sample_order)
star$region_label <- region_label(star$region)
fc$region_label <- region_label(fc$region)

p_star <- ggplot(star, aes(sample_id, `Uniquely mapped reads %`, fill = region_label)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.1f", `Uniquely mapped reads %`)), vjust = -0.35,
            size = FONT_MM, family = "Arial") +
  region_fill_scale +
  coord_cartesian(ylim = c(0, 101), clip = "off") +
  labs(title = "STAR uniquely mapped reads", x = NULL, y = "Uniquely mapped reads (%)", fill = "Region") +
  theme_manuscript(7.5) +
  theme(axis.text.x = element_text(angle = 52, hjust = 1), legend.position = "top")

p_fc <- ggplot(fc, aes(sample_id, assigned_percent, fill = region_label)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.1f", assigned_percent)), vjust = -0.35,
            size = FONT_MM, family = "Arial") +
  region_fill_scale +
  coord_cartesian(ylim = c(0, 101), clip = "off") +
  labs(title = "Gene-level assignment", x = NULL, y = "Assigned fragments (%)", fill = "Region") +
  theme_manuscript(7.5) +
  theme(axis.text.x = element_text(angle = 52, hjust = 1), legend.position = "top")

corr <- read_tsv(file.path(root, "results", "04_featurecounts", "sample_correlation_Pearson.tsv"))
corr_samples <- corr$sample_id
corr_mat <- as.matrix(corr[, setdiff(names(corr), "sample_id"), drop = FALSE])
storage.mode(corr_mat) <- "numeric"
corr_long <- data.frame(
  sample_y = rep(corr_samples, times = ncol(corr_mat)),
  sample_x = rep(colnames(corr_mat), each = nrow(corr_mat)),
  Pearson_r = as.vector(corr_mat),
  stringsAsFactors = FALSE
)
corr_long$sample_x <- factor(corr_long$sample_x, levels = sample_order)
corr_long$sample_y <- factor(corr_long$sample_y, levels = rev(sample_order))

p_corr <- ggplot(corr_long, aes(sample_x, sample_y, fill = Pearson_r)) +
  geom_tile(linewidth = 0.16, colour = "white") +
  scale_fill_gradient(low = "#F4F7FA", high = COL_DARK_BLUE, limits = c(0.8, 1),
                      oob = squish, name = "Pearson r") +
  labs(title = "VST sample correlation") +
  theme_heatmap(6.7) +
  theme(axis.text.x = element_text(angle = 52, hjust = 1), legend.position = "right")

pca <- read_tsv(file.path(root, "results", "04_featurecounts", "PCA_main10_paired_coordinates.tsv"))
pca$region_label <- region_label(pca$region)
p_pca <- ggplot(pca, aes(PC1, PC2)) +
  geom_path(aes(group = patient_id), colour = "#B0B0B0", linewidth = 0.5) +
  geom_point(aes(colour = region_label), size = 2.7) +
  geom_text_repel(aes(label = sample_id), family = "Arial", size = FONT_MM, seed = 1201,
                  min.segment.length = 0, box.padding = 0.25, max.overlaps = Inf,
                  colour = "#222222") +
  scale_colour_manual(values = c("Core" = COL_CORE, "Edge" = COL_EDGE),
                      breaks = c("Core", "Edge"), drop = TRUE) +
  labs(title = "Paired Core-Edge PCA",
       x = sprintf("PC1 (%.1f%%)", unique(pca$PC1_variance_percent)[[1]]),
       y = sprintf("PC2 (%.1f%%)", unique(pca$PC2_variance_percent)[[1]]), colour = "Region") +
  theme_manuscript(7.5) + theme(legend.position = "top")

figure_s1 <- (p_star | p_fc) / (p_corr | p_pca) +
  plot_layout(guides = "collect", heights = c(0.84, 1.16)) +
  plot_annotation(title = "Sequencing, quantification and sample-level quality control", tag_levels = "a",
                  theme = theme(plot.title = element_text(family = "Arial", face = "bold", size = FONT_PT),
                                plot.tag = element_text(family = "Arial", face = "bold", size = FONT_PT))) &
  theme(legend.position = "top")

save_bundle(output_root, "Supplementary_Figure_S1_QualityControl", figure_s1, 183, 142,
            list(STAR_mapping = star, featureCounts_assignment = fc,
                 sample_correlation = corr_long, paired_PCA = pca), script_path)

# -----------------------------------------------------------------------------
# Supplementary Figure S2: Edge-up ORA retained as a small-gene-set analysis
# -----------------------------------------------------------------------------
ora_up <- read_tsv(fig_source("Enrichment", "Figure_Enrichment_04_ORA_Edge_up"))
ora_up_display <- head(ora_up[order(ora_up$padj, -ora_up$foldEnrichment), ], 15L)
ora_up_display$pathway_display_name <- standardize_pathway_name(
  ora_up_display$pathway_label, ora_up_display$collection
)
ora_up_display$label <- ora_up_display$pathway_display_name
duplicate_label <- duplicated(ora_up_display$label) | duplicated(ora_up_display$label, fromLast = TRUE)
ora_up_display$label[duplicate_label] <- paste0(
  ora_up_display$label[duplicate_label], " [", ora_up_display$collection[duplicate_label], "]"
)
ora_up_display$label <- wrap_text(ora_up_display$label, 34)
ora_up_display$label <- factor(ora_up_display$label,
                               levels = rev(ora_up_display$label[order(ora_up_display$padj)]))

p_ora_up <- ggplot(ora_up_display, aes(foldEnrichment, label)) +
  geom_segment(aes(x = 0, xend = foldEnrichment, yend = label), linewidth = 0.55, colour = COL_LIGHT_RED) +
  geom_point(aes(size = overlapSize, colour = minus_log10_padj), alpha = 0.95) +
  scale_colour_gradient(low = "#E9A89C", high = "#9E2B1D", name = "-log10(BH-adjusted P value)") +
  scale_size_continuous(range = c(1.8, 4.8), name = "Overlapping genes") +
  labs(title = "Over-representation analysis: Edge-upregulated DEGs",
       subtitle = "Exploratory analysis of 18 FDR-significant Edge-upregulated genes",
       x = "Fold enrichment", y = NULL) +
  theme_manuscript(8) + theme(legend.position = "right")

save_bundle(output_root, "Supplementary_Figure_S2_EdgeUpORA", p_ora_up, 125, 92,
            list(Edge_up_ORA_display = ora_up_display, Edge_up_ORA_full = ora_up), script_path)

# -----------------------------------------------------------------------------
# Supplementary Figure S3: global STRING context and PDE2A neighbourhood
# -----------------------------------------------------------------------------
ppi <- read_tsv(fig_source("Network_prioritization", "Figure_PPI_01_84DEG_Global"))
nodes <- subset(ppi, record_type == "node")
edges <- subset(ppi, record_type == "edge")
nodes$direction <- factor(nodes$direction, levels = c("Edge-down", "Edge-up"))
top_labels <- head(nodes[order(-nodes$degree, nodes$padj, nodes$gene), ], 6L)

make_global_network <- function(label_pde2a = FALSE) {
  p <- ggplot() +
    geom_segment(data = edges, aes(x, y, xend = xend, yend = yend, alpha = combined_score),
                 linewidth = 0.35, colour = "#A6ADB4") +
    geom_point(data = nodes, aes(x, y, fill = direction, size = degree), shape = 21,
               colour = "white", stroke = 0.25, alpha = 0.95) +
    scale_fill_manual(values = c("Edge-down" = COL_CORE, "Edge-up" = COL_EDGE)) +
    scale_size_continuous(range = c(1.4, 5.6), name = "Degree") +
    scale_alpha_continuous(range = c(0.18, 0.72), guide = "none") +
    geom_text_repel(data = top_labels, aes(x, y, label = gene), family = "Arial", size = FONT_MM,
                    seed = 1203, min.segment.length = 0, box.padding = 0.55,
                    point.padding = 0.35, force = 4, max.iter = 20000,
                    max.overlaps = Inf, colour = "#222222") +
    coord_equal(clip = "off") +
    labs(title = "FDR-significant DEG STRING network", fill = NULL) +
    theme_void(base_family = "Arial", base_size = FONT_PT) +
    theme(text = element_text(family = "Arial", size = FONT_PT),
          plot.title = element_text(face = "bold", hjust = 0, size = FONT_PT),
          legend.position = "right", plot.margin = margin(5, 8, 5, 5))
  if (label_pde2a) {
    p <- p +
      geom_point(data = subset(nodes, gene == "PDE2A"), aes(x, y, size = degree), shape = 21,
                 fill = NA, colour = COL_PDE2A, stroke = 1.1, show.legend = FALSE) +
      geom_text_repel(data = subset(nodes, gene == "PDE2A"), aes(x, y, label = gene),
                      family = "Arial", fontface = "bold", size = FONT_MM, colour = COL_PDE2A,
                      seed = 1204, min.segment.length = 0, show.legend = FALSE)
  }
  p
}

local <- read_tsv(fig_source("Network_prioritization", "Figure_PPI_03_PDE2A_LocalNeighborhood"))
local_nodes <- subset(local, record_type == "node")
local_edges <- subset(local, record_type == "edge")
local_nodes$direction <- factor(local_nodes$direction, levels = c("Edge-down", "Edge-up"))
local_nodes$node_label <- sprintf("%s\nLFC %.2f", local_nodes$gene, local_nodes$shrunken_log2FC)

p_local <- ggplot() +
  geom_segment(data = local_edges, aes(x, y, xend = xend, yend = yend, alpha = combined_score),
               linewidth = 0.65, colour = "#9AA3AB") +
  geom_point(data = local_nodes, aes(x, y, fill = direction, size = degree), shape = 21,
             colour = "white", stroke = 0.35) +
  geom_point(data = subset(local_nodes, gene == "PDE2A"), aes(x, y, size = degree), shape = 21,
             fill = NA, colour = COL_PDE2A, stroke = 1.1, show.legend = FALSE) +
  geom_text_repel(data = local_nodes, aes(x, y, label = node_label), family = "Arial", size = FONT_MM,
                  seed = 1205, min.segment.length = 0, box.padding = 0.35, max.overlaps = Inf) +
  scale_fill_manual(values = c("Edge-down" = COL_CORE, "Edge-up" = COL_EDGE)) +
  scale_size_continuous(range = c(2.6, 6.4), guide = "none") +
  scale_alpha_continuous(range = c(0.35, 0.85), guide = "none") +
  coord_equal(clip = "off") +
  labs(title = "PDE2A first-degree DEG neighbourhood", fill = NULL) +
  theme_void(base_family = "Arial", base_size = FONT_PT) +
  theme(text = element_text(family = "Arial", size = FONT_PT),
        plot.title = element_text(face = "bold", hjust = 0, size = FONT_PT),
        legend.position = "none", plot.margin = margin(5, 12, 5, 8))

ppi_global <- (make_global_network(FALSE) | p_local) +
  plot_layout(widths = c(1.35, 0.65)) +
  plot_annotation(title = "Protein-interaction context of regional DEGs", tag_levels = "a",
                  theme = theme(plot.title = element_text(family = "Arial", face = "bold", size = FONT_PT),
                                plot.tag = element_text(family = "Arial", face = "bold", size = FONT_PT)))
ppi_context <- (make_global_network(TRUE) | p_local) +
  plot_layout(widths = c(1.35, 0.65)) +
  plot_annotation(title = "PDE2A within the frozen regional-DEG interaction network", tag_levels = "a",
                  theme = theme(plot.title = element_text(family = "Arial", face = "bold", size = FONT_PT),
                                plot.tag = element_text(family = "Arial", face = "bold", size = FONT_PT)))

ppi_sources <- list(global_network = ppi, local_neighbourhood = local)
save_bundle(output_root, "Supplementary_Figure_S3_PPI_Global", ppi_global, 183, 92,
            ppi_sources, script_path)
save_bundle(output_root, "Supplementary_Figure_S3_PPI_PDE2A", ppi_context, 183, 92,
            ppi_sources, script_path)

# -----------------------------------------------------------------------------
# Supplementary Figure S4: transparent evidence matrix (no composite score)
# -----------------------------------------------------------------------------
evidence <- read_tsv(fig_source("Network_prioritization", "Figure_Prioritization_01_EvidenceMatrix"))
field_levels <- c("shrunken_log2FC", "patient_consistency", "leading_edge_count",
                  "functional_module_count", "PPI_degree", "enzyme_yes_no")
field_labels <- c("Shrunken LFC", "Patient consistency", "Leading-edge count",
                  "Module count", "PPI degree", "Enzyme")
evidence$evidence_field <- factor(evidence$evidence_field, levels = field_levels, labels = field_labels)
gene_order <- unique(evidence$gene[order(evidence$display_order)])
evidence$gene <- factor(evidence$gene, levels = rev(gene_order))
evidence$value_label <- ifelse(as.character(evidence$evidence_field) == "Shrunken LFC",
                               sprintf("%.2f", evidence$raw_value), sprintf("%.0f", evidence$raw_value))

p_evidence <- ggplot(evidence, aes(evidence_field, gene, fill = within_field_display_scale)) +
  geom_tile(linewidth = 0.18, colour = "white") +
  geom_text(aes(label = value_label), family = "Arial", size = FONT_MM, colour = "#17202A") +
  scale_fill_gradient(low = "#F7FAFC", high = "#3B76AF", name = "Within-column\ndisplay scale") +
  labs(title = "Evidence matrix for Edge-down signaling/regulatory candidates",
       subtitle = "Ordered by leading-edge count, then adjusted P; columns are scaled separately for display",
       caption = "No values were summed or weighted.") +
  theme_heatmap(7.2) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        plot.caption = element_text(hjust = 0, size = FONT_PT, colour = COL_NEUTRAL_DARK),
        legend.position = "right")

save_bundle(output_root, "Supplementary_Figure_S4_EvidenceMatrix", p_evidence, 183, 170,
            list(evidence_matrix = evidence), script_path)

# -----------------------------------------------------------------------------
# Supplementary Figure S5: all-candidate LOPO sensitivity
# -----------------------------------------------------------------------------
lopo <- read_tsv(fig_source("Candidate_nomination", "Figure_Nomination_02_LOPO_Stability"))
lopo <- lopo[order(-lopo$LOPO_negative_count, lopo$LOPO_median_LFC, lopo$gene), ]
lfc_cols <- grep("^LOPO_P(08|09|34|36|38)_LFC$", names(lopo), value = TRUE)
lopo_long <- do.call(rbind, lapply(lfc_cols, function(v) data.frame(
  gene = lopo$gene,
  omitted_patient = gsub("^LOPO_|_LFC$", "", v),
  shrunken_LFC = lopo[[v]],
  stringsAsFactors = FALSE
)))
lopo_long$gene <- factor(lopo_long$gene, levels = rev(lopo$gene))
lopo_long$omitted_patient <- factor(lopo_long$omitted_patient, levels = c("P08", "P09", "P34", "P36", "P38"))

p_lopo_all <- ggplot(lopo_long, aes(omitted_patient, gene, fill = shrunken_LFC)) +
  geom_tile(linewidth = 0.16, colour = "white") +
  scale_fill_gradient2(low = COL_CORE, mid = "white", high = COL_EDGE, midpoint = 0,
                       na.value = "#D9D9D9", name = "Shrunken\nlog2FC") +
  labs(title = "LOPO effects for 35 frozen regulatory candidates",
       subtitle = "Negative values denote lower expression at the tumour Edge; grey = not estimable") +
  theme_heatmap(7) + theme(axis.text.x = element_text(angle = 0), legend.position = "right")

save_bundle(output_root, "Supplementary_Figure_S5_LOPO_AllCandidates", p_lopo_all, 128, 160,
            list(LOPO_candidate_summary = lopo, LOPO_effects = lopo_long), script_path)

# -----------------------------------------------------------------------------
# Supplementary Figure S6: exploratory activity inference, explicitly null-FDR
# -----------------------------------------------------------------------------
activity_root <- file.path(root, "results", "07_figures", "regulatory_activity")
prog <- read_tsv(file.path(activity_root, "Figure_Activity_01_PROGENy_Global_source_data.tsv"))
prog$direction <- ifelse(prog$delta_Edge_minus_Core < 0, "Edge-decreased", "Edge-increased")
prog$pathway <- factor(prog$pathway, levels = prog$pathway[order(prog$delta_Edge_minus_Core)])

p_prog <- ggplot(prog, aes(delta_Edge_minus_Core, pathway)) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "#777777") +
  geom_errorbar(aes(xmin = CI95_low, xmax = CI95_high), orientation = "y", width = 0.18,
                linewidth = 0.45, colour = "#808080") +
  geom_point(aes(colour = direction), size = 2.2) +
  scale_colour_manual(values = c("Edge-decreased" = COL_CORE, "Edge-increased" = COL_EDGE)) +
  labs(title = "Transcriptome-inferred PROGENy activity",
       subtitle = "Exploratory: 0/14 pathways passed BH-FDR < 0.05",
       x = "Inferred activity difference (Edge - Core)", y = NULL, colour = NULL) +
  theme_manuscript(7.5) + theme(legend.position = "bottom")

prog_heat <- read_tsv(file.path(activity_root, "Figure_Activity_02_PROGENy_PairedHeatmap_source_data.tsv"))
sample_levels <- unique(prog_heat$sample_id)
pathway_levels <- unique(prog_heat$pathway)
prog_heat$sample_label <- gsub("_TI$", "\nEdge", gsub("_C$", "\nCore", prog_heat$sample_id))
prog_heat$sample_label <- factor(prog_heat$sample_label, levels = unique(prog_heat$sample_label))
prog_heat$pathway <- factor(prog_heat$pathway, levels = rev(pathway_levels))

p_prog_heat <- ggplot(prog_heat, aes(sample_label, pathway, fill = row_zscore)) +
  geom_tile(linewidth = 0.14, colour = "white") +
  scale_fill_gradient2(low = COL_CORE, mid = "white", high = "#B40426", midpoint = 0,
                       name = "Row z-score") +
  labs(title = "Sample-level inferred pathway activity") +
  theme_heatmap(6.8) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1), legend.position = "right")

tf <- read_tsv(file.path(activity_root, "Figure_Activity_03_TF_Global_source_data.tsv"))
tf$direction <- ifelse(tf$delta_Edge_minus_Core < 0, "Edge-decreased", "Edge-increased")
tf$TF <- factor(tf$TF, levels = tf$TF[order(tf$delta_Edge_minus_Core)])

p_tf <- ggplot(tf, aes(delta_Edge_minus_Core, TF)) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "#777777") +
  geom_errorbar(aes(xmin = CI95_low, xmax = CI95_high), orientation = "y", width = 0.18,
                linewidth = 0.42, colour = "#808080") +
  geom_point(aes(colour = direction), size = 2.0) +
  scale_colour_manual(values = c("Edge-decreased" = COL_CORE, "Edge-increased" = COL_EDGE)) +
  labs(title = "Top TF activities by absolute moderated t",
       subtitle = "Exploratory: 0/717 TF activities passed BH-FDR < 0.05",
       x = "Inferred activity difference (Edge - Core)", y = NULL, colour = NULL) +
  theme_manuscript(7.2) + theme(legend.position = "bottom")

figure_s6 <- p_prog / (p_prog_heat | p_tf) +
  plot_layout(heights = c(0.78, 1.22)) +
  plot_annotation(title = "Exploratory regulatory-activity inference", tag_levels = "a",
                  theme = theme(plot.title = element_text(family = "Arial", face = "bold", size = FONT_PT),
                                plot.tag = element_text(family = "Arial", face = "bold", size = FONT_PT)))

save_bundle(output_root, "Supplementary_Figure_S6_RegulatoryActivity_Exploratory", figure_s6, 183, 158,
            list(PROGENy_effects = prog, PROGENy_heatmap = prog_heat, TF_effects = tf), script_path)

write_tsv(data.frame(
  figure = c("Supplementary_Figure_S1_QualityControl", "Supplementary_Figure_S2_EdgeUpORA",
             "Supplementary_Figure_S3_PPI_Global", "Supplementary_Figure_S3_PPI_PDE2A",
             "Supplementary_Figure_S4_EvidenceMatrix", "Supplementary_Figure_S5_LOPO_AllCandidates",
             "Supplementary_Figure_S6_RegulatoryActivity_Exploratory"),
  scope = c("supplementary_QC", "supplementary_exploratory", "supplementary_global",
            "supplementary_context", "supplementary_transparent_evidence",
            "supplementary_sensitivity", "supplementary_null_FDR_exploratory"),
  batch_visible = FALSE,
  stringsAsFactors = FALSE
), file.path(output_root, "SUPPLEMENTARY_FIGURE_INDEX.tsv"))

cat("Supplementary manuscript redraw completed in:", output_root, "\n")
