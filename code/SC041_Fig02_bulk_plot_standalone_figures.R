#!/usr/bin/env Rscript

# Render one analytical panel per output file. Multi-panel/composite figures and
# faceted multi-collection panels are explicitly prohibited. Frozen source data
# are reused as-is.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 33_plot_manuscript_standalone_v2.R <project_root>")
root <- sub("/+$", "", args[[1]])
script_path <- file.path(root, "scripts", "33_plot_manuscript_standalone_v2.R")

# Build all validated plot objects while suppressing the superseded composite exports.
Sys.setenv(PDE2A_SKIP_BUNDLE_EXPORT = "1")
sys.source(file.path(root, "scripts", "31_plot_manuscript_main_figures.R"), envir = .GlobalEnv)
sys.source(file.path(root, "scripts", "32_plot_manuscript_supplementary_figures.R"), envir = .GlobalEnv)
Sys.unsetenv("PDE2A_SKIP_BUNDLE_EXPORT")
script_path <- file.path(root, "scripts", "33_plot_manuscript_standalone_v2.R")

output_root <- file.path(root, "results", "07_figures", "manuscript_redraw_v4_Arial7pt")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

export_one <- function(order, name, plot, width_mm, height_mm, data, scope,
                       replaces, interpretation_boundary = "none") {
  data$standalone_figure <- name
  save_bundle(output_root, name, plot, width_mm, height_mm,
              list(data = data), script_path)
  data.frame(
    figure_order = order,
    figure = name,
    scope = scope,
    replaces = replaces,
    interpretation_boundary = interpretation_boundary,
    batch_visible = FALSE,
    composite = FALSE,
    stringsAsFactors = FALSE
  )
}

audit_rows <- list()

# Formal PDE2A source data retain both per-sample expression and frozen DESeq2 statistics.
pde_standalone <- pde
pde_standalone$DESeq2_shrunken_log2FC <- pde_stats$log2FoldChange_shrunk[[1]]
pde_standalone$DESeq2_pvalue <- pde_stats$pvalue[[1]]
pde_standalone$DESeq2_padj <- pde_stats$padj[[1]]
pde_standalone$patient_direction_consistency <- pde_stats$direction_consistency_edge_lower[[1]]

audit_rows[[length(audit_rows) + 1L]] <- export_one(
  1, "Figure_01_PairedSamplingDesign", p_design, 78, 58, design_points, "methods_candidate",
  "new standalone sampling schematic")
audit_rows[[length(audit_rows) + 1L]] <- export_one(
  2, "Figure_02_PairedPCA_CoreEdge", p_pca, 92, 78, pca, "main_candidate",
  "Figure_QC_Count_05_PairedPCA_CoreEdge",
  "PCA represents global transcriptomic structure, not PDE2A expression")
audit_rows[[length(audit_rows) + 1L]] <- export_one(
  3, "Figure_03_PDE2A_PairedExpression", p_pde, 82, 78, pde_standalone, "main_candidate",
  "PDE2A_paired_preanalysis_internal",
  "Normalized counts are descriptive; formal inference is the frozen paired DESeq2 model")

audit_rows[[length(audit_rows) + 1L]] <- export_one(
  4, "Figure_04_MA", p_ma, 120, 88, ma, "main_candidate", "Figure_DEG_01_MA")
audit_rows[[length(audit_rows) + 1L]] <- export_one(
  5, "Figure_05_Volcano_Global", p_volcano, 120, 92, volcano, "main_candidate",
  "Figure_DEG_02_Volcano_global", "PDE2A is not specially highlighted")
audit_rows[[length(audit_rows) + 1L]] <- export_one(
  6, "Figure_06_Top30DEG_Heatmap", p_heat, 112, 145, heat, "main_candidate",
  "Figure_DEG_03_TopDEG_heatmap", "PDE2A is not forced into the Top 30")

audit_rows[[length(audit_rows) + 1L]] <- export_one(
  7, "Figure_07_GSEA_Hallmark", p_hallmark, 135, 102, gsea_h, "main_candidate",
  "Figure_Enrichment_01_GSEA_Hallmark")

# GO BP and Reactome are deliberately exported as separate files. Combining
# them with facets would violate the one-analytical-panel-per-file contract.
make_gsea_single_collection <- function(data, title_text) {
  data$pathway_display_name <- standardize_pathway_name(data$pathway_label, data$collection)
  data$label <- wrap_text(data$pathway_display_name, 38)
  data$label <- factor(data$label, levels = data$label[order(data$NES)])
  ggplot(data, aes(NES, label)) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "#777777") +
    geom_segment(aes(x = 0, xend = NES, yend = label, colour = direction), linewidth = 0.5) +
    geom_point(aes(colour = direction, size = minus_log10_padj), alpha = 0.95) +
    scale_colour_manual(values = c("Core-enriched" = COL_CORE, "Edge-enriched" = COL_EDGE)) +
    scale_size_continuous(range = c(1.6, 4.4), name = expression(-log[10]("FDR q-value"))) +
    labs(title = title_text, x = "Normalized Enrichment Score (NES; Edge vs. Core)", y = NULL,
         colour = "Enrichment Direction") +
    theme_manuscript(7) + theme(legend.position = "right")
}
prepare_full_gsea <- function(path) {
  x <- read_tsv(path)
  x$pathway_label <- x$pathway
  x$minus_log10_padj <- -log10(pmax(x$padj, .Machine$double.xmin))
  x
}
gsea_go <- select_balanced(prepare_full_gsea(file.path(
  root, "results", "06_enrichment", "GSEA_GO_BP.tsv"
)), 6L)
gsea_reactome <- select_balanced(prepare_full_gsea(file.path(
  root, "results", "06_enrichment", "GSEA_Reactome.tsv"
)), 6L)
gsea_go$pathway_display_name <- standardize_pathway_name(gsea_go$pathway_label, "GO_BP")
gsea_reactome$pathway_display_name <- standardize_pathway_name(
  gsea_reactome$pathway_label, "Reactome"
)
p_go <- make_gsea_single_collection(gsea_go, "GO Biological Process GSEA")
p_reactome <- make_gsea_single_collection(gsea_reactome, "Reactome pathway GSEA")

audit_rows[[length(audit_rows) + 1L]] <- export_one(
  8, "Figure_08_GSEA_GO_BP", p_go, 145, 122, gsea_go, "main_candidate",
  "Figure_Enrichment_02_GSEA_GO_Reactome",
  "GO BP is kept separate from Reactome to preserve a single analytical panel")
audit_rows[[length(audit_rows) + 1L]] <- export_one(
  9, "Figure_09_GSEA_Reactome", p_reactome, 145, 122, gsea_reactome, "main_candidate",
  "Figure_Enrichment_02_GSEA_GO_Reactome",
  "Reactome is kept separate from GO BP to preserve a single analytical panel")
audit_rows[[length(audit_rows) + 1L]] <- export_one(
  10, "Figure_10_ORA_EdgeDown", p_ora_down, 140, 108, ora_down, "main_candidate",
  "Figure_Enrichment_03_ORA_Edge_down")
audit_rows[[length(audit_rows) + 1L]] <- export_one(
  11, "Figure_11_ORA_EdgeUp", p_ora_up, 140, 112, ora_up_display, "supplementary_candidate",
  "Figure_Enrichment_04_ORA_Edge_up",
  "Exploratory context from 18 Edge-up DEGs")

audit_rows[[length(audit_rows) + 1L]] <- export_one(
  12, "Figure_12_GeneModule_LeadingEdgeHeatmap", p_module, 183, 132, module_heat,
  "main_candidate", "Figure_Convergence_01_GeneModuleHeatmap")
audit_rows[[length(audit_rows) + 1L]] <- export_one(
  13, "Figure_13_Convergence_Global", p_conv_global, 112, 88, conv, "unbiased_global",
  "Figure_Convergence_02_Global", "PDE2A is not specially highlighted")
audit_rows[[length(audit_rows) + 1L]] <- export_one(
  14, "Figure_14_Convergence_PDE2A_Context", p_conv_pde, 112, 88, conv, "main_context_candidate",
  "Figure_Convergence_02_PDE2A_context",
  "Only the PDE2A label/outline differs from the frozen global version")

# Preserve the separate abundance-oriented candidate-context view from Checkpoint 6.
candidate_abundance <- read_tsv(file.path(
  root, "results", "07_figures", "manuscript_candidates", "Enrichment",
  "Figure_Candidate_Context_global", "source_data.tsv"
))
candidate_abundance$direction <- factor(candidate_abundance$direction, levels = c("Edge-down", "Edge-up"))
make_abundance_context <- function(label_pde2a = FALSE) {
  q <- ggplot(candidate_abundance, aes(shrunken_LFC, minus_log10_padj)) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "#777777") +
    geom_point(aes(colour = direction, size = log10_baseMean_plus1), alpha = 0.86) +
    scale_colour_manual(values = c("Edge-down" = COL_CORE, "Edge-up" = COL_EDGE)) +
    scale_size_continuous(range = c(1.5, 5.8), name = "log10(baseMean + 1)") +
    labs(title = "Expression-aware candidate context",
         x = "Shrunken log2 fold change (Edge - Core)",
         y = expression(-log[10](adjusted~italic(P))), colour = NULL) +
    theme_manuscript(8) + theme(legend.position = "right")
  if (label_pde2a) {
    q <- q +
      geom_point(data = subset(candidate_abundance, gene == "PDE2A"),
                 aes(size = log10_baseMean_plus1), shape = 21, fill = NA,
                 colour = COL_PDE2A, stroke = 1.05, show.legend = FALSE) +
      geom_text(data = subset(candidate_abundance, gene == "PDE2A"), aes(label = gene),
                family = "Arial", fontface = "bold", size = FONT_MM, colour = COL_PDE2A,
                nudge_x = -0.18, nudge_y = 0.10, hjust = 1, vjust = 0,
                show.legend = FALSE)
  }
  q
}
p_abundance_global <- make_abundance_context(FALSE)
p_abundance_pde <- make_abundance_context(TRUE)
audit_rows[[length(audit_rows) + 1L]] <- export_one(
  15, "Figure_15_CandidateAbundance_Global", p_abundance_global, 112, 88,
  candidate_abundance, "unbiased_global", "Figure_Candidate_Context_global")
audit_rows[[length(audit_rows) + 1L]] <- export_one(
  16, "Figure_16_CandidateAbundance_PDE2A_Context", p_abundance_pde, 112, 88,
  candidate_abundance, "main_context_candidate", "Figure_Candidate_Context_PDE2A",
  "Only the PDE2A label/outline differs from the frozen global version")

audit_rows[[length(audit_rows) + 1L]] <- export_one(
  17, "Figure_17_NominationFunnel", p_funnel, 88, 120, funnel, "main_candidate",
  "Figure_Nomination_01_Funnel", "Box geometry is non-quantitative; printed counts are authoritative")

# Reconstruct the displayed top-20 LOPO source table after the full supplementary script overwrote lopo_long.
lopo_top20 <- head(lopo[order(-lopo$LOPO_negative_count, lopo$LOPO_median_LFC, lopo$gene), ], 20L)
lopo_cols <- grep("^LOPO_P(08|09|34|36|38)_LFC$", names(lopo_top20), value = TRUE)
lopo_top20_long <- do.call(rbind, lapply(lopo_cols, function(v) data.frame(
  gene = lopo_top20$gene,
  omitted_patient = gsub("^LOPO_|_LFC$", "", v),
  shrunken_LFC = lopo_top20[[v]],
  stringsAsFactors = FALSE
)))
audit_rows[[length(audit_rows) + 1L]] <- export_one(
  18, "Figure_18_LOPO_Stability_Top20", p_lopo, 105, 112, lopo_top20_long,
  "main_candidate", "Figure_Nomination_02_LOPO_Stability",
  "Top 20 display is ordered by median LOPO effect; all 35 candidates are retained separately")
audit_rows[[length(audit_rows) + 1L]] <- export_one(
  19, "Figure_19_LOPO_Stability_All35", p_lopo_all, 105, 160, lopo_long,
  "supplementary_candidate", "full 35-candidate LOPO sensitivity",
  "Near-zero LOPO effects are retained and must not be described as stable magnitudes")
audit_rows[[length(audit_rows) + 1L]] <- export_one(
  20, "Figure_20_LeadingEdge_CoreMap", p_core, 155, 96, core_long, "main_candidate",
  "Figure_Nomination_03_LeadingEdge_CoreMap")
audit_rows[[length(audit_rows) + 1L]] <- export_one(
  21, "Figure_21_Pareto_Global", p_pareto_global, 108, 88, pareto, "unbiased_global",
  "Figure_Nomination_04_Pareto_Global", "No composite score is used")
audit_rows[[length(audit_rows) + 1L]] <- export_one(
  22, "Figure_22_Pareto_PDE2A_Context", p_pareto_pde, 108, 88, pareto,
  "main_context_candidate", "Figure_Nomination_04_PDE2A_Context",
  "Only the PDE2A label/outline differs from the frozen global version; no composite score is used")

p_network_global <- make_global_network(FALSE)
p_network_pde <- make_global_network(TRUE)
audit_rows[[length(audit_rows) + 1L]] <- export_one(
  23, "Figure_23_PPI_Global", p_network_global, 120, 102, ppi, "supplementary_candidate",
  "Figure_PPI_01_84DEG_Global", "Network centrality is descriptive; PDE2A is not a hub")
audit_rows[[length(audit_rows) + 1L]] <- export_one(
  24, "Figure_24_PPI_PDE2A_Context", p_network_pde, 120, 102, ppi,
  "supplementary_context", "Figure_PPI_02_PDE2A_Context",
  "Only the PDE2A label/outline differs from the frozen global version; PDE2A degree = 2")
audit_rows[[length(audit_rows) + 1L]] <- export_one(
  25, "Figure_25_PDE2A_LocalNeighborhood", p_local, 96, 82, local,
  "supplementary_candidate", "Figure_PPI_03_PDE2A_LocalNeighborhood",
  "Only FDR-significant DEG neighbours are shown; no external STRING shell was added")
audit_rows[[length(audit_rows) + 1L]] <- export_one(
  26, "Figure_26_CandidateEvidenceMatrix", p_evidence, 155, 170, evidence,
  "supplementary_candidate", "Figure_Prioritization_01_EvidenceMatrix",
  "Columns are scaled independently for display; no values are summed or weighted")

# Submission-facing QC plots remain standalone and contain no old/new annotation.
audit_rows[[length(audit_rows) + 1L]] <- export_one(
  27, "Figure_S01_STAR_MappingRate", p_star, 120, 82, star, "supplementary_QC",
  "Figure_QC_STAR_01_UniquelyMappedRate")
audit_rows[[length(audit_rows) + 1L]] <- export_one(
  28, "Figure_S02_featureCounts_AssignmentRate", p_fc, 120, 82, fc, "supplementary_QC",
  "Figure_QC_featureCounts_01_AssignmentRate")
audit_rows[[length(audit_rows) + 1L]] <- export_one(
  29, "Figure_S03_SampleCorrelation", p_corr, 108, 104, corr_long, "supplementary_QC",
  "Figure_QC_Count_01_SampleCorrelation")

# Checkpoint 10 is retained only as explicitly null-FDR exploratory evidence.
audit_rows[[length(audit_rows) + 1L]] <- export_one(
  30, "Figure_S04_PROGENy_Exploratory", p_prog, 105, 100, prog,
  "supplementary_null_FDR_exploratory", "Figure_Activity_01_PROGENy_Global",
  "0/14 pathways passed FDR; inferred activity is not biochemical activity")
audit_rows[[length(audit_rows) + 1L]] <- export_one(
  31, "Figure_S05_PROGENy_Heatmap_Exploratory", p_prog_heat, 110, 104, prog_heat,
  "supplementary_null_FDR_exploratory", "Figure_Activity_02_PROGENy_PairedHeatmap",
  "No FDR-significant pathway activity; row z-scores are visualization only")
audit_rows[[length(audit_rows) + 1L]] <- export_one(
  32, "Figure_S06_TF_Activity_Exploratory", p_tf, 110, 125, tf,
  "supplementary_null_FDR_exploratory", "Figure_Activity_03_TF_Global",
  "0/717 TF activities passed FDR; no causal mechanism is inferred")

audit <- do.call(rbind, audit_rows)
write_tsv(audit, file.path(output_root, "STANDALONE_FIGURE_INDEX.tsv"))

# The superseded composite directory is kept for auditability but explicitly excluded.
v1_root <- file.path(root, "results", "07_figures", "manuscript_redraw_v1")
v1_png <- if (dir.exists(v1_root)) list.files(v1_root, pattern = "\\.png$", recursive = TRUE) else character()
v1_status <- data.frame(
  file = v1_png,
  status = rep("EXCLUDED_NO_COMPOSITES", length(v1_png)),
  reason = rep("User requires one analytical panel per figure; multi-panel assemblies are prohibited", length(v1_png)),
  stringsAsFactors = FALSE
)
write_tsv(v1_status, file.path(v1_root, "COMPOSITE_EXCLUSION_STATUS.tsv"))

cat("Standalone manuscript redraw completed:", nrow(audit), "figures in", output_root, "\n")
