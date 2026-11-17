#!/usr/bin/env Rscript

# Checkpoint 5: paired transcriptome-wide DESeq2 analysis, Edge/TI - Core.
# This script stops before enrichment or candidate-family analysis.

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
local_lib <- file.path(project_root, ".Rlib_checkpoint5")
if (dir.exists(local_lib)) .libPaths(c(normalizePath(local_lib), .libPaths()))

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
})

count_file <- file.path(project_root, "results", "04_featurecounts", "count_matrix_gene_id.tsv")
annotation_file <- file.path(project_root, "results", "04_featurecounts", "gene_annotation.tsv")
manifest_file <- file.path(project_root, "manifest_combined_standardized.tsv")
plot_script <- file.path(project_root, "scripts", "16_plot_checkpoint5_deg.R")
for (path in c(count_file, annotation_file, manifest_file, plot_script)) {
  if (!file.exists(path) || file.info(path)$size == 0) stop("Missing or empty required input: ", path)
}

out_dir <- file.path(project_root, "results", "05_deseq2")
fig_root <- file.path(project_root, "results", "07_figures", "manuscript_candidates", "DEG")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_root, recursive = TRUE, showWarnings = FALSE)

sample_order <- c("P08_C", "P08_TI", "P09_C", "P09_TI", "P34_C", "P34_TI", "P36_C", "P36_TI", "P38_C", "P38_TI")
frozen_all12_order <- c("P08_C", "P08_TI", "P08_P", "P09_C", "P09_TI", "P09_P", "P34_C", "P34_TI", "P36_C", "P36_TI", "P38_C", "P38_TI")
patients <- c("P08", "P09", "P34", "P36", "P38")
counts_df <- read.delim(count_file, check.names = FALSE, stringsAsFactors = FALSE)
annotation <- read.delim(annotation_file, check.names = FALSE, stringsAsFactors = FALSE)
manifest <- read.delim(manifest_file, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(counts_df) != 62700L) stop("Frozen total gene count changed: ", nrow(counts_df))
if (anyDuplicated(counts_df$gene_id)) stop("Duplicate gene_id in count matrix")
if (!all(frozen_all12_order %in% names(counts_df))) stop("Frozen count-matrix sample columns are missing")

counts_all <- as.matrix(counts_df[, frozen_all12_order, drop = FALSE])
storage.mode(counts_all) <- "integer"
rownames(counts_all) <- counts_df$gene_id
if (anyNA(counts_all) || any(counts_all < 0L)) stop("Count matrix contains invalid values")
# Reproduce the already-frozen Checkpoint 4 filter on all 12 QC samples.
# Only the retained gene rows and the 10 prespecified Core/Edge columns enter DESeq2.
keep <- rowSums(counts_all >= 10L) >= 3L
if (sum(keep) != 20413L) stop("Low-count gate mismatch: retained ", sum(keep), " rather than frozen 20413")
counts_filtered <- counts_all[keep, sample_order, drop = FALSE]

meta <- manifest[match(sample_order, manifest$sample_id), c("sample_id", "patient_id", "region"), drop = FALSE]
if (anyNA(meta$sample_id) || !identical(meta$sample_id, sample_order)) stop("Frozen manifest mapping failed")
if (any(meta$region == "p")) stop("Paracancerous sample entered the primary model")
meta$patient_id <- factor(meta$patient_id, levels = patients)
meta$region <- factor(meta$region, levels = c("c", "ti"))
rownames(meta) <- meta$sample_id

dds <- DESeqDataSetFromMatrix(countData = counts_filtered, colData = meta, design = ~ patient_id + region)
dds <- DESeq(dds, quiet = TRUE)
coef_names <- resultsNames(dds)
coef_name <- "region_ti_vs_c"
if (!coef_name %in% coef_names) stop("Expected coefficient missing. resultsNames: ", paste(coef_names, collapse = ", "))
res_raw <- results(dds, contrast = c("region", "ti", "c"), alpha = 0.05)
res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")

normalized <- counts(dds, normalized = TRUE)
core_samples <- paste0(patients, "_C")
edge_samples <- paste0(patients, "_TI")
lower_count <- rowSums(normalized[, edge_samples, drop = FALSE] < normalized[, core_samples, drop = FALSE])
higher_count <- rowSums(normalized[, edge_samples, drop = FALSE] > normalized[, core_samples, drop = FALSE])

res <- data.frame(
  gene_id = rownames(res_raw),
  baseMean = res_raw$baseMean,
  log2FoldChange_raw = res_raw$log2FoldChange,
  lfcSE = res_raw$lfcSE,
  stat = res_raw$stat,
  pvalue = res_raw$pvalue,
  padj = res_raw$padj,
  log2FoldChange_shrunk = res_shrunk$log2FoldChange,
  mean_normalized_Core = rowMeans(normalized[, core_samples, drop = FALSE]),
  mean_normalized_Edge = rowMeans(normalized[, edge_samples, drop = FALSE]),
  direction_consistency_edge_lower = as.integer(lower_count),
  direction_consistency_edge_higher = as.integer(higher_count),
  stringsAsFactors = FALSE
)
res <- merge(res, annotation[, c("gene_id", "gene_id_stripped", "gene_name", "gene_type")], by = "gene_id", all.x = TRUE, sort = FALSE)
res <- res[match(rownames(res_raw), res$gene_id), ]
res <- res[, c("gene_id", "gene_id_stripped", "gene_name", "gene_type", "baseMean", "log2FoldChange_raw", "lfcSE", "stat", "pvalue", "padj", "log2FoldChange_shrunk", "mean_normalized_Core", "mean_normalized_Edge", "direction_consistency_edge_lower", "direction_consistency_edge_higher")]
write.table(res, file.path(out_dir, "DESeq2_Edge_vs_Core_all_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

rank_min <- function(x, decreasing = FALSE) {
  value <- if (decreasing) -x else x
  rank(value, na.last = "keep", ties.method = "min")
}
ranks <- res
ranks$rank_padj_ascending <- rank_min(ranks$padj)
ranks$rank_abs_Wald_stat_descending <- rank_min(abs(ranks$stat), decreasing = TRUE)
ranks$rank_Edge_down_stat_ascending <- rank_min(ranks$stat)
ranks$rank_Edge_up_stat_descending <- rank_min(ranks$stat, decreasing = TRUE)
ranks$rank_abs_shrunken_LFC_descending <- rank_min(abs(ranks$log2FoldChange_shrunk), decreasing = TRUE)
write.table(ranks, file.path(out_dir, "GENOMEWIDE_GENE_RANKS.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

direction_table <- res[, c("gene_id", "gene_id_stripped", "gene_name", "gene_type", "direction_consistency_edge_lower", "direction_consistency_edge_higher", "stat", "padj", "log2FoldChange_shrunk", "baseMean")]
names(direction_table)[names(direction_table) == "direction_consistency_edge_lower"] <- "direction_count_edge_lower"
names(direction_table)[names(direction_table) == "direction_consistency_edge_higher"] <- "direction_count_edge_higher"
write.table(direction_table, file.path(out_dir, "GENE_DIRECTION_CONSISTENCY.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
direction_summary <- data.frame(direction_count_edge_lower = 0:5, gene_count = as.integer(table(factor(direction_table$direction_count_edge_lower, levels = 0:5))))
write.table(direction_summary, file.path(out_dir, "GENE_DIRECTION_CONSISTENCY_SUMMARY.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pde_rows <- which(!is.na(ranks$gene_name) & ranks$gene_name == "PDE2A")
if (length(pde_rows) != 1L) stop("PDE2A does not map uniquely in the tested transcriptome: ", length(pde_rows))
pde <- ranks[pde_rows, , drop = FALSE]
pde$genomewide_tested_gene_count <- nrow(ranks)
pde$edge_down_rank_among_negative_stat <- if (is.na(pde$stat) || pde$stat >= 0) NA_integer_ else sum(!is.na(res$stat) & res$stat <= pde$stat)
write.table(pde, file.path(out_dir, "PDE2A_GENOMEWIDE_POSITION.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# Figure source data. PDE2A is not specially labeled or inserted.
ma_dir <- file.path(fig_root, "Figure_DEG_01_MA")
volcano_dir <- file.path(fig_root, "Figure_DEG_02_Volcano_global")
heatmap_dir <- file.path(fig_root, "Figure_DEG_03_TopDEG_heatmap")
for (d in c(ma_dir, volcano_dir, heatmap_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
ma <- res[, c("gene_id", "gene_name", "baseMean", "log2FoldChange_raw", "padj")]
ma$significance <- ifelse(!is.na(ma$padj) & ma$padj < 0.05, "padj < 0.05", "NS")
write.table(ma, file.path(ma_dir, "source_data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

volcano <- res[, c("gene_id", "gene_name", "log2FoldChange_shrunk", "stat", "padj")]
volcano$minus_log10_padj <- ifelse(is.na(volcano$padj), NA_real_, -log10(pmax(volcano$padj, .Machine$double.xmin)))
volcano$class <- "NS"
volcano$class[!is.na(volcano$padj) & volcano$padj < 0.05 & volcano$log2FoldChange_shrunk <= -1] <- "Edge-down"
volcano$class[!is.na(volcano$padj) & volcano$padj < 0.05 & volcano$log2FoldChange_shrunk >= 1] <- "Edge-up"
label_candidates <- which(!is.na(volcano$padj) & !is.na(volcano$stat))
label_order <- label_candidates[order(volcano$padj[label_candidates], -abs(volcano$stat[label_candidates]), na.last = NA)]
label_order <- head(label_order, 15L)
volcano$auto_label <- seq_len(nrow(volcano)) %in% label_order
volcano$display_label <- ifelse(is.na(volcano$gene_name) | volcano$gene_name == "", volcano$gene_id, volcano$gene_name)
write.table(volcano, file.path(volcano_dir, "source_data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

sig_n <- sum(!is.na(res$padj) & res$padj < 0.05)
if (sig_n >= 30L) {
  top_idx <- order(res$padj, -abs(res$stat), na.last = NA)[seq_len(30L)]
  heatmap_rule <- "Top 30 genes by adjusted P value"
} else {
  top_idx <- order(-abs(res$stat), res$padj, na.last = NA)[seq_len(30L)]
  heatmap_rule <- "Top 30 genes by absolute Wald statistic because fewer than 30 genes had adjusted P < 0.05"
}
vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)[res$gene_id[top_idx], sample_order, drop = FALSE]
row_z <- t(scale(t(vst_mat)))
labels <- ifelse(is.na(res$gene_name[top_idx]) | res$gene_name[top_idx] == "", res$gene_id[top_idx], res$gene_name[top_idx])
labels <- make.unique(labels)
heatmap_long <- do.call(rbind, lapply(seq_along(top_idx), function(i) {
  data.frame(gene_id = res$gene_id[top_idx[i]], gene_name = res$gene_name[top_idx[i]], gene_label = labels[i],
             sample_id = sample_order, row_z = as.numeric(row_z[i, sample_order]),
             patient_id = as.character(meta[sample_order, "patient_id"]),
             region_label = ifelse(as.character(meta[sample_order, "region"]) == "c", "Core", "Edge"),
             selection_rule = heatmap_rule, stringsAsFactors = FALSE)
}))
write.table(heatmap_long, file.path(heatmap_dir, "source_data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

for (d in c(ma_dir, volcano_dir, heatmap_dir)) {
  file.copy(plot_script, file.path(d, "plotting_script.R"), overwrite = TRUE)
  status <- system2("Rscript", c(shQuote(file.path(d, "plotting_script.R")), shQuote(d)))
  if (status != 0L) stop("Plotting failed for ", d)
}

sig <- !is.na(res$padj) & res$padj < 0.05
strict <- sig & abs(res$log2FoldChange_shrunk) >= 1
tested_n <- sum(!is.na(res$pvalue))
padj_na_n <- sum(is.na(res$padj))
up_n <- sum(sig & res$log2FoldChange_shrunk > 0)
down_n <- sum(sig & res$log2FoldChange_shrunk < 0)
strict_up <- sum(strict & res$log2FoldChange_shrunk > 0)
strict_down <- sum(strict & res$log2FoldChange_shrunk < 0)
top_down <- head(res[order(res$stat, na.last = NA), c("gene_name", "gene_id", "stat", "padj", "log2FoldChange_shrunk")], 10)
top_up <- head(res[order(-res$stat, na.last = NA), c("gene_name", "gene_id", "stat", "padj", "log2FoldChange_shrunk")], 10)
write.table(top_down, file.path(out_dir, "TOP10_EDGE_DOWN_BY_WALD_STAT.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(top_up, file.path(out_dir, "TOP10_EDGE_UP_BY_WALD_STAT.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

summary_lines <- c(
  "# DESeq2 Edge vs Core summary", "",
  paste0("- Total genes before filtering: ", nrow(counts_all)),
  paste0("- Genes retained by `rowSums(counts >= 10) >= 3`: ", nrow(counts_filtered)),
  paste0("- Genes removed: ", nrow(counts_all) - nrow(counts_filtered)),
  paste0("- Genes with a DESeq2 P value: ", tested_n),
  paste0("- Genes with adjusted P unavailable after independent filtering: ", padj_na_n),
  paste0("- Genes with adjusted P < 0.05: ", sum(sig)),
  paste0("- Edge-up at adjusted P < 0.05: ", up_n),
  paste0("- Edge-down at adjusted P < 0.05: ", down_n),
  paste0("- Edge-up at adjusted P < 0.05 and |shrunken LFC| >= 1: ", strict_up),
  paste0("- Edge-down at adjusted P < 0.05 and |shrunken LFC| >= 1: ", strict_down),
  paste0("- Heatmap selection: ", heatmap_rule)
)
writeLines(summary_lines, file.path(out_dir, "DESEQ2_SUMMARY.md"))

pde_natural <- if (!is.na(pde$rank_abs_Wald_stat_descending) && pde$rank_abs_Wald_stat_descending <= 30) "PDE2A naturally falls within the top 30 genes by absolute Wald statistic." else "PDE2A does not naturally fall within the top 30 genes by absolute Wald statistic."
audit_lines <- c(
  "# Checkpoint 5 paired differential-expression audit", "",
  paste0("Date: ", Sys.Date()), "",
  paste0("1. DESeq2 version: ", packageVersion("DESeq2")),
  paste0("2. apeglm version: ", packageVersion("apeglm")),
  "3. Input matrix: `results/04_featurecounts/count_matrix_gene_id.tsv` (raw non-negative integer featureCounts matrix)",
  "4. Analysis samples: 10; patients: 5; paracancerous samples: excluded from this model",
  "5. Model: `~ patient_id + region`",
  "6. Contrast: `region ti vs c`, permanently interpreted as Edge - Core",
  paste0("7. Low-count filter: 62,700 -> ", nrow(counts_filtered), " retained; frozen Checkpoint 4 count matched"),
  paste0("8. Adjusted P < 0.05: ", sum(sig), " total; ", up_n, " Edge-up; ", down_n, " Edge-down"),
  paste0("9. Adjusted P < 0.05 and |shrunken LFC| >= 1: ", strict_up, " Edge-up; ", strict_down, " Edge-down"),
  "10. PCA and sequencing-batch annotations were not included in the differential-expression model; no batch correction was applied.",
  "11. No sample was excluded.",
  "12. Filtering, FDR, effect-size, ranking, labeling, and heatmap rules were fixed before locating PDE2A.",
  paste0("13. PDE2A: raw LFC=", signif(pde$log2FoldChange_raw, 6), ", shrunken LFC=", signif(pde$log2FoldChange_shrunk, 6), ", Wald stat=", signif(pde$stat, 6), ", P=", signif(pde$pvalue, 6), ", adjusted P=", signif(pde$padj, 6), ", Edge-lower consistency=", pde$direction_consistency_edge_lower, "/5."),
  paste0("14. PDE2A ranks: adjusted-P rank=", pde$rank_padj_ascending, "; absolute-Wald rank=", pde$rank_abs_Wald_stat_descending, "; Edge-down rank=", pde$edge_down_rank_among_negative_stat, "."),
  paste0("15. ", pde_natural, " It is nevertheless within the adjusted-P < 0.05 set."),
  paste0("16. Heatmap rule: ", heatmap_rule),
  paste0("17. Adjusted P was unavailable for ", padj_na_n, " genes after DESeq2 independent filtering; all ", tested_n, " retained genes had a raw P value."),
  "18. GSEA, GO, KEGG, Reactome, PDE-family analysis, and candidate selection were not run."
)
writeLines(audit_lines, file.path(out_dir, "QC_CHECKPOINT5_DEG_AUDIT.md"))

saveRDS(dds, file.path(out_dir, "DESeq2_paired_dds.rds"))
writeLines(c(paste0("DESeq2=", packageVersion("DESeq2")), paste0("apeglm=", packageVersion("apeglm")),
             paste0("coefficient=", coef_name), paste0("resultsNames=", paste(coef_names, collapse = ","))),
           file.path(out_dir, "SOFTWARE_AND_COEFFICIENTS.txt"))
