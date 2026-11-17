#!/usr/bin/env Rscript

# PDE2A pre-analysis provenance audit.
# Deliberately limited to annotation lookup, DESeq2 size-factor normalization,
# patient-wise direction summaries, and one internal development PNG.

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
out_dir <- file.path(project_root, "results", "05_pde2a_audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

annotation_file <- file.path(project_root, "results", "04_featurecounts", "gene_annotation.tsv")
count_file <- file.path(project_root, "results", "04_featurecounts", "count_matrix_gene_id.tsv")
manifest_file <- file.path(project_root, "manifest_combined_standardized.tsv")
for (path in c(annotation_file, count_file, manifest_file)) {
  if (!file.exists(path) || file.info(path)$size == 0) stop("Missing or empty required input: ", path)
}

annotation <- read.delim(annotation_file, check.names = FALSE, stringsAsFactors = FALSE)
pde2a_annotation <- annotation[!is.na(annotation$gene_name) & annotation$gene_name == "PDE2A", , drop = FALSE]
if (nrow(pde2a_annotation) != 1L) stop("PDE2A annotation is not unique: found ", nrow(pde2a_annotation), " entries")
write.table(pde2a_annotation[, c("gene_id", "gene_id_stripped", "gene_name", "gene_type")],
            file.path(out_dir, "PDE2A_annotation_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

counts_all <- read.delim(count_file, check.names = FALSE, stringsAsFactors = FALSE)
if (anyDuplicated(counts_all$gene_id)) stop("Duplicate gene_id values in frozen count matrix")
row <- counts_all[counts_all$gene_id == pde2a_annotation$gene_id, , drop = FALSE]
if (nrow(row) != 1L) stop("Resolved PDE2A gene_id does not map uniquely to frozen count matrix")

sample_order <- c("P08_C", "P08_TI", "P09_C", "P09_TI", "P34_C", "P34_TI", "P36_C", "P36_TI", "P38_C", "P38_TI")
if (!all(sample_order %in% names(counts_all))) stop("One or more frozen main-analysis sample columns are missing")
count_matrix <- matrix(as.integer(unlist(counts_all[, sample_order], use.names = FALSE)), nrow = nrow(counts_all),
                       dimnames = list(counts_all$gene_id, sample_order))
if (anyNA(count_matrix) || any(count_matrix < 0L)) stop("Frozen count matrix contains invalid counts")

manifest <- read.delim(manifest_file, check.names = FALSE, stringsAsFactors = FALSE)
meta <- manifest[match(sample_order, manifest$sample_id), c("sample_id", "patient_id", "region"), drop = FALSE]
if (anyNA(meta$sample_id) || !identical(meta$sample_id, sample_order)) stop("Manifest mapping failed for the frozen sample order")
meta$patient_id <- factor(meta$patient_id)
meta$region <- factor(meta$region, levels = c("c", "ti"))
rownames(meta) <- meta$sample_id

dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = meta, design = ~ patient_id + region)
dds <- estimateSizeFactors(dds)
normalized <- counts(dds, normalized = TRUE)
sf <- sizeFactors(dds)
gene_id <- pde2a_annotation$gene_id[[1]]

normalized_long <- data.frame(
  patient_id = as.character(meta$patient_id),
  region = as.character(meta$region),
  sample_id = sample_order,
  raw_count = as.integer(count_matrix[gene_id, sample_order]),
  size_factor = as.numeric(sf[sample_order]),
  normalized_count = as.numeric(normalized[gene_id, sample_order]),
  log2_normalized_count_plus1 = log2(as.numeric(normalized[gene_id, sample_order]) + 1),
  stringsAsFactors = FALSE
)
write.table(normalized_long, file.path(out_dir, "PDE2A_normalized_counts_preanalysis.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

patients <- c("P08", "P09", "P34", "P36", "P38")
direction_rows <- lapply(patients, function(patient) {
  core <- normalized_long[normalized_long$patient_id == patient & normalized_long$region == "c", , drop = FALSE]
  edge <- normalized_long[normalized_long$patient_id == patient & normalized_long$region == "ti", , drop = FALSE]
  if (nrow(core) != 1L || nrow(edge) != 1L) stop("Pairing failure for ", patient)
  data.frame(
    patient_id = patient,
    core_raw = core$raw_count,
    edge_raw = edge$raw_count,
    core_normalized = core$normalized_count,
    edge_normalized = edge$normalized_count,
    delta_normalized_edge_minus_core = edge$normalized_count - core$normalized_count,
    log2_ratio_edge_vs_core = log2((edge$normalized_count + 1) / (core$normalized_count + 1)),
    edge_lower_than_core = ifelse(edge$normalized_count < core$normalized_count, "YES", "NO"),
    stringsAsFactors = FALSE
  )
})
direction <- do.call(rbind, direction_rows)
write.table(direction, file.path(out_dir, "PDE2A_patientwise_direction_audit.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

raw_pairs <- data.frame(
  patient_id = direction$patient_id,
  core_sample = paste0(direction$patient_id, "_C"),
  edge_sample = paste0(direction$patient_id, "_TI"),
  core_raw_count = direction$core_raw,
  edge_raw_count = direction$edge_raw,
  raw_edge_minus_core = direction$edge_raw - direction$core_raw,
  raw_edge_div_core = ifelse(direction$core_raw == 0, NA_real_, direction$edge_raw / direction$core_raw),
  stringsAsFactors = FALSE
)
write.table(raw_pairs, file.path(out_dir, "PDE2A_raw_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

plot_data <- normalized_long
plot_data$region_label <- factor(ifelse(plot_data$region == "c", "Core", "Edge"), levels = c("Core", "Edge"))
write.table(plot_data, file.path(out_dir, "PDE2A_paired_preanalysis_internal_source_data.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

p <- ggplot(plot_data, aes(x = region_label, y = log2_normalized_count_plus1, group = patient_id, color = patient_id)) +
  geom_line(linewidth = 0.65, alpha = 0.8) +
  geom_point(size = 2.8) +
  geom_text(aes(label = patient_id), nudge_x = 0.055, hjust = 0, size = 3.2, show.legend = FALSE) +
  scale_color_manual(values = c(P08 = "#4E79A7", P09 = "#59A14F", P34 = "#F28E2B", P36 = "#B07AA1", P38 = "#E15759")) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "log2(DESeq2 normalized count + 1)", title = "PDE2A paired pre-analysis audit") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", plot.margin = margin(5.5, 34, 5.5, 5.5))
ggsave(file.path(out_dir, "PDE2A_paired_preanalysis_internal.png"), p, width = 5.4, height = 4.4, dpi = 300, bg = "white")

writeLines(c(
  "This script intentionally does not call DESeq(), results(), lfcShrink(), a Wald test, or any DEG procedure.",
  paste0("DESeq2 version: ", as.character(packageVersion("DESeq2"))),
  paste0("Edge lower than Core: ", sum(direction$edge_lower_than_core == "YES"), "/5")
), file.path(out_dir, "PREANALYSIS_EXECUTION_NOTE.txt"))
