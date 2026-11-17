#!/usr/bin/env Rscript

# Checkpoint 4 count-level QC only. No differential-expression result is computed.
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(dendextend)
  library(RColorBrewer)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: 12_checkpoint4_qc.R <count_matrix.tsv> <metadata.tsv> <assignment_audit.tsv> <outdir> <figure_root>")
}
count_file <- args[1]
metadata_file <- args[2]
assignment_file <- args[3]
outdir <- args[4]
figure_root <- args[5]
for (path in c(count_file, metadata_file, assignment_file)) {
  if (!file.exists(path) || file.info(path)$size == 0) stop("Missing/empty input: ", path)
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)

sample_order <- c("P08_C","P08_TI","P08_P","P09_C","P09_TI","P09_P",
                  "P34_C","P34_TI","P36_C","P36_TI","P38_C","P38_TI")
main10 <- c("P08_C","P08_TI","P09_C","P09_TI","P34_C","P34_TI",
            "P36_C","P36_TI","P38_C","P38_TI")
counts_df <- read.delim(count_file, check.names = FALSE, stringsAsFactors = FALSE)
metadata <- read.delim(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
assignment <- read.delim(assignment_file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(identical(colnames(counts_df), c("gene_id", sample_order)))
stopifnot(identical(metadata$sample_id, sample_order))
rownames(metadata) <- metadata$sample_id
counts <- as.matrix(counts_df[, sample_order])
rownames(counts) <- counts_df$gene_id
storage.mode(counts) <- "integer"
if (any(counts < 0L) || any(is.na(counts))) stop("Counts are not non-negative integers.")

write_figure <- function(stem, source_data, draw_fun, width = 9, height = 6) {
  section <- if (grepl("^Figure_QC_featureCounts_", stem)) "featureCounts" else "counts"
  fig_dir <- file.path(figure_root, section, stem)
  expected <- c("source_data.tsv", "plotting_script.R", paste0(stem, c(".pdf", ".png", ".tiff")))
  if (dir.exists(fig_dir) && all(file.exists(file.path(fig_dir, expected))) &&
      all(file.info(file.path(fig_dir, expected))$size > 0)) {
    message("Complete figure already exists; keeping it: ", stem)
    return(invisible(NULL))
  }
  if (dir.exists(fig_dir) && length(list.files(fig_dir)) > 0) {
    archived <- paste0(fig_dir, ".failed.", format(Sys.time(), "%Y%m%dT%H%M%S"))
    if (!file.rename(fig_dir, archived)) stop("Could not archive partial figure directory: ", fig_dir)
  }
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(source_data, file.path(fig_dir, "source_data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  file.copy(normalizePath(sys.frame(1)$ofile %||% "scripts/12_checkpoint4_qc.R", mustWork = TRUE),
            file.path(fig_dir, "plotting_script.R"), overwrite = FALSE)
  pdf(file.path(fig_dir, paste0(stem, ".pdf")), width = width, height = height, useDingbats = FALSE)
  draw_fun(); dev.off()
  png(file.path(fig_dir, paste0(stem, ".png")), width = width, height = height, units = "in", res = 300)
  draw_fun(); dev.off()
  tiff(file.path(fig_dir, paste0(stem, ".tiff")), width = width, height = height, units = "in", res = 600, compression = "lzw")
  draw_fun(); dev.off()
}
`%||%` <- function(x, y) if (is.null(x)) y else x

theme_pub <- theme_classic(base_size = 11) +
  theme(axis.text = element_text(color = "black"), plot.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
batch_colors <- c(old = "#0072B2", new = "#D55E00")
region_colors <- c(c = "#0072B2", ti = "#D55E00", p = "#009E73")
region_labels <- c(c = "Core", ti = "Edge", p = "Paracancerous")

# featureCounts figures
assign_source <- assignment[, c("sample_id","patient_id","region","batch","assigned","assigned_percent")]
p_assign <- ggplot(assign_source, aes(factor(sample_id, levels = sample_order), assigned_percent, fill = batch)) +
  geom_col(color = "black", linewidth = 0.25) +
  scale_fill_manual(values = batch_colors) +
  labs(title = "featureCounts assignment rate", x = "Sample", y = "Assigned fragments (%)", fill = "Sequencing set") +
  theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1))
write_figure("Figure_QC_featureCounts_01_AssignmentRate", assign_source, function() print(p_assign), 10, 5.8)

composition_cols <- c("assigned", "unassigned_no_features", "unassigned_multimapping",
                      "unassigned_unmapped", "unassigned_ambiguity", "unassigned_chimera",
                      "unassigned_mapping_quality", "other_unassigned_if_present")
comp <- assignment[, c("sample_id","patient_id","region","batch",composition_cols)]
comp_long <- reshape(comp, varying = composition_cols, v.names = "fragments", timevar = "category",
                     times = composition_cols, direction = "long")
comp_long$percent <- 100 * comp_long$fragments / assignment$total_fragments_considered[match(comp_long$sample_id, assignment$sample_id)]
comp_check <- aggregate(percent ~ sample_id, comp_long, sum)
if (any(abs(comp_check$percent - 100) > 0.01)) stop("Assignment composition does not sum to 100%.")
comp_long$category <- factor(comp_long$category, levels = composition_cols,
  labels = c("Assigned","Unassigned_NoFeatures","Unassigned_MultiMapping","Unassigned_Unmapped",
             "Unassigned_Ambiguity","Unassigned_Chimera","Unassigned_MappingQuality","Other unassigned"))
comp_palette <- c("#0072B2","#E69F00","#56B4E9","#999999","#CC79A7","#D55E00","#009E73","#F0E442")
p_comp <- ggplot(comp_long, aes(factor(sample_id, levels = sample_order), percent, fill = category)) +
  geom_col(width = 0.82) + scale_fill_manual(values = comp_palette) +
  labs(title = "featureCounts assignment composition", x = "Sample", y = "Fragments (%)", fill = "Category") +
  theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
write_figure("Figure_QC_featureCounts_02_AssignmentComposition", comp_long, function() print(p_comp), 11.5, 6.2)

library_source <- data.frame(sample_id = sample_order,
                             patient_id = metadata[sample_order,"patient_id"],
                             region = metadata[sample_order,"region"],
                             batch = metadata[sample_order,"batch"],
                             gene_level_library_size = colSums(counts))
p_lib <- ggplot(library_source, aes(factor(sample_id, levels = sample_order), gene_level_library_size / 1e6, fill = batch)) +
  geom_col(color = "black", linewidth = 0.25) + scale_fill_manual(values = batch_colors) +
  labs(title = "Gene-level library size", x = "Sample", y = "Sum of raw gene counts (millions)", fill = "Sequencing set") +
  theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1))
write_figure("Figure_QC_featureCounts_03_LibrarySize", library_source, function() print(p_lib), 10, 5.8)

# Fixed low-count filter for transformation/QC only.
keep <- rowSums(counts >= 10L) >= 3L
filter_audit <- data.frame(rule = "rowSums(counts >= 10) >= 3",
                           total_genes_before_filtering = nrow(counts),
                           genes_retained = sum(keep), genes_removed = sum(!keep))
write.table(filter_audit, file.path(outdir, "LOW_COUNT_FILTER_AUDIT.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
filtered_counts <- counts[keep, , drop = FALSE]

dds_qc <- DESeqDataSetFromMatrix(countData = filtered_counts, colData = metadata, design = ~1)
vsd_qc <- vst(dds_qc, blind = TRUE)
saveRDS(vsd_qc, file.path(outdir, "VST_all12_QC.rds"))
vst_mat <- assay(vsd_qc)

# Pearson correlation and hierarchical clustering using distance = 1 - Pearson correlation.
cor_mat <- cor(vst_mat, method = "pearson")
write.table(cbind(sample_id = rownames(cor_mat), as.data.frame(cor_mat, check.names = FALSE)),
            file.path(outdir, "sample_correlation_Pearson.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
ann <- metadata[sample_order, c("region","batch","patient_id"), drop = FALSE]
colnames(ann) <- c("Region","Batch","Patient")
ann$Region <- factor(ann$Region, levels = c("c","ti","p"), labels = c("Core","Edge","Paracancerous"))
ann$Batch <- factor(ann$Batch, levels = c("old","new"))
ann_colors <- list(Region = c(Core="#0072B2", Edge="#D55E00", Paracancerous="#009E73"),
                   Batch = c(old="#56B4E9", new="#E69F00"),
                   Patient = setNames(brewer.pal(5, "Set2"), c("P08","P09","P34","P36","P38")))
cor_source <- data.frame(sample_id = rownames(cor_mat), cor_mat, check.names = FALSE)
draw_cor <- function() pheatmap(cor_mat, cluster_rows = TRUE, cluster_cols = TRUE,
  annotation_row = ann, annotation_col = ann, annotation_colors = ann_colors,
  color = colorRampPalette(c("white", "#56B4E9", "#0072B2"))(100),
  border_color = NA, main = "Pearson correlation of VST expression")
write_figure("Figure_QC_Count_01_SampleCorrelation", cor_source, draw_cor, 10, 9)

dist_mat <- as.dist(1 - cor_mat)
hc <- hclust(dist_mat, method = "average")
distance_full <- as.matrix(dist_mat)
hc_source <- merge(data.frame(sample_id = sample_order, distance_full, check.names = FALSE),
                   metadata, by = "sample_id", sort = FALSE)
dummy <- matrix(seq_len(12), nrow = 12, ncol = 1, dimnames = list(sample_order, ""))
draw_hc <- function() pheatmap(dummy, cluster_rows = hc, cluster_cols = FALSE,
  annotation_row = ann, annotation_colors = ann_colors, show_colnames = FALSE,
  color = rep("white", 100), breaks = seq(0.5, 12.5, length.out = 101),
  legend = FALSE, border_color = NA, cellwidth = 2, treeheight_row = 100,
  main = "Hierarchical clustering: 1 - Pearson correlation")
write_figure("Figure_QC_Count_02_HierarchicalClustering", hc_source, draw_hc, 6.5, 7)

# PCA on the all-12 blind VST matrix.
pca_all <- prcomp(t(vst_mat), center = TRUE, scale. = FALSE)
var_all <- 100 * pca_all$sdev^2 / sum(pca_all$sdev^2)
pca_coords <- data.frame(sample_id = rownames(pca_all$x), PC1 = pca_all$x[,1], PC2 = pca_all$x[,2])
pca_coords <- merge(pca_coords, metadata, by = "sample_id", sort = FALSE)
pca_coords <- pca_coords[match(sample_order, pca_coords$sample_id), ]
pca_coords$PC1_variance_percent <- var_all[1]
pca_coords$PC2_variance_percent <- var_all[2]
write.table(pca_coords, file.path(outdir, "PCA_all12_coordinates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

p_region <- ggplot(pca_coords, aes(PC1, PC2, color = region, shape = batch, label = sample_id)) +
  geom_point(size = 3) + geom_text_repel(size = 3, max.overlaps = Inf) +
  scale_color_manual(values = region_colors, labels = region_labels) +
  labs(title = "PCA of all 12 samples by region", x = sprintf("PC1 (%.1f%%)", var_all[1]),
       y = sprintf("PC2 (%.1f%%)", var_all[2]), color = "Region", shape = "Sequencing set") + theme_pub
write_figure("Figure_QC_Count_03_PCA_Region", pca_coords, function() print(p_region), 8.5, 6.5)

p_batch <- ggplot(pca_coords, aes(PC1, PC2, color = batch, shape = region, label = sample_id)) +
  geom_point(size = 3) + geom_text_repel(size = 3, max.overlaps = Inf) +
  scale_color_manual(values = batch_colors) +
  scale_shape_manual(values = c(c=16, ti=17, p=15), labels = region_labels) +
  labs(title = "PCA by sequencing batch annotation", x = sprintf("PC1 (%.1f%%)", var_all[1]),
       y = sprintf("PC2 (%.1f%%)", var_all[2]), color = "Sequencing set", shape = "Region") + theme_pub
write_figure("Figure_QC_Count_04_PCA_Batch", pca_coords, function() print(p_batch), 8.5, 6.5)

# Paired PCA uses the same VST matrix, restricted to the frozen 10 Core/Edge samples.
pca_main <- prcomp(t(vst_mat[, main10, drop = FALSE]), center = TRUE, scale. = FALSE)
var_main <- 100 * pca_main$sdev^2 / sum(pca_main$sdev^2)
paired <- data.frame(sample_id = rownames(pca_main$x), PC1 = pca_main$x[,1], PC2 = pca_main$x[,2])
paired <- merge(paired, metadata, by = "sample_id", sort = FALSE)
paired <- paired[match(main10, paired$sample_id), ]
paired$PC1_variance_percent <- var_main[1]
paired$PC2_variance_percent <- var_main[2]
write.table(paired, file.path(outdir, "PCA_main10_paired_coordinates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
p_pair <- ggplot(paired, aes(PC1, PC2, color = region, shape = batch, label = patient_id, group = patient_id)) +
  geom_line(color = "grey55", linewidth = 0.45) + geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps = Inf) +
  scale_color_manual(values = region_colors[c("c","ti")], labels = region_labels[c("c","ti")]) +
  labs(title = "Paired PCA of Core and Edge samples", x = sprintf("PC1 (%.1f%%)", var_main[1]),
       y = sprintf("PC2 (%.1f%%)", var_main[2]), color = "Region", shape = "Sequencing set") + theme_pub
write_figure("Figure_QC_Count_05_PairedPCA_CoreEdge", paired, function() print(p_pair), 8.5, 6.5)

# Objective descriptors for the final audit; no batch correction is performed.
median_corr <- sapply(sample_order, function(s) median(cor_mat[s, setdiff(sample_order, s)]))
nearest_distance <- apply(distance_full + diag(Inf, 12), 1, min)
robust_flags <- function(x, direction = "high") {
  med <- median(x); mad0 <- median(abs(x - med))
  z <- if (mad0 == 0) rep(0, length(x)) else 0.6745 * (x - med) / mad0
  if (direction == "high") abs(z) > 3.5 & z > 0 else abs(z) > 3.5 & z < 0
}
corr_outlier <- robust_flags(median_corr, "low")
cluster_outlier <- robust_flags(nearest_distance, "high")
batch_pc1_ranges <- tapply(pca_coords$PC1, pca_coords$batch, range)
batch_pc2_ranges <- tapply(pca_coords$PC2, pca_coords$batch, range)
range_overlap <- function(a, b) max(a[1], b[1]) <= min(a[2], b[2])
metrics <- data.frame(
  metric = c("genes_before_filter","genes_retained","genes_removed","PC1_variance_percent","PC2_variance_percent",
             "correlation_outlier_count","clustering_outlier_count","old_new_PC1_ranges_overlap","old_new_PC2_ranges_overlap",
             "batch_correction_applied"),
  value = c(nrow(counts), sum(keep), sum(!keep), var_all[1], var_all[2], sum(corr_outlier), sum(cluster_outlier),
            range_overlap(batch_pc1_ranges[["old"]], batch_pc1_ranges[["new"]]),
            range_overlap(batch_pc2_ranges[["old"]], batch_pc2_ranges[["new"]]), FALSE)
)
write.table(metrics, file.path(outdir, "COUNT_QC_METRICS.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
outliers <- data.frame(sample_id = sample_order, median_Pearson_to_others = median_corr,
                       nearest_neighbor_distance_1_minus_Pearson = nearest_distance,
                       correlation_outlier = corr_outlier, clustering_outlier = cluster_outlier)
write.table(outliers, file.path(outdir, "COUNT_QC_SAMPLE_OUTLIER_AUDIT.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Checkpoint 4 count-level QC completed without differential-expression testing.\n")
