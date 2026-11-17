suppressPackageStartupMessages({
  library(DESeq2)
  library(fgsea)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript run_deseq2_gsea.R OUT_DIR MANIFEST GENE_MAP GENE_SETS_GMT")
}
out_dir <- normalizePath(args[[1]], mustWork = TRUE)
manifest_file <- normalizePath(args[[2]], mustWork = TRUE)
gene_map_file <- normalizePath(args[[3]], mustWork = TRUE)
gmt_file <- normalizePath(args[[4]], mustWork = TRUE)
dir.create(file.path(out_dir, "deseq2"), showWarnings = FALSE)
de_dir <- file.path(out_dir, "deseq2")

manifest <- read.delim(manifest_file, check.names = FALSE, stringsAsFactors = FALSE)
gene_map <- unique(read.delim(gene_map_file, stringsAsFactors = FALSE))

read_fc <- function(strand) {
  path <- file.path(out_dir, "featurecounts", paste0("gene_counts_s", strand, ".txt"))
  tab <- read.delim(path, comment.char = "#", check.names = FALSE)
  counts <- as.matrix(tab[, 7:ncol(tab), drop = FALSE])
  storage.mode(counts) <- "integer"
  rownames(counts) <- tab$Geneid
  if (ncol(counts) != nrow(manifest)) stop("featureCounts column count does not match manifest")
  colnames(counts) <- manifest$sample
  counts
}

count_sets <- lapply(0:2, read_fc)
totals <- sapply(count_sets, sum)
names(totals) <- paste0("s", 0:2)
directional_ratio <- max(totals[["s1"]], totals[["s2"]]) / totals[["s0"]]
selected_strand <- if (directional_ratio >= 0.60) {
  if (totals[["s1"]] >= totals[["s2"]]) 1L else 2L
} else {
  0L
}
strand_qc <- data.frame(
  strand = 0:2,
  assigned_fragments = as.numeric(totals),
  selected = 0:2 == selected_strand,
  directional_to_unstranded_ratio = directional_ratio
)
write.table(strand_qc, file.path(de_dir, "strandedness_selection.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
counts <- count_sets[[selected_strand + 1L]]

run_contrast <- function(comparison, numerator, reference) {
  keep_samples <- manifest$comparison == comparison
  meta <- manifest[keep_samples, c("sample", "condition"), drop = FALSE]
  meta$condition <- relevel(factor(meta$condition), ref = reference)
  rownames(meta) <- meta$sample
  subcounts <- counts[, meta$sample, drop = FALSE]
  keep_genes <- rowSums(subcounts) >= 10
  dds <- DESeqDataSetFromMatrix(subcounts[keep_genes, , drop = FALSE], meta, ~ condition)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", numerator, reference))
  res_df <- data.frame(gene_id = rownames(res), as.data.frame(res), check.names = FALSE)
  res_df <- merge(res_df, gene_map, by = "gene_id", all.x = TRUE, sort = FALSE)
  res_df$comparison <- comparison
  write.table(res_df, file.path(de_dir, paste0(comparison, "_DESeq2.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)

  vst_mat <- assay(vst(dds, blind = TRUE))
  cors <- cor(vst_mat, method = "pearson")
  write.table(cors, file.path(de_dir, paste0(comparison, "_VST_correlations.tsv")), sep = "\t", quote = FALSE)
  pca <- plotPCA(vst(dds, blind = TRUE), intgroup = "condition", returnData = TRUE)
  percent_var <- round(100 * attr(pca, "percentVar"))
  p <- ggplot(pca, aes(PC1, PC2, color = condition, label = name)) +
    geom_point(size = 3) + geom_text(vjust = -0.7, show.legend = FALSE) +
    xlab(paste0("PC1: ", percent_var[1], "%")) + ylab(paste0("PC2: ", percent_var[2], "%")) +
    theme_classic(base_size = 11)
  ggsave(file.path(de_dir, paste0(comparison, "_PCA.pdf")), p, width = 5.2, height = 4.5)

  isolation <- lapply(meta$sample, function(sample) {
    peers <- meta$sample[meta$condition == meta[sample, "condition"] & meta$sample != sample]
    within <- cors[sample, peers]
    data.frame(
      sample = sample,
      condition = as.character(meta[sample, "condition"]),
      median_within_r = median(within),
      peer_r_1 = within[[1]],
      peer_r_2 = within[[2]],
      severe_isolation = FALSE
    )
  })
  isolation <- do.call(rbind, isolation)
  for (condition_name in unique(isolation$condition)) {
    idx <- which(isolation$condition == condition_name)
    for (i in idx) {
      peer_medians <- isolation$median_within_r[setdiff(idx, i)]
      isolation$severe_isolation[i] <- isolation$median_within_r[i] < 0.80 &&
        all(peer_medians - isolation$median_within_r[i] >= 0.10)
    }
  }
  write.table(isolation, file.path(de_dir, paste0(comparison, "_outlier_flags.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  list(dds = dds, result = res_df, isolation = isolation)
}

kd <- run_contrast("KD_vs_NC", "KD", "NC")
oe <- run_contrast("OE_vs_Vector", "OE", "Vector")

pde_rows <- function(x) {
  rows <- x$result[x$result$gene_name == "PDE2A", c("comparison", "gene_id", "gene_name", "baseMean", "log2FoldChange", "stat", "pvalue", "padj")]
  rows[order(rows$padj, -rows$baseMean, na.last = TRUE), , drop = FALSE]
}
pde <- rbind(pde_rows(kd), pde_rows(oe))
write.table(pde, file.path(de_dir, "PDE2A_perturbation_gate.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

best_pde <- function(rows, comparison) {
  z <- rows[rows$comparison == comparison & !is.na(rows$padj), , drop = FALSE]
  if (!nrow(z)) return(NULL)
  z[order(z$padj, -z$baseMean), , drop = FALSE][1, , drop = FALSE]
}
kd_pde <- best_pde(pde, "KD_vs_NC")
oe_pde <- best_pde(pde, "OE_vs_Vector")
pde_gate <- !is.null(kd_pde) && !is.null(oe_pde) &&
  kd_pde$padj < 0.05 && kd_pde$log2FoldChange < 0 &&
  oe_pde$padj < 0.05 && oe_pde$log2FoldChange > 0
outlier_gate <- !any(kd$isolation$severe_isolation) && !any(oe$isolation$severe_isolation)

decision_file <- file.path(de_dir, "bridge_decision.txt")
if (!pde_gate || !outlier_gate) {
  reasons <- c(
    if (!pde_gate) "STOP: PDE2A perturbation did not pass the prespecified direction-and-FDR gate.",
    if (!outlier_gate) "STOP: at least one sample met the prespecified severe-isolation rule."
  )
  writeLines(c(reasons, "GSEA was not run."), decision_file)
  quit(save = "no", status = 0)
}

merge_cols <- c("gene_id", "gene_name", "log2FoldChange", "padj", "stat")
rev <- merge(kd$result[, merge_cols], oe$result[, merge_cols], by = "gene_id", suffixes = c("_KD", "_OE"))
rev$gene_name <- ifelse(!is.na(rev$gene_name_KD), rev$gene_name_KD, rev$gene_name_OE)
rev$reversible <- !is.na(rev$padj_KD) & !is.na(rev$padj_OE) &
  rev$padj_KD < 0.05 & rev$padj_OE < 0.05 &
  abs(rev$log2FoldChange_KD) >= 0.5 & abs(rev$log2FoldChange_OE) >= 0.5 &
  sign(rev$log2FoldChange_KD) == -sign(rev$log2FoldChange_OE)
write.table(rev, file.path(de_dir, "KD_OE_reversible_genes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

collapse_rank <- function(res_df) {
  x <- res_df[!is.na(res_df$gene_name) & !is.na(res_df$stat), c("gene_name", "gene_id", "stat")]
  x <- x[order(x$gene_name, -abs(x$stat), x$gene_id), ]
  x <- x[!duplicated(x$gene_name), ]
  stats <- x$stat
  names(stats) <- x$gene_name
  stats <- stats[order(-stats, names(stats))]
  tied <- duplicated(stats) | duplicated(stats, fromLast = TRUE)
  for (value in unique(stats[tied])) {
    idx <- which(stats == value)
    stats[idx] <- stats[idx] + rev(seq_along(idx)) * 1e-12
  }
  stats[order(-stats, names(stats))]
}
pathways <- gmtPathways(gmt_file)
ranks <- list(KD_vs_NC = collapse_rank(kd$result), OE_vs_Vector = collapse_rank(oe$result))
set.seed(20260805)
gsea <- do.call(rbind, lapply(names(ranks), function(comparison) {
  z <- fgseaMultilevel(pathways = pathways, stats = ranks[[comparison]], minSize = 5, maxSize = Inf, eps = 0)
  z$comparison <- comparison
  z
}))
gsea <- as.data.frame(gsea)
gsea$padj_global <- p.adjust(gsea$pval, method = "BH")
gsea$leadingEdge_genes <- vapply(gsea$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
gsea$leadingEdge <- NULL
write.table(gsea, file.path(de_dir, "preset_GSEA.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

get_row <- function(comparison, pathway) gsea[gsea$comparison == comparison & gsea$pathway == pathway, , drop = FALSE]
kd_full <- get_row("KD_vs_NC", "EDGE_CSC_155")
oe_full <- get_row("OE_vs_Vector", "EDGE_CSC_155")
kd_net <- get_row("KD_vs_NC", "STRICT_NETWORK_15")
oe_net <- get_row("OE_vs_Vector", "STRICT_NETWORK_15")
kd_le <- strsplit(kd_full$leadingEdge_genes, ";", fixed = TRUE)[[1]]
oe_le <- strsplit(oe_full$leadingEdge_genes, ";", fixed = TRUE)[[1]]
le_overlap <- intersect(kd_le, oe_le)
le_opposite <- length(le_overlap) > 0 && all(sign(ranks$KD_vs_NC[le_overlap]) == -sign(ranks$OE_vs_Vector[le_overlap]))

bridge_pass <- kd_full$padj_global < 0.05 && oe_full$padj_global < 0.05 &&
  sign(kd_full$NES) == -sign(oe_full$NES) && le_opposite &&
  sign(kd_net$NES) == sign(kd_full$NES) && sign(oe_net$NES) == sign(oe_full$NES)
gate_table <- data.frame(
  criterion = c("PDE2A perturbation", "no severe outlier", "full set FDR both", "opposite full-set NES", "leading-edge concordance", "strict-network direction", "overall bridge"),
  passed = c(pde_gate, outlier_gate, kd_full$padj_global < 0.05 && oe_full$padj_global < 0.05,
             sign(kd_full$NES) == -sign(oe_full$NES), le_opposite,
             sign(kd_net$NES) == sign(kd_full$NES) && sign(oe_net$NES) == sign(oe_full$NES), bridge_pass),
  detail = c("KD down and OE up, both FDR < 0.05", "prespecified VST correlation rule",
             paste("global FDR:", signif(kd_full$padj_global, 3), signif(oe_full$padj_global, 3)),
             paste("NES:", signif(kd_full$NES, 3), signif(oe_full$NES, 3)),
             paste(length(le_overlap), "shared leading-edge genes; all opposite =", le_opposite),
             paste("network NES:", signif(kd_net$NES, 3), signif(oe_net$NES, 3)),
             if (bridge_pass) "PASS" else "STOP")
)
write.table(gate_table, file.path(de_dir, "bridge_gate_table.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(if (bridge_pass) "PASS: prespecified reversible Edge CSC bridge criteria met." else "STOP: prespecified reversible Edge CSC bridge criteria not met.", decision_file)
