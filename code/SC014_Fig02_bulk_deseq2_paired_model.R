#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(apeglm)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
  stop("Usage: 07_deseq2_paired_analysis.R <counts.tsv> <manifest.tsv> <gene_map.tsv> <outdir>")
}
count_file <- args[1]
manifest_file <- args[2]
gene_map_file <- args[3]
outdir <- args[4]
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

meta <- read.delim(manifest_file, check.names=FALSE, stringsAsFactors=FALSE)
counts <- read.delim(count_file, check.names=FALSE, stringsAsFactors=FALSE)
gene_map <- read.delim(gene_map_file, check.names=FALSE, stringsAsFactors=FALSE)

stopifnot("Geneid" %in% colnames(counts))
rownames(counts) <- counts$Geneid
counts$Geneid <- NULL

if (!all(meta$sample_id %in% colnames(counts))) {
  stop("Some manifest sample_id values are absent from count matrix.")
}
counts <- counts[, meta$sample_id, drop=FALSE]
counts <- as.matrix(counts)
storage.mode(counts) <- "integer"

rownames(meta) <- meta$sample_id
meta$patient_id <- factor(meta$patient_id)
meta$region <- factor(meta$region, levels=c("c","ti","p"))
meta$batch <- factor(meta$batch)

# ----------------------
# 1) All-12 QC / PCA
# ----------------------
keep_all <- rowSums(counts >= 10) >= 3
dds_all <- DESeqDataSetFromMatrix(
  countData = counts[keep_all, , drop=FALSE],
  colData = meta,
  design = ~ 1
)
vsd_all <- vst(dds_all, blind=TRUE)
saveRDS(vsd_all, file.path(outdir, "VST_all12.rds"))

pc <- prcomp(t(assay(vsd_all)), scale.=FALSE)
pct <- 100 * (pc$sdev^2 / sum(pc$sdev^2))
pca <- data.frame(
  sample_id=rownames(pc$x),
  PC1=pc$x[,1],
  PC2=pc$x[,2],
  meta[rownames(pc$x), c("patient_id","region","batch"), drop=FALSE]
)
write.csv(pca, file.path(outdir, "PCA_all12_coordinates.csv"), row.names=FALSE)

p <- ggplot(pca, aes(PC1, PC2, color=region, shape=batch, label=sample_id)) +
  geom_point(size=3.2) +
  geom_text_repel(size=3, max.overlaps=Inf) +
  theme_classic(base_size=13) +
  labs(
    x=sprintf("PC1 (%.1f%%)", pct[1]),
    y=sprintf("PC2 (%.1f%%)", pct[2]),
    title="All 12 samples: VST PCA"
  )
ggsave(file.path(outdir, "PCA_all12.pdf"), p, width=7.2, height=5.8)
ggsave(file.path(outdir, "PCA_all12.png"), p, width=7.2, height=5.8, dpi=300)

cors <- cor(assay(vsd_all), method="pearson")
pdf(file.path(outdir, "sample_correlation_all12.pdf"), width=9, height=8)
pheatmap(cors, main="Sample correlation (all 12)")
dev.off()

# ----------------------
# 2) Main paired C/TI analysis: 5 patients, 10 samples
# ----------------------
meta_main <- meta[meta$region %in% c("c","ti"), , drop=FALSE]
meta_main$region <- droplevels(meta_main$region)
meta_main$region <- relevel(meta_main$region, ref="c")
meta_main$patient_id <- droplevels(meta_main$patient_id)
counts_main <- counts[, rownames(meta_main), drop=FALSE]

pairtab <- table(meta_main$patient_id, meta_main$region)
if (nrow(pairtab) != 5 || any(pairtab[, "c"] != 1) || any(pairtab[, "ti"] != 1)) {
  stop("Main dataset is not exactly five complete Core/Edge pairs.")
}

keep <- rowSums(counts_main >= 10) >= 3
dds <- DESeqDataSetFromMatrix(
  countData = counts_main[keep, , drop=FALSE],
  colData = meta_main,
  design = ~ patient_id + region
)
dds <- DESeq(dds)

res_raw <- results(dds, contrast=c("region","ti","c"), alpha=0.05)
coef_name <- grep("^region_ti_vs_c$", resultsNames(dds), value=TRUE)
if (length(coef_name) != 1) {
  stop(paste("Could not find coefficient region_ti_vs_c. Results names:", paste(resultsNames(dds), collapse=", ")))
}
res_shr <- lfcShrink(dds, coef=coef_name, type="apeglm")

rr <- as.data.frame(res_raw)
rs <- as.data.frame(res_shr)
rr$gene_id <- rownames(rr)
rr$log2FoldChange_shrunk <- rs[rownames(rr), "log2FoldChange"]

gm <- gene_map
gm$gene_id <- sub("\\.\\d+$", "", gm$gene_id)
rr$gene_id <- sub("\\.\\d+$", "", rr$gene_id)
rr <- merge(rr, gm, by="gene_id", all.x=TRUE, sort=FALSE)
rr <- rr[order(rr$pvalue), ]
write.csv(rr, file.path(outdir, "DESeq2_Edge_vs_Core_all_genes.csv"), row.names=FALSE)

# Main 10-sample VST and paired PCA
vsd <- vst(dds, blind=FALSE)
saveRDS(dds, file.path(outdir, "dds_main_5pairs.rds"))
saveRDS(vsd, file.path(outdir, "vst_main_5pairs.rds"))

pc2 <- prcomp(t(assay(vsd)), scale.=FALSE)
pct2 <- 100 * (pc2$sdev^2 / sum(pc2$sdev^2))
pca2 <- data.frame(
  sample_id=rownames(pc2$x),
  PC1=pc2$x[,1],
  PC2=pc2$x[,2],
  meta_main[rownames(pc2$x), c("patient_id","region","batch"), drop=FALSE]
)
write.csv(pca2, file.path(outdir, "PCA_main10_coordinates.csv"), row.names=FALSE)

p2 <- ggplot(pca2, aes(PC1, PC2, color=region, shape=batch)) +
  geom_line(aes(group=patient_id), color="grey70", linewidth=0.5) +
  geom_point(size=3.5) +
  geom_text_repel(aes(label=sample_id), size=3, max.overlaps=Inf) +
  theme_classic(base_size=13) +
  labs(
    x=sprintf("PC1 (%.1f%%)", pct2[1]),
    y=sprintf("PC2 (%.1f%%)", pct2[2]),
    title="Five paired patients: Core vs Edge"
  )
ggsave(file.path(outdir, "PCA_main10_paired.pdf"), p2, width=7.2, height=5.8)
ggsave(file.path(outdir, "PCA_main10_paired.png"), p2, width=7.2, height=5.8, dpi=300)

# ----------------------
# 3) PDE2A prespecified target analysis
# ----------------------
pde_rows <- rr[!is.na(rr$gene_name) & rr$gene_name == "PDE2A", , drop=FALSE]
if (nrow(pde_rows) < 1) {
  stop("PDE2A not found in gene annotation/results. Check GTF gene_name mapping.")
}
pde_rows <- pde_rows[order(pde_rows$pvalue), , drop=FALSE]
write.csv(pde_rows, file.path(outdir, "PDE2A_DESeq2_result.csv"), row.names=FALSE)

pde_gid <- pde_rows$gene_id[1]
norm <- counts(dds, normalized=TRUE)
if (!(pde_gid %in% rownames(norm))) {
  # Ensembl version stripping fallback
  gid_match <- rownames(norm)[sub("\\.\\d+$","",rownames(norm)) == pde_gid]
  if (length(gid_match) < 1) stop("PDE2A gene_id not found in normalized counts.")
  pde_gid_use <- gid_match[1]
} else {
  pde_gid_use <- pde_gid
}

pde <- data.frame(
  sample_id=colnames(norm),
  normalized_count=as.numeric(norm[pde_gid_use, colnames(norm)]),
  stringsAsFactors=FALSE
)
pde <- cbind(pde, meta_main[pde$sample_id, c("patient_id","region","batch"), drop=FALSE])
pde$log2_norm_plus1 <- log2(pde$normalized_count + 1)
write.csv(pde, file.path(outdir, "PDE2A_normalized_counts_main10.csv"), row.names=FALSE)

wide <- reshape(
  pde[, c("patient_id","region","normalized_count","log2_norm_plus1")],
  idvar="patient_id", timevar="region", direction="wide"
)
wide$delta_log2_Edge_minus_Core <- wide$log2_norm_plus1.ti - wide$log2_norm_plus1.c
wide$edge_lower_than_core <- wide$delta_log2_Edge_minus_Core < 0
write.csv(wide, file.path(outdir, "PDE2A_patientwise_direction.csv"), row.names=FALSE)

# Descriptive paired tests on normalized expression.
wil <- wilcox.test(
  wide$log2_norm_plus1.ti,
  wide$log2_norm_plus1.c,
  paired=TRUE,
  exact=TRUE,
  alternative="two.sided"
)
direction_n <- sum(wide$edge_lower_than_core)
sign_two <- binom.test(direction_n, nrow(wide), p=0.5, alternative="two.sided")
sign_one <- binom.test(direction_n, nrow(wide), p=0.5, alternative="greater")

stats_txt <- c(
  paste0("Primary DESeq2 paired model: design = ~ patient_id + region"),
  paste0("Contrast: Edge(TI) - Core(C)"),
  paste0("PDE2A log2FC raw = ", signif(pde_rows$log2FoldChange[1], 5)),
  paste0("PDE2A log2FC shrunk = ", signif(pde_rows$log2FoldChange_shrunk[1], 5)),
  paste0("PDE2A Wald P = ", signif(pde_rows$pvalue[1], 5)),
  paste0("PDE2A BH-adjusted P = ", signif(pde_rows$padj[1], 5)),
  paste0("Direction consistency Edge<Core = ", direction_n, "/", nrow(wide)),
  paste0("Paired exact Wilcoxon two-sided P = ", signif(wil$p.value, 5)),
  paste0("Exact sign test two-sided P = ", signif(sign_two$p.value, 5)),
  paste0("Exact sign test one-sided P (prespecified Edge<Core direction only) = ", signif(sign_one$p.value, 5)),
  "Interpret the one-sided sign test only as a secondary direction-consistency test, not as the primary differential-expression test."
)
writeLines(stats_txt, file.path(outdir, "PDE2A_summary.txt"))

pp <- ggplot(pde, aes(x=region, y=log2_norm_plus1, group=patient_id)) +
  geom_line(alpha=0.65, linewidth=0.7) +
  geom_point(size=3.2) +
  theme_classic(base_size=13) +
  scale_x_discrete(labels=c(c="Core", ti="Edge")) +
  labs(
    x=NULL,
    y="PDE2A log2(normalized count + 1)",
    title="PDE2A expression: paired Core vs Edge"
  )
ggsave(file.path(outdir, "PDE2A_paired_Core_vs_Edge.pdf"), pp, width=5.2, height=5.5)
ggsave(file.path(outdir, "PDE2A_paired_Core_vs_Edge.png"), pp, width=5.2, height=5.5, dpi=300)

# ----------------------
# 4) Volcano plot
# ----------------------
vr <- rr
vr$neglog10_padj <- -log10(pmax(vr$padj, .Machine$double.xmin))
vr$sig <- !is.na(vr$padj) & vr$padj < 0.05 & abs(vr$log2FoldChange_shrunk) >= 1
vr$label <- ifelse(!is.na(vr$gene_name) & vr$gene_name == "PDE2A", "PDE2A", "")

pv <- ggplot(vr, aes(log2FoldChange_shrunk, neglog10_padj)) +
  geom_point(aes(alpha=sig), size=1.3) +
  scale_alpha_manual(values=c(`FALSE`=0.35, `TRUE`=0.9), guide="none") +
  geom_vline(xintercept=c(-1,1), linetype=2, linewidth=0.4) +
  geom_hline(yintercept=-log10(0.05), linetype=2, linewidth=0.4) +
  geom_text_repel(
    data=vr[vr$label != "", ],
    aes(label=label),
    size=4, min.segment.length=0
  ) +
  theme_classic(base_size=13) +
  labs(
    x="Shrunken log2FC (Edge - Core)",
    y="-log10(BH-adjusted P)",
    title="Paired transcriptome: Edge vs Core"
  )
ggsave(file.path(outdir, "Volcano_Edge_vs_Core.pdf"), pv, width=7, height=6)
ggsave(file.path(outdir, "Volcano_Edge_vs_Core.png"), pv, width=7, height=6, dpi=300)

# ----------------------
# 5) Heatmap: top DE genes + PDE2A
# ----------------------
ok <- rr[!is.na(rr$padj), , drop=FALSE]
top <- head(ok$gene_id[order(ok$padj)], 30)
top <- unique(c(top, pde_gid))
hm_ids <- rownames(vsd)[sub("\\.\\d+$","",rownames(vsd)) %in% top]
hm <- assay(vsd)[hm_ids, , drop=FALSE]
hm <- t(scale(t(hm)))
hm[is.na(hm)] <- 0

lookup <- gm
rownames(lookup) <- lookup$gene_id
hm_symbols <- lookup[sub("\\.\\d+$","",rownames(hm)), "gene_name"]
rownames(hm) <- ifelse(is.na(hm_symbols) | hm_symbols=="", rownames(hm), hm_symbols)

ann <- data.frame(
  Patient=meta_main[colnames(hm), "patient_id"],
  Region=meta_main[colnames(hm), "region"],
  Batch=meta_main[colnames(hm), "batch"]
)
rownames(ann) <- colnames(hm)

pdf(file.path(outdir, "Heatmap_top30_plus_PDE2A.pdf"), width=9, height=10)
pheatmap(hm, annotation_col=ann, cluster_cols=TRUE,
         show_colnames=TRUE, fontsize_row=7,
         main="Top paired DE genes + PDE2A")
dev.off()

# ----------------------
# 6) Ranked list for GSEA
# ----------------------
rank <- rr[!is.na(rr$stat) & !is.na(rr$gene_name) & rr$gene_name != "", c("gene_name","stat")]
rank <- rank[order(abs(rank$stat), decreasing=TRUE), ]
rank <- rank[!duplicated(rank$gene_name), ]
rank <- rank[order(rank$stat, decreasing=TRUE), ]
write.table(rank, file.path(outdir, "GSEA_rank_geneSymbol_WaldStat.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# ----------------------
# 7) P8/P9 paracancerous descriptive PDE2A plot
# ----------------------
if (all(c("P08_P","P09_P") %in% meta$sample_id)) {
  norm_all <- counts(dds_all, normalized=TRUE)
  pde_all_gid <- rownames(norm_all)[sub("\\.\\d+$","",rownames(norm_all)) == pde_gid]
  if (length(pde_all_gid) >= 1) {
    tmp <- data.frame(
      sample_id=colnames(norm_all),
      normalized_count=as.numeric(norm_all[pde_all_gid[1], ]),
      stringsAsFactors=FALSE
    )
    tmp <- cbind(tmp, meta[tmp$sample_id, c("patient_id","region","batch"), drop=FALSE])
    tmp <- tmp[tmp$patient_id %in% c("P08","P09"), ]
    tmp$log2_norm_plus1 <- log2(tmp$normalized_count + 1)
    write.csv(tmp, file.path(outdir, "PDE2A_P8_P9_Core_Edge_Paracancerous_descriptive.csv"), row.names=FALSE)
    ps <- ggplot(tmp, aes(region, log2_norm_plus1, group=patient_id)) +
      geom_line(alpha=0.65) +
      geom_point(size=3) +
      theme_classic(base_size=13) +
      scale_x_discrete(labels=c(c="Core", ti="Edge", p="Paracancerous")) +
      labs(x=NULL, y="PDE2A log2(normalized count + 1)",
           title="P8/P9 descriptive three-region PDE2A")
    ggsave(file.path(outdir, "PDE2A_P8_P9_three_region_descriptive.pdf"), ps, width=5.8, height=5.2)
  }
}

message("DESeq2 paired analysis finished.")
