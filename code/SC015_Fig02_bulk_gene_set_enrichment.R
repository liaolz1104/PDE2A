#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) stop("Usage: 08_enrichment.R <DESeq2_csv> <outdir>")
res_file <- args[1]
outdir <- args[2]
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

x <- read.csv(res_file, stringsAsFactors=FALSE)

# Threshold-based ORA is secondary in this n=5 paired cohort.
# If too few DEGs pass the threshold, the script skips rather than relaxing thresholds.
deg <- x[!is.na(x$padj) & x$padj < 0.05 & abs(x$log2FoldChange_shrunk) >= 1 & !is.na(x$gene_name), ]
up <- unique(deg$gene_name[deg$log2FoldChange_shrunk > 0])
down <- unique(deg$gene_name[deg$log2FoldChange_shrunk < 0])

run_go <- function(genes, prefix) {
  if (length(genes) < 10) {
    writeLines(
      paste0("Skipped: only ", length(genes), " genes passed the prespecified threshold."),
      file.path(outdir, paste0(prefix, "_GO_skipped.txt"))
    )
    return(invisible(NULL))
  }
  ego <- enrichGO(
    gene=genes,
    OrgDb=org.Hs.eg.db,
    keyType="SYMBOL",
    ont="BP",
    pAdjustMethod="BH",
    pvalueCutoff=0.05,
    qvalueCutoff=0.2,
    readable=TRUE
  )
  write.csv(as.data.frame(ego), file.path(outdir, paste0(prefix, "_GO_BP.csv")), row.names=FALSE)
}

run_go(up, "Edge_up")
run_go(down, "Edge_down")

# GSEA: preferred pathway-level analysis for this small paired cohort.
if (!requireNamespace("msigdbr", quietly=TRUE)) {
  writeLines(
    "Install msigdbr to run Hallmark GSEA: install.packages('msigdbr')",
    file.path(outdir, "GSEA_msigdbr_not_installed.txt")
  )
  quit(save="no", status=0)
}

rank <- x[!is.na(x$stat) & !is.na(x$gene_name) & x$gene_name != "", c("gene_name","stat")]
rank <- rank[order(abs(rank$stat), decreasing=TRUE), ]
rank <- rank[!duplicated(rank$gene_name), ]
geneList <- rank$stat
names(geneList) <- rank$gene_name
geneList <- sort(geneList, decreasing=TRUE)

get_msigdbr <- function(collection, subcollection=NULL) {
  f <- get("msigdbr", asNamespace("msigdbr"))
  fm <- names(formals(f))
  if ("collection" %in% fm) {
    if (is.null(subcollection)) {
      f(species="Homo sapiens", collection=collection)
    } else {
      f(species="Homo sapiens", collection=collection, subcollection=subcollection)
    }
  } else {
    if (is.null(subcollection)) {
      f(species="Homo sapiens", category=collection)
    } else {
      f(species="Homo sapiens", category=collection, subcategory=subcollection)
    }
  }
}

hall <- get_msigdbr("H")
term2gene <- unique(hall[, c("gs_name","gene_symbol")])
colnames(term2gene) <- c("term","gene")

g <- GSEA(
  geneList=geneList,
  TERM2GENE=term2gene,
  minGSSize=10,
  maxGSSize=500,
  pvalueCutoff=1,
  pAdjustMethod="BH",
  verbose=FALSE,
  seed=TRUE
)
write.csv(as.data.frame(g), file.path(outdir, "GSEA_Hallmark_all.csv"), row.names=FALSE)

message("Enrichment analysis finished.")
message("For Edge-Core orientation: NES > 0 = enriched in Edge; NES < 0 = enriched in Core.")
