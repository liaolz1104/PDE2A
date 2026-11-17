args <- commandArgs(trailingOnly = TRUE)
expr_file <- args[1]
out_dir <- args[2]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(DESeq2))

x <- read.delim(expr_file, check.names = FALSE, stringsAsFactors = FALSE)
count_cols <- grep('^read_count_', colnames(x), value = TRUE)
stopifnot(length(count_cols) == 12)
genes <- as.character(x$gene_symbol)
keep <- !is.na(genes) & genes != ''
x <- x[keep, , drop = FALSE]
genes <- genes[keep]
counts <- as.matrix(x[, count_cols, drop = FALSE])
mode(counts) <- 'numeric'
counts <- round(counts)
rownames(counts) <- make.unique(genes)
colnames(counts) <- sub('^read_count_', '', count_cols)
counts <- counts[rowSums(counts) > 0, , drop = FALSE]

condition <- ifelse(grepl('^HN6_KD', colnames(counts)), 'KD',
                    ifelse(grepl('^HN6_NC', colnames(counts)), 'NC',
                    ifelse(grepl('^HN6_OE', colnames(counts)), 'OE', 'Vector')))
meta <- data.frame(condition=factor(condition, levels=c('NC','KD','Vector','OE')))
rownames(meta) <- colnames(counts)

run_contrast <- function(numerator, denominator, prefix) {
  keep_gene <- rowSums(counts[, rownames(meta)[meta$condition %in% c(numerator, denominator)], drop=FALSE] >= 10) >= 3
  sub_counts <- counts[keep_gene, rownames(meta)[meta$condition %in% c(numerator, denominator)], drop=FALSE]
  sub_meta <- meta[colnames(sub_counts), , drop=FALSE]
  sub_meta$condition <- droplevels(sub_meta$condition)
  dds <- DESeqDataSetFromMatrix(countData=sub_counts, colData=sub_meta, design=~condition)
  dds <- DESeq(dds, quiet=TRUE)
  res <- results(dds, contrast=c('condition', numerator, denominator), independentFiltering=TRUE)
  tab <- as.data.frame(res)
  tab$gene <- rownames(tab)
  tab <- tab[, c('gene','log2FoldChange','pvalue','padj','stat')]
  colnames(tab) <- c('gene','logFC','P.Value','adj.P.Val','statistic')
  tab$direction <- ifelse(tab$logFC > 0, 'up', ifelse(tab$logFC < 0, 'down', 'neutral'))
  tab$contrast <- paste0('PDE2A-', numerator, ' vs ', denominator)
  write.table(tab, file.path(out_dir, paste0(prefix, '.tsv')), sep='\t', row.names=FALSE, quote=FALSE, na='NA')
}

run_contrast('KD', 'NC', 'PDE2A_KD_vs_Control_DEG_recomputed')
run_contrast('OE', 'Vector', 'PDE2A_OE_vs_Control_DEG_recomputed')
write.table(meta, file.path(out_dir, 'PDE2A_KD_OE_sample_metadata.tsv'), sep='\t', row.names=TRUE, quote=FALSE)
writeLines(c('count_source=vendor read_count expected counts', 'preprocessing=rounded to integers before DESeq2', 'design=~condition', 'replicates=3 per group'), file.path(out_dir, 'PDE2A_KD_OE_recomputed_metadata.txt'))
