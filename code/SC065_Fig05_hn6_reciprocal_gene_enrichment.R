suppressPackageStartupMessages({library(clusterProfiler); library(org.Hs.eg.db); library(readr); library(dplyr)})
infile <- '/data/liaolz/PDE2A/results/reversible_signature/PDE2A_KD_OE_reversible_signature.tsv'
out <- '/data/liaolz/PDE2A/results/PDE2A_reversible_enrichment'
dir.create(out, showWarnings=FALSE, recursive=TRUE)
x <- read_tsv(infile, show_col_types=FALSE)
for (cls in c('loss_up','loss_down')) {
  genes <- unique(x$gene[x$loss_class==cls])
  eg <- bitr(genes, fromType='SYMBOL', toType='ENTREZID', OrgDb=org.Hs.eg.db)
  if (!nrow(eg)) next
  e <- enrichGO(eg$ENTREZID, OrgDb=org.Hs.eg.db, ont='BP', pAdjustMethod='BH', pvalueCutoff=0.2, qvalueCutoff=0.5, readable=TRUE)
  write_tsv(as.data.frame(e), file.path(out, paste0('PDE2A_', cls, '_GO_BP.tsv')))
}
writeLines(c('Input: PDE2A KD/OE reversible signature','Method: clusterProfiler enrichGO, BP, BH correction','Threshold: p<0.2, q<0.5'), file.path(out,'enrichment_metadata.txt'))
