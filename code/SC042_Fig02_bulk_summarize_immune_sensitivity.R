base <- "results/12_gene_set_immune_sensitivity"
rd <- function(f) read.delim(file.path(base, f), check.names = FALSE,
                             stringsAsFactors = FALSE)
show_sig <- function(x, pcol, n = 12L) {
  cat("n=", nrow(x), " sig=", sum(x[[pcol]] < 0.05, na.rm = TRUE), "\n", sep = "")
  print(utils::head(x[order(x[[pcol]], na.last = TRUE), , drop = FALSE], n),
        row.names = FALSE)
}

cat("=== CAMERA ===\n")
x <- rd("CAMERA_EDGE_VS_CORE.tsv")
print(with(x, tapply(FDR, database, function(z) sum(z < 0.05, na.rm = TRUE))))
show_sig(x, "FDR")

cat("=== FRY ===\n")
x <- rd("FRY_EDGE_VS_CORE.tsv")
print(with(x, tapply(FDR, database, function(z) sum(z < 0.05, na.rm = TRUE))))
show_sig(x, "FDR")

cat("=== PADOG ===\n")
x <- rd("PADOG_PAIRED_EDGE_VS_CORE.tsv")
print(with(x, tapply(BH_FDR, database, function(z) sum(z < 0.05, na.rm = TRUE))))
show_sig(x, "BH_FDR")

cat("=== DIRAC ===\n")
show_sig(rd("DIRAC_HALLMARK_EXPLORATORY.tsv"), "BH_FDR")

cat("=== SSGSEA ===\n")
show_sig(rd("SSGSEA_HALLMARK_PAIRED_LIMMA.tsv"), "BH_FDR")

cat("=== MCP ===\n")
show_sig(rd("MPCOUNTER_PAIRED_LIMMA.tsv"), "BH_FDR", 20)

cat("=== ESTIMATE ===\n")
show_sig(rd("ESTIMATE_PAIRED_LIMMA.tsv"), "BH_FDR", 20)

cat("=== IPS ===\n")
show_sig(rd("IPS_PAIRED_LIMMA.tsv"), "BH_FDR", 20)

cat("=== TIDE ===\n")
print(rd("TIDE_EXECUTION_STATUS.tsv"), row.names = FALSE)

cat("=== IPS mapping ===\n")
print(table(rd("IPS_GENE_MAPPING_AUDIT.tsv")$status))

cat("=== CAMERA significant direction by database ===\n")
x <- rd("CAMERA_EDGE_VS_CORE.tsv")
print(with(subset(x, FDR < 0.05), table(database, Direction)))
cat("CAMERA top Edge-enriched\n")
print(head(subset(x, FDR < 0.05 & Direction == "Up")[order(subset(x, FDR < 0.05 & Direction == "Up")$FDR),
      c("database", "pathway", "NGenes", "FDR")], 15), row.names = FALSE)
cat("CAMERA top Core-enriched\n")
print(head(subset(x, FDR < 0.05 & Direction == "Down")[order(subset(x, FDR < 0.05 & Direction == "Down")$FDR),
      c("database", "pathway", "NGenes", "FDR")], 15), row.names = FALSE)

pat <- "cyclic|cAMP|cGMP|adenylate|guanylate|GPCR|G protein|MAPK|growth factor|vascular|endothelial|extracellular matrix|ECM"
cat("=== CAMERA selected signaling/context terms ===\n")
sel <- x[grepl(pat, x$pathway, ignore.case = TRUE), ]
print(head(sel[order(sel$FDR), c("database", "pathway", "Direction", "PValue", "FDR")], 40), row.names = FALSE)

cat("=== PADOG significant descriptive directions ===\n")
x <- rd("PADOG_PAIRED_EDGE_VS_CORE.tsv")
print(with(subset(x, BH_FDR < 0.05), table(database, descriptive_direction)))
cat("=== PADOG selected signaling/context terms ===\n")
sel <- x[grepl(pat, x$pathway, ignore.case = TRUE), ]
print(head(sel[order(sel$BH_FDR), c("database", "pathway", "Ppadog", "BH_FDR",
                                   "mean_member_gene_delta", "descriptive_direction")], 40), row.names = FALSE)
