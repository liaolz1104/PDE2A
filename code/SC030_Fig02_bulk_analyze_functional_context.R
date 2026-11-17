#!/usr/bin/env Rscript

# Checkpoint 6: transcriptome-wide enrichment and transparent candidate context.
# Uses frozen Checkpoint 5 results only; does not rerun DESeq2 or calculate a composite score.

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
local_lib <- file.path(project_root, ".Rlib_checkpoint5")
if (dir.exists(local_lib)) .libPaths(c(normalizePath(local_lib), .libPaths()))

suppressPackageStartupMessages({
  library(fgsea)
  library(ggplot2)
  library(ggrepel)
})

de_file <- file.path(project_root, "results", "05_deseq2", "DESeq2_Edge_vs_Core_all_genes.tsv")
gene_set_dir <- file.path(project_root, "results", "06_enrichment", "gene_sets")
plot_script <- file.path(project_root, "scripts", "19_plot_checkpoint6_enrichment.R")
out_dir <- file.path(project_root, "results", "06_enrichment")
fig_root <- file.path(project_root, "results", "07_figures", "manuscript_candidates", "Enrichment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_root, recursive = TRUE, showWarnings = FALSE)

libraries <- c(
  Hallmark = "MSigDB_Hallmark_2020.gmt",
  Reactome = "Reactome_2022.gmt",
  GO_BP = "GO_Biological_Process_2023.gmt",
  KEGG = "KEGG_2021_Human.gmt"
)
for (path in c(de_file, plot_script, file.path(gene_set_dir, libraries))) {
  if (!file.exists(path) || file.info(path)$size == 0) stop("Missing or empty frozen input: ", path)
}

de <- read.delim(de_file, check.names = FALSE, stringsAsFactors = FALSE, na.strings = "NA")
if (nrow(de) != 20413L || sum(!is.na(de$padj) & de$padj < 0.05) != 84L) stop("Frozen Checkpoint 5 result gate failed")
if (sum(!is.na(de$padj) & de$padj < 0.05 & de$log2FoldChange_shrunk < 0) != 66L) stop("Frozen Edge-down DEG count changed")
if (sum(!is.na(de$padj) & de$padj < 0.05 & de$log2FoldChange_shrunk > 0) != 18L) stop("Frozen Edge-up DEG count changed")

# Resolve duplicate symbols with one fixed rule for all libraries: greatest |Wald statistic|,
# then lowest raw P value, then lexical gene_id. No pathway or candidate is consulted.
rank_input <- de[!is.na(de$stat) & !is.na(de$gene_name) & nzchar(de$gene_name), ]
rank_input$gene_symbol <- toupper(rank_input$gene_name)
rank_input <- rank_input[order(-abs(rank_input$stat), rank_input$pvalue, rank_input$gene_id), ]
rank_input <- rank_input[!duplicated(rank_input$gene_symbol), ]
ranks <- rank_input$stat
names(ranks) <- rank_input$gene_symbol
ranks <- sort(ranks, decreasing = TRUE)
if (anyDuplicated(names(ranks)) || anyNA(ranks)) stop("Invalid GSEA ranking vector")
write.table(data.frame(gene = names(ranks), Wald_stat = as.numeric(ranks)),
            file.path(out_dir, "GSEA_WALD_STAT_RANK.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pathways <- lapply(libraries, function(filename) {
  x <- gmtPathways(file.path(gene_set_dir, filename))
  lapply(x, function(g) unique(toupper(g[nzchar(g)])))
})

run_gsea <- function(collection, pathway_list) {
  set.seed(20260829)
  fg <- fgseaMultilevel(pathways = pathway_list, stats = ranks, minSize = 15, maxSize = 500, eps = 0)
  fg <- as.data.frame(fg)
  fg <- fg[order(fg$padj, -abs(fg$NES), na.last = TRUE), ]
  data.frame(
    pathway = fg$pathway,
    ES = fg$ES,
    NES = fg$NES,
    pvalue = fg$pval,
    padj = fg$padj,
    geneSetSize = fg$size,
    leadingEdge = vapply(fg$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1)),
    direction = ifelse(fg$NES > 0, "Edge-enriched", "Core-enriched"),
    evidence_status = ifelse(fg$padj < 0.05, "significant_FDR_lt_0.05", ifelse(fg$padj < 0.25, "exploratory_FDR_lt_0.25", "not_significant")),
    collection = collection,
    stringsAsFactors = FALSE
  )
}

gsea <- lapply(names(pathways), function(nm) run_gsea(nm, pathways[[nm]]))
names(gsea) <- names(pathways)
gsea_files <- c(Hallmark = "GSEA_Hallmark.tsv", Reactome = "GSEA_Reactome.tsv", GO_BP = "GSEA_GO_BP.tsv", KEGG = "GSEA_KEGG.tsv")
for (nm in names(gsea)) write.table(gsea[[nm]], file.path(out_dir, gsea_files[[nm]]), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

universe <- names(ranks)
deg <- de[!is.na(de$padj) & de$padj < 0.05 & !is.na(de$gene_name) & nzchar(de$gene_name), ]
deg$gene_symbol <- toupper(deg$gene_name)
down_genes <- unique(deg$gene_symbol[deg$log2FoldChange_shrunk < 0])
up_genes <- unique(deg$gene_symbol[deg$log2FoldChange_shrunk > 0])

run_ora <- function(query, collection, pathway_list) {
  query <- intersect(unique(query), universe)
  rows <- lapply(names(pathway_list), function(term) {
    genes <- intersect(pathway_list[[term]], universe)
    if (length(genes) < 10L || length(genes) > 500L) return(NULL)
    overlap <- intersect(query, genes)
    u <- length(universe); m <- length(genes); n <- length(query); k <- length(overlap)
    p <- phyper(k - 1L, m, u - m, n, lower.tail = FALSE)
    expected <- n * m / u
    data.frame(pathway = term, pvalue = p, geneSetSize = m, querySize = n,
               overlapSize = k, foldEnrichment = ifelse(expected > 0, k / expected, NA_real_),
               genes = paste(overlap, collapse = ";"), collection = collection, stringsAsFactors = FALSE)
  })
  ans <- do.call(rbind, rows)
  ans$padj <- p.adjust(ans$pvalue, method = "BH")
  ans <- ans[order(ans$padj, -ans$foldEnrichment, na.last = TRUE), ]
  ans[, c("pathway", "pvalue", "padj", "geneSetSize", "querySize", "overlapSize", "foldEnrichment", "genes", "collection")]
}

ora_collections <- c("GO_BP", "KEGG", "Reactome")
ora_down <- lapply(ora_collections, function(nm) run_ora(down_genes, nm, pathways[[nm]])); names(ora_down) <- ora_collections
ora_up <- lapply(ora_collections, function(nm) run_ora(up_genes, nm, pathways[[nm]])); names(ora_up) <- ora_collections
ora_names <- c(GO_BP = "GO", KEGG = "KEGG", Reactome = "Reactome")
for (nm in ora_collections) {
  write.table(ora_down[[nm]], file.path(out_dir, paste0("ORA_Edge_down_", ora_names[[nm]], ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  write.table(ora_up[[nm]], file.path(out_dir, paste0("ORA_Edge_up_", ora_names[[nm]], ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

classify_program <- function(terms) {
  text <- toupper(paste(terms, collapse = " | "))
  categories <- character()
  add <- function(label, pattern) if (grepl(pattern, text, perl = TRUE)) categories <<- c(categories, label)
  add("vascular/endothelial", "ANGIO|VASCUL|ENDOTHEL|BLOOD VESSEL|VEGF|LYMPHATIC")
  add("stromal/ECM", "EXTRACELLULAR MATRIX|COLLAGEN|FIBROBLAST|FOCAL ADHESION|MATRIX ORGAN|ECM|CONNECTIVE TISSUE")
  add("immune/inflammatory", "IMMUNE|INTERFERON|CYTOKINE|CHEMOKINE|LEUKOCYTE|LYMPHOCYTE|COMPLEMENT|INFLAM|ANTIGEN|T CELL|B CELL|NATURAL KILLER")
  add("epithelial/differentiation", "EPITHEL|KERATIN|CORNIF|DIFFERENTIATION|CELL JUNCTION|ADHERENS|DESMOSOM|MUCIN")
  add("cell cycle/proliferation", "CELL CYCLE|MITOTIC|DNA REPLICATION|E2F|G2M|PROLIFER|CHROMOSOME SEGREGATION")
  add("metabolism", "METABOL|GLYCOL|OXIDATIVE PHOSPHORYLATION|FATTY ACID|AMINO ACID|MITOCHON|RESPIRATION|BIOSYNTH")
  add("stress/hypoxia", "HYPOX|STRESS|UNFOLDED PROTEIN|APOPT|REACTIVE OXYGEN")
  add("signaling/regulatory", "SIGNAL|MAPK|PI3K|WNT|TGF|NOTCH|CAMP|RAP1|KINASE|TRANSCRIPTION|SECOND MESSENGER")
  if (!length(categories)) categories <- "other"
  paste(unique(categories), collapse = ";")
}

# All-library annotation membership for every DEG.
membership <- function(gene, pathway_list) names(pathway_list)[vapply(pathway_list, function(x) gene %in% x, logical(1))]
sig_gsea <- do.call(rbind, lapply(gsea, function(x) x[!is.na(x$padj) & x$padj < 0.05, ]))
sig_leading_map <- lapply(unique(deg$gene_symbol), function(gene) {
  if (!nrow(sig_gsea)) return(character())
  sig_gsea$pathway[vapply(strsplit(sig_gsea$leadingEdge, ";", fixed = TRUE), function(x) gene %in% x, logical(1))]
}); names(sig_leading_map) <- unique(deg$gene_symbol)

collapse_terms <- function(x) if (length(x)) paste(sort(unique(x)), collapse = ";") else ""
context_rows <- lapply(seq_len(nrow(deg)), function(i) {
  gene <- deg$gene_symbol[i]
  go <- membership(gene, pathways$GO_BP)
  kegg <- membership(gene, pathways$KEGG)
  reactome <- membership(gene, pathways$Reactome)
  leading <- sig_leading_map[[gene]]
  all_terms <- unique(c(go, kegg, reactome, leading))
  data.frame(
    gene_id = deg$gene_id[i], gene_name = deg$gene_name[i], gene_type = deg$gene_type[i],
    baseMean = deg$baseMean[i], raw_LFC = deg$log2FoldChange_raw[i], shrunken_LFC = deg$log2FoldChange_shrunk[i],
    Wald_stat = deg$stat[i], pvalue = deg$pvalue[i], padj = deg$padj[i],
    direction = ifelse(deg$log2FoldChange_shrunk[i] < 0, "Edge-down", "Edge-up"),
    direction_consistency_edge_lower = deg$direction_consistency_edge_lower[i],
    direction_consistency_edge_higher = deg$direction_consistency_edge_higher[i],
    GO_annotation = collapse_terms(go), KEGG_annotation = collapse_terms(kegg), Reactome_annotation = collapse_terms(reactome),
    significant_GSEA_leading_edge = ifelse(length(leading) > 0, "yes", "no"),
    leading_edge_pathway_count = length(unique(leading)), leading_edge_pathways = collapse_terms(leading),
    major_biological_programs = classify_program(all_terms), stringsAsFactors = FALSE)
})
functional_context <- do.call(rbind, context_rows)
write.table(functional_context, file.path(out_dir, "DEG_FUNCTIONAL_CONTEXT.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

functional_context$patient_direction_consistency <- ifelse(functional_context$direction == "Edge-down", functional_context$direction_consistency_edge_lower, functional_context$direction_consistency_edge_higher)
functional_context$consistency_group <- ifelse(functional_context$patient_direction_consistency == 5, "5/5", ifelse(functional_context$patient_direction_consistency == 4, "4/5", "<=3/5"))
consistency_summary <- aggregate(gene_name ~ direction + consistency_group, functional_context, length)
names(consistency_summary)[3] <- "gene_count"
write.table(consistency_summary, file.path(out_dir, "DEG_DIRECTION_CONSISTENCY_SUMMARY.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

candidate <- functional_context[, c("gene_name", "direction", "padj", "shrunken_LFC", "Wald_stat", "baseMean", "patient_direction_consistency", "leading_edge_pathway_count", "major_biological_programs")]
names(candidate)[1] <- "gene"
names(candidate)[9] <- "major_functional_category"
candidate$rank_padj <- rank(candidate$padj, ties.method = "min")
candidate$rank_abs_Wald_stat <- rank(-abs(candidate$Wald_stat), ties.method = "min")
candidate$rank_abs_shrunken_LFC <- rank(-abs(candidate$shrunken_LFC), ties.method = "min")
candidate <- candidate[order(candidate$rank_padj, candidate$rank_abs_Wald_stat), ]
write.table(candidate, file.path(out_dir, "CANDIDATE_CONTEXT_TABLE.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# Program summary: one transparent row per significant/exploratory GSEA result or significant ORA result.
program_rows <- list()
for (nm in names(gsea)) {
  x <- gsea[[nm]][gsea[[nm]]$padj < 0.25, ]
  if (nrow(x)) program_rows[[paste0("GSEA_", nm)]] <- data.frame(
    program = vapply(x$pathway, classify_program, FUN.VALUE = character(1)), direction = x$direction,
    source = paste0("GSEA_", nm), pathway = x$pathway, NES_or_enrichment_score = x$NES, padj = x$padj,
    representative_genes = x$leadingEdge, leading_edge_genes = x$leadingEdge,
    evidence_status = x$evidence_status,
    interpretation = paste0(x$direction, " transcriptomic program; pathway-level association, not causality."), stringsAsFactors = FALSE)
}
for (direction_name in c("Edge_down", "Edge_up")) {
  collection_list <- if (direction_name == "Edge_down") ora_down else ora_up
  for (nm in names(collection_list)) {
    x <- collection_list[[nm]][collection_list[[nm]]$padj < 0.05, ]
    if (nrow(x)) program_rows[[paste0("ORA_", direction_name, "_", nm)]] <- data.frame(
      program = vapply(x$pathway, classify_program, FUN.VALUE = character(1)),
      direction = ifelse(direction_name == "Edge_down", "Edge-down DEG over-representation", "Edge-up DEG over-representation"),
      source = paste0("ORA_", nm), pathway = x$pathway, NES_or_enrichment_score = x$foldEnrichment, padj = x$padj,
      representative_genes = x$genes, leading_edge_genes = "", evidence_status = "significant_FDR_lt_0.05",
      interpretation = paste0(direction_name, " DEG over-representation; pathway-level association, not causality."), stringsAsFactors = FALSE)
  }
}
program_summary <- if (length(program_rows)) do.call(rbind, program_rows) else data.frame()
write.table(program_summary, file.path(out_dir, "EDGE_CORE_BIOLOGICAL_PROGRAM_SUMMARY.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

top_down <- de[!is.na(de$padj) & de$padj < 0.05 & de$log2FoldChange_shrunk < 0, ]
top_down <- head(top_down[order(top_down$stat), ], 10)
top_class <- functional_context[match(top_down$gene_id, functional_context$gene_id), ]
top_class <- data.frame(
  gene = top_class$gene_name, raw_LFC = top_class$raw_LFC, shrunken_LFC = top_class$shrunken_LFC,
  Wald_stat = top_class$Wald_stat, padj = top_class$padj,
  functional_categories = top_class$major_biological_programs,
  representative_annotation_terms = vapply(seq_len(nrow(top_class)), function(i) {
    terms <- unique(c(strsplit(top_class$GO_annotation[i], ";", fixed = TRUE)[[1]], strsplit(top_class$Reactome_annotation[i], ";", fixed = TRUE)[[1]], strsplit(top_class$KEGG_annotation[i], ";", fixed = TRUE)[[1]]))
    terms <- terms[nzchar(terms)]
    paste(head(terms, 12), collapse = ";")
  }, FUN.VALUE = character(1)),
  significant_leading_edge = top_class$significant_GSEA_leading_edge,
  leading_edge_pathways = top_class$leading_edge_pathways,
  stringsAsFactors = FALSE
)
write.table(top_class, file.path(out_dir, "TOP_EDGE_DOWN_FUNCTIONAL_CLASSIFICATION.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

pde <- functional_context[functional_context$gene_name == "PDE2A", , drop = FALSE]
if (nrow(pde) != 1L) stop("PDE2A does not map uniquely in the 84-DEG context")
pde$patient_direction_consistency <- pde$direction_consistency_edge_lower
pde$significant_pathway_context_status <- ifelse(pde$significant_GSEA_leading_edge == "yes", "IDENTIFIED_IN_SIGNIFICANT_GSEA_LEADING_EDGE", "NOT_IDENTIFIED_IN_SIGNIFICANT_GSEA_LEADING_EDGE")
pde$significant_ORA_GO_terms <- collapse_terms(ora_down$GO_BP$pathway[ora_down$GO_BP$padj < 0.05 & vapply(strsplit(ora_down$GO_BP$genes, ";", fixed = TRUE), function(x) "PDE2A" %in% x, logical(1))])
pde$significant_ORA_KEGG_terms <- collapse_terms(ora_down$KEGG$pathway[ora_down$KEGG$padj < 0.05 & vapply(strsplit(ora_down$KEGG$genes, ";", fixed = TRUE), function(x) "PDE2A" %in% x, logical(1))])
pde$significant_ORA_Reactome_terms <- collapse_terms(ora_down$Reactome$pathway[ora_down$Reactome$padj < 0.05 & vapply(strsplit(ora_down$Reactome$genes, ";", fixed = TRUE), function(x) "PDE2A" %in% x, logical(1))])
write.table(pde, file.path(out_dir, "PDE2A_FUNCTIONAL_CONTEXT.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# Development figures: PNG 300 dpi + source data + plotting script only.
make_fig_dir <- function(name, source_data) {
  d <- file.path(fig_root, name); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  write.table(source_data, file.path(d, "source_data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  file.copy(plot_script, file.path(d, "plotting_script.R"), overwrite = TRUE)
  status <- system2("Rscript", c(shQuote(file.path(d, "plotting_script.R")), shQuote(d)))
  if (status != 0L) stop("Plotting failed: ", name)
}
prepare_gsea_plot <- function(x, n = 20L) {
  original <- x
  x <- x[x$padj < 0.25, ]
  if (!nrow(x)) x <- head(original[order(original$padj), ], n)
  x <- head(x[order(x$padj, -abs(x$NES)), ], n)
  x$pathway_label <- make.unique(x$pathway)
  x$minus_log10_padj <- -log10(pmax(x$padj, .Machine$double.xmin))
  x
}
hallmark_plot <- prepare_gsea_plot(gsea$Hallmark, 20)
go_react <- rbind(gsea$GO_BP, gsea$Reactome)
prepare_gsea_bidirectional <- function(x, n_each = 12L) {
  x <- x[x$padj < 0.25, ]
  selected <- do.call(rbind, lapply(c("Core-enriched", "Edge-enriched"), function(direction_name) {
    y <- x[x$direction == direction_name, ]
    head(y[order(y$padj, -abs(y$NES)), ], n_each)
  }))
  selected$pathway_label <- make.unique(selected$pathway)
  selected$minus_log10_padj <- -log10(pmax(selected$padj, .Machine$double.xmin))
  selected
}
go_react <- prepare_gsea_bidirectional(go_react, 12)
make_fig_dir("Figure_Enrichment_01_GSEA_Hallmark", hallmark_plot)
make_fig_dir("Figure_Enrichment_02_GSEA_GO_Reactome", go_react)

prepare_ora_plot <- function(xlist, n_each = 10L) {
  out <- do.call(rbind, lapply(names(xlist), function(nm) {
    x <- xlist[[nm]]
    sig <- x[x$padj < 0.05 & x$overlapSize > 0, ]
    if (!nrow(sig)) sig <- x[x$overlapSize > 0, ]
    head(sig[order(sig$padj, -sig$foldEnrichment), ], n_each)
  }))
  out$pathway_label <- make.unique(out$pathway)
  out$minus_log10_padj <- -log10(pmax(out$padj, .Machine$double.xmin))
  out
}
make_fig_dir("Figure_Enrichment_03_ORA_Edge_down", prepare_ora_plot(ora_down))
make_fig_dir("Figure_Enrichment_04_ORA_Edge_up", prepare_ora_plot(ora_up))

candidate_plot <- candidate
candidate_plot$minus_log10_padj <- -log10(candidate_plot$padj)
candidate_plot$log10_baseMean_plus1 <- log10(candidate_plot$baseMean + 1)
make_fig_dir("Figure_Candidate_Context_global", candidate_plot)
make_fig_dir("Figure_Candidate_Context_PDE2A", candidate_plot)

sig_counts <- vapply(gsea, function(x) sum(!is.na(x$padj) & x$padj < 0.05), integer(1))
ora_down_sig <- vapply(ora_down, function(x) sum(!is.na(x$padj) & x$padj < 0.05), integer(1))
ora_up_sig <- vapply(ora_up, function(x) sum(!is.na(x$padj) & x$padj < 0.05), integer(1))
category_tokens <- unlist(strsplit(functional_context$major_biological_programs, ";", fixed = TRUE))
category_counts <- sort(table(category_tokens), decreasing = TRUE)
top_programs <- function(x, direction_sign) {
  y <- x[x$NES * direction_sign > 0 & x$padj < 0.05, ]
  paste(head(y$pathway[order(y$padj, -abs(y$NES))], 8), collapse = "; ")
}
pde_link <- if (pde$significant_GSEA_leading_edge == "yes" || nzchar(pde$significant_ORA_GO_terms) || nzchar(pde$significant_ORA_KEGG_terms) || nzchar(pde$significant_ORA_Reactome_terms)) {
  "PDE2A participates in at least one FDR-significant enriched context; this is associative transcriptomic context, not a causal mechanism."
} else {
  "PDE2A is an FDR-significant DEG with a large negative effect and 5/5 directional consistency, but no direct FDR-significant pathway/leading-edge context was identified."
}
audit <- c(
  "# Checkpoint 6 functional interpretation audit", "", paste0("Date: ", Sys.Date()), "",
  paste0("1. GSEA input: ", length(ranks), " unique gene symbols with valid DESeq2 Wald statistic (from 20,413 tested genes)."),
  paste0("2. Hallmark FDR < 0.05 pathways: ", sig_counts[["Hallmark"]], "."),
  paste0("3. Reactome FDR < 0.05 pathways: ", sig_counts[["Reactome"]], "."),
  paste0("4. GO Biological Process FDR < 0.05 pathways: ", sig_counts[["GO_BP"]], "."),
  paste0("5. KEGG FDR < 0.05 pathways: ", sig_counts[["KEGG"]], "."),
  paste0("6. Edge-down ORA FDR < 0.05 terms (GO/KEGG/Reactome): ", paste(ora_down_sig, collapse = "/"), "."),
  paste0("7. Edge-up ORA FDR < 0.05 terms (GO/KEGG/Reactome): ", paste(ora_up_sig, collapse = "/"), "."),
  "7a. Edge-up ORA is based on 18 genes; its significant GO/Reactome terms are largely driven by MMP1, MMP13, and MMP11 and therefore require cautious interpretation.",
  paste0("8. Functional-category annotations among 84 DEGs (non-exclusive): ", paste(names(category_counts), as.integer(category_counts), sep = "=", collapse = "; "), "."),
  paste0("9. Patient direction consistency by direction: ", paste(apply(consistency_summary, 1, paste, collapse = ":"), collapse = "; "), "."),
  paste0("10. PDE2A pathway status: ", pde$significant_pathway_context_status, "; leading-edge pathways=", ifelse(nzchar(pde$leading_edge_pathways), pde$leading_edge_pathways, "none"), "."),
  paste0("11. Top Core-enriched/Edge-down GSEA programs: ", top_programs(do.call(rbind, gsea), -1), "."),
  paste0("12. Top Edge-enriched/Edge-up GSEA programs: ", top_programs(do.call(rbind, gsea), 1), "."),
  "12a. Top Edge-down genes span vascular/endothelial and angiogenesis, epithelial/junction, immune/chemokine, ECM, and regulatory categories. This is compatible with regional cellular/tissue composition plus state differences, but does not by itself establish contamination.",
  paste0("13. ", pde_link),
  "14. Patient consistency was retained as an annotation only and was not used as a filter.",
  "15. No composite score, PDE-family analysis, WGCNA, PPI, transcription-factor network, or DEG-threshold modification was performed.",
  "16. All development figures are PNG 300 dpi with source data and plotting script; PDF/TIFF were not generated."
)
writeLines(audit, file.path(out_dir, "QC_CHECKPOINT6_FUNCTIONAL_AUDIT.md"))
writeLines(c(
  paste0("fgsea=", packageVersion("fgsea")),
  "GSEA algorithm=fgseaMultilevel; seed=20260829; minSize=15; maxSize=500; eps=0",
  "ORA=one-sided hypergeometric test; tested-symbol universe; minSize=10; maxSize=500; BH within each collection and direction",
  "Libraries=MSigDB_Hallmark_2020; Reactome_2022; GO_Biological_Process_2023; KEGG_2021_Human"
), file.path(out_dir, "ENRICHMENT_METHODS_AND_VERSIONS.txt"))
