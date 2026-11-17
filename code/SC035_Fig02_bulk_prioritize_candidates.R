#!/usr/bin/env Rscript

# Checkpoint 11: leave-one-patient-out robustness, NES<0 leading-edge core,
# and non-weighted Pareto prioritization for the frozen 35-candidate universe.
suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 25_checkpoint11_candidate_prioritization.R <project_root>")
root <- normalizePath(args[[1]], mustWork = TRUE)
out <- file.path(root, "results", "11_candidate_prioritization")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path, ...) {
  if (!file.exists(path) || file.info(path)$size == 0) stop("Missing or empty input: ", path)
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, ...)
}
write_tsv <- function(x, name) {
  write.table(x, file.path(out, name), sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA")
}
truth <- function(x) tolower(as.character(x)) %in% c("true", "t", "1", "yes")
entropy <- function(x) {
  x <- x[x > 0]
  if (length(x) <= 1L) return(0)
  p <- x / sum(x)
  -sum(p * log(p))
}

count_file <- file.path(root, "results", "04_featurecounts", "count_matrix_gene_id.tsv")
annotation_file <- file.path(root, "results", "04_featurecounts", "gene_annotation.tsv")
manifest_file <- file.path(root, "manifest_combined_standardized.tsv")
de_file <- file.path(root, "results", "05_deseq2", "DESeq2_Edge_vs_Core_all_genes.tsv")
candidate_file <- file.path(root, "results", "08_network_prioritization", "EDGE_DOWN_REGULATORY_CANDIDATES.tsv")
module_file <- file.path(root, "results", "07_candidate_convergence", "PATHWAY_TO_FUNCTIONAL_MODULE.tsv")
gsea_files <- c(
  Hallmark = file.path(root, "results", "06_enrichment", "GSEA_Hallmark.tsv"),
  Reactome = file.path(root, "results", "06_enrichment", "GSEA_Reactome.tsv"),
  GO_BP = file.path(root, "results", "06_enrichment", "GSEA_GO_BP.tsv"),
  KEGG = file.path(root, "results", "06_enrichment", "GSEA_KEGG.tsv")
)

counts <- read_tsv(count_file)
annotation <- read_tsv(annotation_file)
manifest <- read_tsv(manifest_file)
de <- read_tsv(de_file)
candidates <- read_tsv(candidate_file)
modules <- read_tsv(module_file)
if (nrow(candidates) != 35L || anyDuplicated(candidates$gene)) {
  stop("Frozen candidate universe gate failed: expected 35 unique genes")
}
if (nrow(de) != 20413L || sum(!is.na(de$padj) & de$padj < 0.05) != 84L ||
    sum(!is.na(de$padj) & de$padj < 0.05 & de$log2FoldChange_shrunk < 0) != 66L) {
  stop("Frozen Checkpoint 5 universe gate failed")
}

main_samples <- c("P08_C", "P08_TI", "P09_C", "P09_TI", "P34_C", "P34_TI",
                  "P36_C", "P36_TI", "P38_C", "P38_TI")
patients <- c("P08", "P09", "P34", "P36", "P38")
if (!all(main_samples %in% names(counts))) stop("Frozen count matrix lacks main samples")
meta <- manifest[match(main_samples, manifest$sample_id), c("sample_id", "patient_id", "region")]
if (anyNA(meta$sample_id) || any(!meta$region %in% c("c", "ti"))) stop("Main metadata gate failed")
meta$patient_id <- factor(meta$patient_id, levels = patients)
meta$region <- factor(meta$region, levels = c("c", "ti"))
rownames(meta) <- meta$sample_id

count_mat <- as.matrix(counts[, main_samples, drop = FALSE])
storage.mode(count_mat) <- "integer"
rownames(count_mat) <- counts$gene_id
if (any(count_mat < 0L) || any(count_mat != round(count_mat))) stop("Counts are not non-negative integers")

gene_lookup <- de[, c("gene_id", "gene_name")]
candidate_ids <- gene_lookup$gene_id[match(candidates$gene, gene_lookup$gene_name)]
if (anyNA(candidate_ids)) stop("Candidate gene-to-ID mapping failed: ", paste(candidates$gene[is.na(candidate_ids)], collapse = ", "))
if (anyDuplicated(candidate_ids)) stop("Candidate gene IDs are not unique")

# Five independent four-patient sensitivity fits. These do not redefine the candidate set.
lopo_lfc <- matrix(NA_real_, nrow = nrow(candidates), ncol = length(patients),
                   dimnames = list(candidates$gene, paste0("LOPO_", patients, "_LFC")))
lopo_filter <- data.frame(analysis = character(), omitted_patient = character(),
                          samples = integer(), genes_before = integer(), genes_retained = integer(),
                          genes_removed = integer(), coefficient = character(), stringsAsFactors = FALSE)

for (i in seq_along(patients)) {
  omitted <- patients[[i]]
  keep_sample <- meta$patient_id != omitted
  sub_meta <- droplevels(meta[keep_sample, , drop = FALSE])
  sub_counts <- count_mat[, rownames(sub_meta), drop = FALSE]
  keep_gene <- rowSums(sub_counts >= 10L) >= 3L
  filtered <- sub_counts[keep_gene, , drop = FALSE]
  dds <- DESeqDataSetFromMatrix(countData = filtered, colData = sub_meta,
                                design = ~ patient_id + region)
  dds <- DESeq(dds, quiet = TRUE)
  coef_name <- grep("^region_.*ti.*_vs_.*c$", resultsNames(dds), value = TRUE)
  if (length(coef_name) != 1L) stop("Could not identify region coefficient for ", omitted,
                                    ": ", paste(resultsNames(dds), collapse = ", "))
  shr <- lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)
  vals <- shr$log2FoldChange[match(candidate_ids, rownames(shr))]
  lopo_lfc[, i] <- vals
  lopo_filter <- rbind(lopo_filter, data.frame(
    analysis = paste0("LOPO_", omitted), omitted_patient = omitted,
    samples = ncol(filtered), genes_before = nrow(sub_counts), genes_retained = nrow(filtered),
    genes_removed = nrow(sub_counts) - nrow(filtered), coefficient = coef_name,
    stringsAsFactors = FALSE
  ))
}
write_tsv(lopo_filter, "LOPO_LOW_COUNT_FILTER_AUDIT.tsv")

negative_count <- rowSums(lopo_lfc < 0, na.rm = TRUE)
available_count <- rowSums(!is.na(lopo_lfc))
rank_mat <- apply(lopo_lfc, 2, function(x) rank(x, ties.method = "min", na.last = "keep"))
rank_iqr <- apply(rank_mat, 1, IQR, na.rm = TRUE)
lopo <- data.frame(
  gene = candidates$gene,
  full_5pair_shrunken_LFC = candidates$shrunken_log2FC,
  lopo_lfc,
  LOPO_available_count = available_count,
  LOPO_negative_count = negative_count,
  LOPO_median_LFC = apply(lopo_lfc, 1, median, na.rm = TRUE),
  LOPO_min_LFC = apply(lopo_lfc, 1, min, na.rm = TRUE),
  LOPO_max_LFC = apply(lopo_lfc, 1, max, na.rm = TRUE),
  LOPO_sign_stability = ifelse(available_count < 5L, "INCOMPLETE",
                               ifelse(negative_count == 5L, "STABLE_NEGATIVE",
                                      ifelse(negative_count == 4L, "MOSTLY_NEGATIVE", "UNSTABLE"))),
  median_LOPO_effect_rank = apply(rank_mat, 1, median, na.rm = TRUE),
  best_LOPO_effect_rank = apply(rank_mat, 1, min, na.rm = TRUE),
  worst_LOPO_effect_rank = apply(rank_mat, 1, max, na.rm = TRUE),
  rank_IQR = rank_iqr,
  rank_universe = "35 frozen Edge-down signaling/regulatory candidates",
  stringsAsFactors = FALSE, check.names = FALSE
)
write_tsv(lopo, "LOPO_CANDIDATE_STABILITY.tsv")
# A candidate can legitimately fail the frozen low-count rule in one four-patient
# subset. It remains in the 35-gene universe with NA and INCOMPLETE status; the
# filter is not relaxed and the candidate is not deleted.

# Build the all-direction and primary NES<0 gene-pathway edge tables from frozen GSEA.
all_edges <- list()
for (db in names(gsea_files)) {
  g <- read_tsv(gsea_files[[db]])
  g <- g[!is.na(g$padj) & g$padj < 0.05, , drop = FALSE]
  for (j in seq_len(nrow(g))) {
    genes <- unique(trimws(strsplit(g$leadingEdge[[j]], ";", fixed = TRUE)[[1]]))
    genes <- intersect(genes, candidates$gene)
    if (!length(genes)) next
    all_edges[[length(all_edges) + 1L]] <- data.frame(
      database = db, pathway = g$pathway[[j]], NES = g$NES[[j]], padj = g$padj[[j]],
      direction = ifelse(g$NES[[j]] < 0, "Core-enriched / Edge-down", "Edge-enriched / Edge-up"),
      gene = genes, stringsAsFactors = FALSE
    )
  }
}
all_edges <- if (length(all_edges)) do.call(rbind, all_edges) else data.frame()
if (!nrow(all_edges)) stop("No candidate-leading-edge edges were recovered")
module_key <- paste(modules$database, modules$pathway, sep = "\r")
all_edges$functional_module <- modules$functional_module[match(paste(all_edges$database, all_edges$pathway, sep = "\r"), module_key)]
all_edges$functional_module[is.na(all_edges$functional_module) | !nzchar(all_edges$functional_module)] <- "Other / unresolved"
all_edges <- all_edges[order(all_edges$database, all_edges$pathway, all_edges$gene), ]
write_tsv(all_edges, "LEADING_EDGE_BIPARTITE_ALL_DIRECTIONS.tsv")
primary_edges <- all_edges[all_edges$NES < 0, , drop = FALSE]
write_tsv(primary_edges, "LEADING_EDGE_BIPARTITE_NES_NEGATIVE.tsv")

# Candidate-level primary pathway-core metrics.
core_rows <- lapply(candidates$gene, function(gene) {
  x <- primary_edges[primary_edges$gene == gene, , drop = FALSE]
  module_counts <- table(x$functional_module)
  data.frame(
    gene = gene,
    edge_down_leading_edge_count = nrow(x),
    functional_module_count = length(unique(x$functional_module)),
    database_support_count = length(unique(x$database)),
    module_participation_count = length(unique(x$functional_module)),
    module_entropy = entropy(as.numeric(module_counts)),
    bipartite_gene_degree = nrow(x),
    primary_functional_modules = paste(sort(unique(x$functional_module)), collapse = ";"),
    primary_database_support = paste(sort(unique(x$database)), collapse = ";"),
    stringsAsFactors = FALSE
  )
})
core <- do.call(rbind, core_rows)
write_tsv(core, "PATHWAY_CORE_CANDIDATE_METRICS.tsv")

driver <- merge(candidates, lopo, by.x = "gene", by.y = "gene", all.x = TRUE, sort = FALSE)
driver <- merge(driver, core, by = "gene", all.x = TRUE, sort = FALSE)
driver <- driver[match(candidates$gene, driver$gene), ]
driver_out <- data.frame(
  gene = driver$gene,
  full_shrunken_LFC = driver$shrunken_log2FC,
  full_padj = driver$padj,
  patient_consistency = driver$patient_consistency,
  LOPO_negative_count = driver$LOPO_negative_count,
  LOPO_median_LFC = driver$LOPO_median_LFC,
  edge_down_leading_edge_count = driver$edge_down_leading_edge_count,
  functional_module_count = driver$functional_module_count.y,
  database_support_count = driver$database_support_count,
  module_entropy = driver$module_entropy,
  functional_class = driver$functional_class,
  enzyme_yes_no = driver$enzyme_yes_no,
  receptor_yes_no = driver$receptor_yes_no,
  structural_marker_yes_no = driver$structural_marker_yes_no,
  PPI_degree = driver$PPI_degree,
  PPI_community = driver$PPI_community,
  primary_functional_modules = driver$primary_functional_modules,
  stringsAsFactors = FALSE
)
write_tsv(driver_out, "PATHWAY_CORE_DRIVER_CANDIDATES.tsv")

# Gene-by-module pathway-count matrix for all 35 candidates.
module_levels <- sort(unique(primary_edges$functional_module))
module_matrix <- matrix(0L, nrow = nrow(candidates), ncol = length(module_levels),
                        dimnames = list(candidates$gene, module_levels))
for (i in seq_len(nrow(primary_edges))) {
  module_matrix[primary_edges$gene[[i]], primary_edges$functional_module[[i]]] <-
    module_matrix[primary_edges$gene[[i]], primary_edges$functional_module[[i]]] + 1L
}
module_matrix_out <- data.frame(gene = rownames(module_matrix), module_matrix,
                                check.names = FALSE, row.names = NULL)
write_tsv(module_matrix_out, "GENE_FUNCTIONAL_MODULE_NES_NEGATIVE_MATRIX.tsv")

# Pareto non-dominated sorting: all criteria are transformed so larger is better.
pareto_fronts <- function(df, criteria) {
  values <- as.matrix(df[, criteria, drop = FALSE])
  if (anyNA(values)) stop("Pareto criteria contain NA")
  remaining <- seq_len(nrow(values))
  fronts <- integer(nrow(values))
  front <- 1L
  while (length(remaining)) {
    non_dominated <- vapply(remaining, function(i) {
      !any(vapply(setdiff(remaining, i), function(j) {
        all(values[j, ] >= values[i, ]) && any(values[j, ] > values[i, ])
      }, logical(1)))
    }, logical(1))
    selected <- remaining[non_dominated]
    fronts[selected] <- front
    remaining <- setdiff(remaining, selected)
    front <- front + 1L
  }
  fronts
}

pareto <- driver_out
pareto$criterion_edge_down_effect <- -pareto$full_shrunken_LFC
pareto$criterion_LOPO_robustness <- pareto$LOPO_negative_count
pareto$criterion_leading_edge_recurrence <- pareto$edge_down_leading_edge_count
pareto$criterion_module_breadth <- pareto$functional_module_count
pareto$criterion_database_breadth <- pareto$database_support_count
criteria <- c("criterion_edge_down_effect", "criterion_LOPO_robustness",
              "criterion_leading_edge_recurrence", "criterion_module_breadth",
              "criterion_database_breadth")
pareto$pareto_front <- pareto_fronts(pareto, criteria)
pareto$rank_universe <- "35 frozen Edge-down signaling/regulatory candidates"
pareto <- pareto[order(pareto$pareto_front, pareto$full_padj), ]
write_tsv(pareto, "CANDIDATE_PARETO_FRONTS.tsv")

sensitivity <- data.frame(gene = pareto$gene, stringsAsFactors = FALSE)
for (criterion in criteria) {
  label <- paste0("front1_without_", sub("^criterion_", "", criterion))
  sensitivity[[label]] <- pareto_fronts(pareto, setdiff(criteria, criterion)) == 1L
}
sensitivity$pareto_front1_frequency <- rowSums(sensitivity[, -1, drop = FALSE])
sensitivity$sensitivity_analysis_count <- 5L
sensitivity$rank_universe <- "35 frozen Edge-down signaling/regulatory candidates"
write_tsv(sensitivity, "PARETO_SENSITIVITY.tsv")

natural <- c("ADCY4", "KDR", "CAMK2B", "GNG11", "PDE1A", "PDE2A", "NR4A1")
comparison_genes <- unique(c(intersect(natural, pareto$gene), pareto$gene[pareto$pareto_front == 1L]))
comparison <- merge(pareto[pareto$gene %in% comparison_genes, ], sensitivity, by = "gene", sort = FALSE)
comparison <- comparison[match(comparison_genes, comparison$gene), ]
write_tsv(comparison, "PARETO_TOP_CANDIDATE_COMPARISON.tsv")

# Transparent functional-role annotation, derived from frozen annotations only.
role <- function(cls, enzyme, receptor, gene) {
  cls_l <- tolower(ifelse(is.na(cls), "", cls))
  roles <- character()
  if (truth(enzyme)) roles <- c(roles, "intracellular enzyme")
  if (grepl("kinase", cls_l)) roles <- c(roles, "kinase")
  if (truth(receptor)) roles <- c(roles, "receptor")
  if (gene %in% c("GNG11")) roles <- c(roles, "G-protein/signaling component")
  if (grepl("transcription", cls_l) || gene %in% c("NR4A1", "KLF6", "SOX6", "SOX18")) roles <- c(roles, "transcription regulator")
  if (grepl("vascular|endothelial", cls_l)) roles <- c(roles, "vascular/endothelial-associated regulator")
  if (truth(driver_out$structural_marker_yes_no[match(gene, driver_out$gene)])) roles <- c(roles, "structural/marker")
  if (!length(roles)) roles <- "other"
  paste(unique(roles), collapse = ";")
}
role_table <- pareto[pareto$pareto_front %in% c(1L, 2L),
                     c("gene", "pareto_front", "functional_class", "enzyme_yes_no", "receptor_yes_no")]
role_table$functional_role <- mapply(role, role_table$functional_class,
                                     role_table$enzyme_yes_no, role_table$receptor_yes_no,
                                     role_table$gene, USE.NAMES = FALSE)
write_tsv(role_table, "PARETO_FRONT1_2_FUNCTIONAL_ROLES.tsv")

# PDE2A is extracted only after the complete 35-candidate tables are frozen.
pde <- merge(pareto[pareto$gene == "PDE2A", ], sensitivity[sensitivity$gene == "PDE2A", ], by = "gene")
pde_lopo <- lopo[lopo$gene == "PDE2A", ]
pde <- merge(pde, pde_lopo[, c("gene", colnames(lopo_lfc), "LOPO_median_LFC",
                               "best_LOPO_effect_rank", "worst_LOPO_effect_rank", "rank_IQR")], by = "gene")
write_tsv(pde, "PDE2A_CHECKPOINT11_CONTEXT.tsv")

# Development figure source data; plotting is delegated to one frozen script.
fig_root <- file.path(root, "results", "07_figures", "manuscript_candidates", "Candidate_nomination")
fig_specs <- c("Figure_Nomination_01_Funnel", "Figure_Nomination_02_LOPO_Stability",
               "Figure_Nomination_03_LeadingEdge_CoreMap", "Figure_Nomination_04_Pareto_Global",
               "Figure_Nomination_04_PDE2A_Context")
for (d in fig_specs) dir.create(file.path(fig_root, d), recursive = TRUE, showWarnings = FALSE)
funnel <- data.frame(stage = c("Tested genes", "FDR-significant DEGs", "Edge-down DEGs",
                               "Edge-down signaling/regulatory candidates", "Pareto Front 1 candidates", "PDE2A prioritized"),
                     n = c(20413L, 84L, 66L, 35L, sum(pareto$pareto_front == 1L), 1L),
                     stage_order = 1:6)
write.table(funnel, file.path(fig_root, fig_specs[[1]], "source_data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lopo, file.path(fig_root, fig_specs[[2]], "source_data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
display_genes <- unique(c(pareto$gene[pareto$pareto_front %in% c(1L, 2L)], intersect(natural, pareto$gene)))
heat <- module_matrix_out[module_matrix_out$gene %in% display_genes, , drop = FALSE]
heat <- merge(heat, pareto[, c("gene", "pareto_front", "full_shrunken_LFC", "LOPO_negative_count", "patient_consistency")], by = "gene")
write.table(heat, file.path(fig_root, fig_specs[[3]], "source_data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
plot_data <- pareto[, c("gene", "full_shrunken_LFC", "edge_down_leading_edge_count", "functional_module_count", "pareto_front", "full_padj")]
write.table(plot_data, file.path(fig_root, fig_specs[[4]], "source_data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(plot_data, file.path(fig_root, fig_specs[[5]], "source_data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

plot_script <- file.path(root, "scripts", "26_plot_checkpoint11_nomination.R")
for (d in fig_specs) file.copy(plot_script, file.path(fig_root, d, "plotting_script.R"), overwrite = TRUE)
for (d in fig_specs) {
  status <- system2("Rscript", c(shQuote(file.path(fig_root, d, "plotting_script.R")), shQuote(file.path(fig_root, d))))
  if (status != 0L) stop("Plotting failed: ", d)
}

stability_counts <- table(factor(lopo$LOPO_sign_stability,
                                 levels = c("STABLE_NEGATIVE", "MOSTLY_NEGATIVE", "UNSTABLE", "INCOMPLETE")))
front1 <- pareto$gene[pareto$pareto_front == 1L]
front2 <- pareto$gene[pareto$pareto_front == 2L]
pde_p <- pareto[pareto$gene == "PDE2A", ]
pde_s <- sensitivity[sensitivity$gene == "PDE2A", ]
pde_l <- lopo[lopo$gene == "PDE2A", ]
nomination_supported <- pde_p$pareto_front <= 2L && pde_l$LOPO_negative_count >= 4L &&
  pde_p$edge_down_leading_edge_count > 0L && truth(pde_p$enzyme_yes_no)

audit <- c(
  "# Checkpoint 11 final nomination audit", "",
  "## Frozen scope", "",
  "- Candidate universe: 35 Edge-down DEGs previously annotated as signaling/regulatory functional-class candidates in Checkpoint 8.",
  "- Reference universe: 20,413 tested genes -> 84 FDR-significant DEGs -> 66 Edge-down DEGs -> 35 annotated Edge-down signaling/regulatory candidates.",
  "- No candidate was added, removed, or reclassified in Checkpoint 11.", "",
  "## Leave-one-patient-out robustness", "",
  "- Five LOPO paired DESeq2 analyses completed: P08, P09, P34, P36, and P38 omitted in turn.",
  "- Each model used raw integer counts, `design = ~ patient_id + region`, and the frozen Edge/TI - Core contrast.",
  "- Each subset independently applied `rowSums(counts >= 10) >= 3`; retained-gene counts are in `LOPO_LOW_COUNT_FILTER_AUDIT.tsv`.",
  paste0("- 5/5 LOPO negative: ", stability_counts[["STABLE_NEGATIVE"]], " candidates."),
  paste0("- 4/5 LOPO negative: ", stability_counts[["MOSTLY_NEGATIVE"]], " candidates."),
  paste0("- <=3/5 LOPO negative: ", stability_counts[["UNSTABLE"]], " candidates."),
  paste0("- Incomplete under the frozen low-count rule: ", stability_counts[["INCOMPLETE"]],
         " candidate (retained in the 35-gene universe with NA; no filter relaxation)."),
  paste0("- PDE2A LOPO LFCs: ", paste(sprintf("%.4f", as.numeric(pde_l[, colnames(lopo_lfc)])), collapse = ", "), "."),
  paste0("- PDE2A LOPO negative count: ", pde_l$LOPO_negative_count, "/5; median LFC = ", sprintf("%.4f", pde_l$LOPO_median_LFC),
         "; candidate-universe rank range = ", pde_l$best_LOPO_effect_rank, "-", pde_l$worst_LOPO_effect_rank,
         "; rank IQR = ", sprintf("%.2f", pde_l$rank_IQR), "."), "",
  "## NES<0 leading-edge pathway core", "",
  paste0("- Primary bipartite network: ", length(unique(primary_edges$gene)), " of 35 candidate genes, ",
         length(unique(paste(primary_edges$database, primary_edges$pathway))), " FDR-significant NES<0 pathways, and ", nrow(primary_edges), " gene-pathway edges."),
  paste0("- Functional modules represented: ", length(unique(primary_edges$functional_module)), "."),
  paste0("- Highest pathway-core recurrence genes: ", paste(head(driver_out$gene[order(-driver_out$edge_down_leading_edge_count, driver_out$full_padj)], 10), collapse = ", "), "."), "",
  "## Pareto prioritization", "",
  "- Primary criteria: more negative full shrunken LFC, higher LOPO negative count, more NES<0 significant leading-edge memberships, broader functional-module participation, and broader database support.",
  "- PPI degree and patient consistency were annotations, not Pareto criteria.",
  paste0("- Pareto Front 1: ", paste(front1, collapse = ", "), "."),
  paste0("- Pareto Front 2: ", paste(front2, collapse = ", "), "."),
  paste0("- PDE2A Pareto front: ", pde_p$pareto_front, "."),
  paste0("- PDE2A Front-1 frequency across five leave-one-criterion-out analyses: ", pde_s$pareto_front1_frequency, "/5."), "",
  "## Nomination decision", "",
  if (nomination_supported) "PDE2A can be described as a robust and functionally convergent Edge-down regulatory candidate within a transparent multidimensional candidate set."
  else "PDE2A does not satisfy the prespecified evidence pattern for the strongest nomination wording; its measured evidence remains reported without changing criteria.",
  "PDE2A was not required to be the top DEG, the top PPI hub, or the unique candidate.",
  "No weighted or composite score, machine learning, post-hoc threshold change, PPI rerun, or mechanism analysis was performed.",
  "Checkpoint 9 external discordance and Checkpoint 10 null-FDR activity results remain unchanged and were not used as positive nomination evidence."
)
writeLines(audit, file.path(out, "QC_CHECKPOINT11_FINAL_NOMINATION_AUDIT.md"))

write_tsv(data.frame(field = c("run_time", "R_version", "DESeq2_version", "apeglm_version",
                               "candidate_universe_n", "tested_gene_reference_n", "FDR_DEG_reference_n", "Edge_down_reference_n"),
                     value = c(as.character(Sys.time()), R.version.string, as.character(packageVersion("DESeq2")),
                               as.character(packageVersion("apeglm")), 35L, 20413L, 84L, 66L)),
          "CHECKPOINT11_RUN_METADATA.tsv")
message("Checkpoint 11 completed.")
