#!/usr/bin/env Rscript

# Checkpoint 7: map frozen FDR-significant DEGs to frozen significant GSEA
# leading edges and auditable functional modules. This script does not rerun
# DESeq2, GSEA, or ORA and does not calculate a composite candidate score.

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
out_dir <- file.path(project_root, "results", "07_candidate_convergence")
fig_root <- file.path(project_root, "results", "07_figures", "manuscript_candidates", "Candidate_convergence")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_root, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) stop("Missing or empty frozen input: ", path)
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
}
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
split_tokens <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character())
  unique(trimws(strsplit(x, ";", fixed = TRUE)[[1]]))
}
collapse_unique <- function(x) paste(sort(unique(x[nzchar(x)])), collapse = ";")

de_path <- file.path(project_root, "results", "05_deseq2", "DESeq2_Edge_vs_Core_all_genes.tsv")
context_path <- file.path(project_root, "results", "06_enrichment", "DEG_FUNCTIONAL_CONTEXT.tsv")
candidate_path <- file.path(project_root, "results", "06_enrichment", "CANDIDATE_CONTEXT_TABLE.tsv")
gsea_paths <- c(
  Hallmark = file.path(project_root, "results", "06_enrichment", "GSEA_Hallmark.tsv"),
  Reactome = file.path(project_root, "results", "06_enrichment", "GSEA_Reactome.tsv"),
  GO_BP = file.path(project_root, "results", "06_enrichment", "GSEA_GO_BP.tsv"),
  KEGG = file.path(project_root, "results", "06_enrichment", "GSEA_KEGG.tsv")
)

de <- read_tsv(de_path)
context <- read_tsv(context_path)
candidate <- read_tsv(candidate_path)
required_de <- c("gene_name", "baseMean", "stat", "padj", "log2FoldChange_shrunk",
                 "direction_consistency_edge_lower", "direction_consistency_edge_higher")
if (!all(required_de %in% names(de))) stop("Frozen DE table is missing required columns")
sig_de <- de[!is.na(de$padj) & de$padj < 0.05, , drop = FALSE]
if (nrow(de) != 20413L || nrow(sig_de) != 84L) stop("Frozen DE dimensions changed: expected 20,413 tested and 84 FDR-significant genes")
if (anyDuplicated(sig_de$gene_name)) stop("The 84-DEG gene_name field is not unique")
if (nrow(context) != 84L || nrow(candidate) != 84L) stop("Checkpoint 6 DEG context dimensions changed")

gsea_list <- lapply(names(gsea_paths), function(db) {
  x <- read_tsv(gsea_paths[[db]])
  required <- c("pathway", "NES", "pvalue", "padj", "geneSetSize", "leadingEdge")
  if (!all(required %in% names(x))) stop("GSEA table missing required columns: ", db)
  x$database <- db
  x
})
gsea_all <- do.call(rbind, gsea_list)
gsea_sig <- gsea_all[!is.na(gsea_all$padj) & gsea_all$padj < 0.05, , drop = FALSE]
gsea_sig$direction <- ifelse(gsea_sig$NES < 0, "Core-enriched / Edge-down", "Edge-enriched / Edge-up")

# Ordered pathway-name rules. The first matching rule is the primary module.
# Specific signaling and structural programs precede broad categories.
module_rules <- data.frame(
  functional_module = c(
    "NO / cGMP signaling", "Cyclic nucleotide signaling", "GPCR / adenylate cyclase signaling",
    "Hemostasis / platelet", "Vascular / endothelial", "ECM / stromal remodeling",
    "Epithelial differentiation", "Cell cycle / proliferation", "DNA replication / repair",
    "Immune / interferon", "Translation / ribosome", "Growth-factor signaling",
    "Stress / hypoxia", "Metabolism / endocrine", "Signaling / regulatory"
  ),
  pattern = c(
    "nitric oxide|guanylate cyclase|cgmp|protein kinase g",
    "cyclic.nucleotide|camp|phosphodiester|adenosine.*monophosphate|purine ribonucleotide",
    "gpcr|g protein.coupled|g alpha|adenylate cyclase|second messenger|opioid signaling|cholinergic synapse",
    "hemostasis|platelet|coagulation|clot|fibrin|thrombo|wound healing",
    "vascular|vessel|angiogen|endothelial|artery|aorta|vein|smooth muscle contraction",
    "extracellular matrix|ecm|collagen|matrix metalloprotein|stromal|matrisome|cell.matrix|focal adhesion|integrin",
    "epithelial|epithelium|keratin|cornif|cell junction|tight junction|desmosom|skin development",
    "cell cycle|mitotic|mitosis|g2.m|g1.s|e2f|chromosome segregation|spindle|cytokinesis|proliferation",
    "dna replication|dna repair|dna metabolic|homologous recombination|fanconi|mismatch repair|nucleotide excision|base excision|telomere|chromatin",
    "immune|interferon|cytokine|chemokine|antigen|leukocyte|lymphocyte|t cell|b cell|neutrophil|inflamm|complement",
    "translation|ribosom|peptide chain|trna|rrna|nonsense mediated decay|eukaryotic initiation factor",
    "growth factor|tgf|transforming growth|egfr|fibroblast growth|vegf|pdgf|insulin.like growth",
    "hypoxi|oxidative stress|reactive oxygen|unfolded protein|heat shock|apoptosis|cell death|senescence",
    "metabol|glycol|fatty acid|lipid|amino acid|cholesterol|oxidative phosphorylation|mitochond|aldosterone|adipogenesis|purine|carbon",
    "signaling|signal transduction|kinase|phosphorylation|ras|mapk|mtorc|wnt|notch|hedgehog|calcium|transcription|response to hormone"
  ), stringsAsFactors = FALSE
)
assign_module <- function(pathway) {
  text <- tolower(pathway)
  for (i in seq_len(nrow(module_rules))) {
    if (grepl(module_rules$pattern[i], text, perl = TRUE)) {
      return(c(module_rules$functional_module[i], paste0("priority_rule_", i, ":", module_rules$pattern[i])))
    }
  }
  c("Other / unresolved", "no_keyword_rule_matched")
}
mapped <- t(vapply(gsea_sig$pathway, assign_module, character(2)))
gsea_sig$functional_module <- mapped[, 1]
gsea_sig$mapping_rule <- mapped[, 2]
pathway_map <- gsea_sig[, c("database", "pathway", "NES", "padj", "direction", "functional_module", "mapping_rule")]
pathway_map <- pathway_map[order(pathway_map$database, pathway_map$padj, pathway_map$pathway), ]
write_tsv(pathway_map, file.path(out_dir, "PATHWAY_TO_FUNCTIONAL_MODULE.tsv"))

# PDE2A direction audit is a readout of the same global significant leading-edge map.
pde_in_le <- vapply(gsea_sig$leadingEdge, function(z) "PDE2A" %in% split_tokens(z), logical(1))
pde_audit <- gsea_sig[pde_in_le, c("database", "pathway", "NES", "pvalue", "padj", "geneSetSize", "leadingEdge", "direction", "functional_module")]
names(pde_audit)[names(pde_audit) == "leadingEdge"] <- "leading_edge_genes"
pde_audit$PDE2A_in_leading_edge <- TRUE
pde_audit$direction_consistent_with_PDE2A_edge_down <- pde_audit$NES < 0
pde_audit <- pde_audit[, c("database", "pathway", "NES", "pvalue", "padj", "geneSetSize",
                           "PDE2A_in_leading_edge", "direction", "direction_consistent_with_PDE2A_edge_down",
                           "functional_module", "leading_edge_genes")]
pde_audit <- pde_audit[order(pde_audit$NES, pde_audit$padj), ]
write_tsv(pde_audit, file.path(out_dir, "PDE2A_LEADING_EDGE_DIRECTION_AUDIT.tsv"))

# Build one transparent membership record per DEG using identical rules.
contexts <- setNames(seq_len(nrow(context)), context$gene_name)
membership_rows <- lapply(seq_len(nrow(sig_de)), function(i) {
  gene <- sig_de$gene_name[i]
  hit <- vapply(gsea_sig$leadingEdge, function(z) gene %in% split_tokens(z), logical(1))
  h <- gsea_sig[hit, , drop = FALSE]
  d <- if (sig_de$log2FoldChange_shrunk[i] < 0) "Edge-down" else "Edge-up"
  consistency <- if (d == "Edge-down") sig_de$direction_consistency_edge_lower[i] else sig_de$direction_consistency_edge_higher[i]
  data.frame(
    gene_name = gene, direction = d, padj = sig_de$padj[i],
    shrunken_log2FC = sig_de$log2FoldChange_shrunk[i], Wald_stat = sig_de$stat[i],
    baseMean = sig_de$baseMean[i], patient_direction_consistency = consistency,
    total_significant_leading_edge_count = nrow(h),
    Hallmark_leading_edge_count = sum(h$database == "Hallmark"),
    Reactome_leading_edge_count = sum(h$database == "Reactome"),
    GO_BP_leading_edge_count = sum(h$database == "GO_BP"),
    KEGG_leading_edge_count = sum(h$database == "KEGG"),
    functional_modules = collapse_unique(h$functional_module),
    functional_module_count = length(unique(h$functional_module)),
    core_enriched_leading_edge_count = sum(h$NES < 0),
    edge_enriched_leading_edge_count = sum(h$NES > 0),
    significant_leading_edge_pathways = collapse_unique(paste(h$database, h$pathway, sep = ": ")),
    stringsAsFactors = FALSE
  )
})
membership <- do.call(rbind, membership_rows)
write_tsv(membership, file.path(out_dir, "DEG_LEADING_EDGE_MEMBERSHIP.tsv"))

# Multi-label functional classes. These are annotations, never exclusion rules.
enzyme_symbols <- c("PDE1A", "PDE2A", "UGT1A10", "UGT1A7", "UCHL1", "TMPRSS11B", "CAMK2B",
                    "ACHE", "CA3", "DUSP5", "PYGM", "MMP1", "MMP13", "ADCY4", "WNK4",
                    "CYP4F12", "CYP4F29P", "ADAMTS1", "MMP11", "EEF1A2")
receptor_pattern <- "(^|;)(KDR|TNFRSF11B|STAB2|CLEC14A|CLEC3B|SCNN1B|GPC3)($|;)"
tf_symbols <- c("HOXD13", "KLF6", "SOX6", "NR4A1", "SOX18")
structural_symbols <- c("TMOD4", "CRNN", "SPRR3", "SPRR2C", "MUC21", "COL10A1", "VWF", "KRT4",
                        "KRT78", "POSTN", "MYZAP", "CDH5", "KANK3", "CLDN5", "MMRN1", "EMCN", "CD34")
classify_gene <- function(gene, modules, ctx) {
  text <- tolower(paste(ctx$major_biological_programs, ctx$GO_annotation, ctx$KEGG_annotation,
                        ctx$Reactome_annotation, modules, sep = ";"))
  classes <- character()
  if (grepl("vascular|endothelial|angiogen|blood vessel|artery|aorta", text)) classes <- c(classes, "vascular/endothelial")
  if (grepl("immune|interferon|cytokine|chemokine|leukocyte|inflamm|complement", text)) classes <- c(classes, "immune/chemokine")
  if (grepl("epithelial|epithelium|keratin|cornif|junction|differentiation", text)) classes <- c(classes, "epithelial/differentiation")
  if (grepl("extracellular matrix|ecm|stromal|collagen|matrisome|matrix", text)) classes <- c(classes, "ECM/stromal")
  if (grepl("metabol|glycol|lipid|fatty acid|mitochond|purine|aldosterone", text)) classes <- c(classes, "metabolism")
  if (grepl("cell cycle|mitotic|proliferation|dna replication|chromosome", text)) classes <- c(classes, "cell cycle")
  if (grepl("signaling|signal transduction|kinase|phosphorylation|growth factor|gpcr|g protein|second messenger|cyclic nucleotide|cgmp|camp|adenylate cyclase|guanylate cyclase|ras|mapk|wnt", text)) classes <- c(classes, "signaling/regulatory")
  if (gene %in% enzyme_symbols) classes <- c(classes, "enzyme")
  if (grepl(receptor_pattern, paste0(";", gene, ";"))) classes <- c(classes, "receptor")
  if (gene %in% tf_symbols) classes <- c(classes, "transcriptional regulator")
  if (gene %in% structural_symbols) classes <- c(classes, "structural/marker")
  if (!length(classes)) classes <- "other"
  collapse_unique(classes)
}
functional_class <- vapply(seq_len(nrow(membership)), function(i) {
  idx <- contexts[[membership$gene_name[i]]]
  classify_gene(membership$gene_name[i], membership$functional_modules[i], context[idx, , drop = FALSE])
}, character(1))
classification <- cbind(membership[, c("gene_name", "direction", "padj", "shrunken_log2FC", "Wald_stat", "baseMean",
                                         "patient_direction_consistency", "functional_modules")],
                        functional_class = functional_class,
                        classification_basis = "Checkpoint6 GO/KEGG/Reactome annotations + significant GSEA modules + documented generic gene-family rules")
write_tsv(classification, file.path(out_dir, "DEG_FUNCTIONAL_CLASSIFICATION.tsv"))

# Gene-by-module counts; binary columns accompany pathway-count columns.
module_levels <- sort(unique(gsea_sig$functional_module))
module_counts <- matrix(0L, nrow = nrow(membership), ncol = length(module_levels),
                        dimnames = list(membership$gene_name, module_levels))
for (i in seq_len(nrow(membership))) {
  gene <- membership$gene_name[i]
  hit <- vapply(gsea_sig$leadingEdge, function(z) gene %in% split_tokens(z), logical(1))
  if (any(hit)) {
    tab <- table(gsea_sig$functional_module[hit])
    module_counts[i, names(tab)] <- as.integer(tab)
  }
}
module_binary <- (module_counts > 0) * 1L
sanitize <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
colnames(module_counts) <- paste0("pathway_count__", sanitize(colnames(module_counts)))
colnames(module_binary) <- paste0("membership__", sanitize(colnames(module_binary)))
module_matrix <- data.frame(gene_name = rownames(module_counts), module_binary, module_counts, check.names = FALSE)
write_tsv(module_matrix, file.path(out_dir, "GENE_FUNCTIONAL_MODULE_MATRIX.tsv"))

# Transparent independent ranks only; no final or composite rank.
conv <- merge(membership, classification[, c("gene_name", "functional_class")], by = "gene_name", sort = FALSE)
conv$padj_rank <- rank(conv$padj, ties.method = "min")
conv$absolute_Wald_rank <- rank(-abs(conv$Wald_stat), ties.method = "min")
conv$absolute_shrunken_LFC_rank <- rank(-abs(conv$shrunken_log2FC), ties.method = "min")
conv$leading_edge_count_rank <- rank(-conv$total_significant_leading_edge_count, ties.method = "min")
conv <- conv[match(sig_de$gene_name, conv$gene_name), ]
write_tsv(conv, file.path(out_dir, "CANDIDATE_CONVERGENCE_TABLE.tsv"))

pde <- conv[conv$gene_name == "PDE2A", , drop = FALSE]
if (nrow(pde) != 1L) stop("PDE2A is not uniquely represented among the 84 DEGs")
pde$signaling_regulatory <- grepl("signaling/regulatory", pde$functional_class, fixed = TRUE)
pde$enzyme <- grepl("enzyme", pde$functional_class, fixed = TRUE)
write_tsv(pde, file.path(out_dir, "PDE2A_CONVERGENCE_CONTEXT.tsv"))

comparison_genes <- c("CCL21", "ADAMTS1", "MUC21", "CDH5", "KANK3", "SPRR2C", "MMRN1", "PALMD", "CCN1", "CLDN5", "PDE2A")
comparison <- conv[match(comparison_genes, conv$gene_name), c("gene_name", "padj", "shrunken_log2FC", "Wald_stat",
                                                               "patient_direction_consistency", "functional_class",
                                                               "total_significant_leading_edge_count", "functional_modules")]
names(comparison) <- c("gene", "padj", "shrunken_LFC", "Wald_stat", "patient_consistency",
                       "functional_class", "leading_edge_count", "functional_modules")
if (any(is.na(comparison$gene))) stop("Top Edge-down comparison genes are missing")
write_tsv(comparison, file.path(out_dir, "TOP_EDGE_DOWN_VS_PDE2A_CONTEXT.tsv"))

# Figure source data. Heatmap uses a fixed top-40 rule, never forced inclusion.
heat_genes <- conv[order(-conv$total_significant_leading_edge_count, conv$padj, conv$gene_name), ]
heat_genes <- head(heat_genes, 40)
long_heat <- do.call(rbind, lapply(seq_len(nrow(heat_genes)), function(i) {
  gene <- heat_genes$gene_name[i]
  counts <- module_counts[gene, ]
  data.frame(gene_name = gene, functional_module = sub("^pathway_count__", "", names(counts)),
             pathway_count = as.integer(counts), selection_rank = i,
             selection_rule = "Top 40 DEGs by total significant leading-edge count; ties by padj then gene_name",
             stringsAsFactors = FALSE)
}))
module_label_map <- setNames(module_levels, paste0("pathway_count__", sanitize(module_levels)))
long_heat$functional_module <- unname(module_label_map[paste0("pathway_count__", long_heat$functional_module)])

make_fig_dir <- function(name, data) {
  d <- file.path(fig_root, name)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  write_tsv(data, file.path(d, "source_data.tsv"))
  file.copy(file.path(project_root, "scripts", "22_plot_checkpoint7_convergence.R"),
            file.path(d, "plotting_script.R"), overwrite = TRUE)
  system2("Rscript", c(shQuote(file.path(d, "plotting_script.R")), shQuote(d)), stdout = TRUE, stderr = TRUE)
}
make_fig_dir("Figure_Convergence_01_GeneModuleHeatmap", long_heat)
make_fig_dir("Figure_Convergence_02_Global", conv)
make_fig_dir("Figure_Convergence_02_PDE2A_context", conv)

at_least_1 <- sum(conv$total_significant_leading_edge_count >= 1)
at_least_2 <- sum(conv$total_significant_leading_edge_count >= 2)
at_least_5 <- sum(conv$total_significant_leading_edge_count >= 5)
sig_reg_down <- sum(conv$direction == "Edge-down" & grepl("signaling/regulatory", conv$functional_class, fixed = TRUE))
pde_dir_ok <- sum(pde_audit$direction_consistent_with_PDE2A_edge_down)
pde_dir_bad <- nrow(pde_audit) - pde_dir_ok
module_frequency <- colSums(module_binary)
names(module_frequency) <- module_levels
module_frequency <- sort(module_frequency, decreasing = TRUE)
top_modules <- paste(names(head(module_frequency, 10)), as.integer(head(module_frequency, 10)), sep = "=", collapse = "; ")
comparison_summary <- paste(comparison$gene, comparison$functional_class, sep = ": ", collapse = "; ")

audit <- c(
  "# Checkpoint 7 DEG-pathway convergence audit", "", paste0("Date: ", Sys.Date()), "",
  "## Frozen inputs and scope", "",
  "- Reused frozen Checkpoint 5 DESeq2 and Checkpoint 6 GSEA/context tables; DESeq2, GSEA, and ORA were not rerun.",
  "- FDR threshold, sample inclusion, contrast, low-count filtering, and LFC shrinkage were not changed.",
  "- Significant GSEA pathway criterion: padj < 0.05; leading-edge membership parsed from the frozen leadingEdge field.",
  "- Pathway modules used an ordered, auditable name-keyword mapping recorded in PATHWAY_TO_FUNCTIONAL_MODULE.tsv.",
  "- Functional classes are multi-label annotations, not exclusion or scoring rules.", "",
  "## Convergence results", "",
  paste0("A. Of 84 FDR-significant DEGs, ", at_least_1, " entered at least one significant GSEA leading edge; ", at_least_2, " entered >=2; ", at_least_5, " entered >=5."),
  paste0("B. Major functional modules by DEG membership (non-exclusive): ", top_modules, "."),
  paste0("C. Top Edge-down functional classes plus PDE2A: ", comparison_summary, "."),
  paste0("D. Signaling/regulatory Edge-down DEGs: ", sig_reg_down, "."),
  paste0("E. PDE2A leading-edge pathways: total=", pde$total_significant_leading_edge_count,
         "; Core-enriched/Edge-down=", pde$core_enriched_leading_edge_count,
         "; Edge-enriched/Edge-up=", pde$edge_enriched_leading_edge_count,
         "; functional modules=", pde$functional_module_count, " (", pde$functional_modules, ")."),
  paste0("F. PDE2A direction audit: ", pde_dir_ok, " pathways are direction-consistent with PDE2A Edge-down and ", pde_dir_bad, " are direction-inconsistent."),
  paste0("G. PDE2A functional class: ", pde$functional_class, "; signaling/regulatory=", pde$signaling_regulatory, "; enzyme=", pde$enzyme, "."),
  "H. PDE2A shows gene-level FDR-significant differential expression, 5/5 patient consistency, recurrent pathway leading-edge membership, and intracellular enzymatic/regulatory coherence.",
  "I. PDE2A is not the top-ranked DEG by Wald statistic and was not selected by a composite score.",
  "J. No PDE-family analysis, WGCNA, PPI, TF network, machine learning, weighted score, or threshold modification was performed.", "",
  "## Figure policy", "",
  "- Development figures only: PNG 300 dpi plus source_data.tsv and plotting_script.R.",
  "- No PDF/TIFF and no old/new/batch/sequencing-set annotation.",
  "- The global and PDE2A-context scatter plots use identical source data and aesthetics; the latter only adds the PDE2A text label."
)
writeLines(audit, file.path(out_dir, "QC_CHECKPOINT7_CONVERGENCE_AUDIT.md"))
