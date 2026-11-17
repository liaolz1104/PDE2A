#!/usr/bin/env Rscript

# Checkpoint 10: paired transcriptome-inferred pathway and TF activity.
# Global activity results are written and frozen before any PDE2A-context extraction.

suppressPackageStartupMessages({
  library(DESeq2)
  library(decoupleR)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 10_regulatory_activity_inference.R <project_root>")
root <- normalizePath(args[[1]], mustWork = TRUE)
out <- file.path(root, "results", "10_regulatory_activity")
fig_dir <- file.path(root, "results", "07_figures", "regulatory_activity")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

expected_samples <- c("P08_C", "P08_TI", "P09_C", "P09_TI", "P34_C", "P34_TI",
                      "P36_C", "P36_TI", "P38_C", "P38_TI")
dds_file <- file.path(root, "results", "05_deseq2", "DESeq2_paired_dds.rds")
ann_file <- file.path(root, "results", "04_featurecounts", "gene_annotation.tsv")
progeny_file <- file.path(out, "resources", "PROGENy_human_top500_resource.rds")
collectri_file <- file.path(out, "resources", "CollecTRI_human_resource.rds")
resource_meta_file <- file.path(out, "resources", "ACTIVITY_RESOURCE_METADATA.tsv")
for (p in c(dds_file, ann_file, progeny_file, collectri_file, resource_meta_file)) {
  if (!file.exists(p)) stop("Required frozen input is missing: ", p)
}

# 1. Frozen 10-sample VST input.
dds <- readRDS(dds_file)
if (!identical(colnames(dds), expected_samples)) stop("DDS sample order differs from the frozen 10-sample order.")
if (ncol(dds) != 10L || nrow(dds) != 20413L) stop("Unexpected DDS dimensions: ", nrow(dds), " x ", ncol(dds))
if (paste(deparse(design(dds)), collapse = "") != "~patient_id + region") stop("Unexpected DDS design.")

meta <- as.data.frame(colData(dds)) %>%
  transmute(sample_id = as.character(sample_id), patient_id = factor(patient_id),
            region = factor(region, levels = c("c", "ti")))
rownames(meta) <- meta$sample_id
if (anyNA(meta$patient_id) || anyNA(meta$region)) stop("Missing patient or region metadata.")
if (!all(table(meta$patient_id, meta$region) == 1L)) stop("Each patient must have exactly one Core and one Edge sample.")

vsd <- vst(dds, blind = FALSE)
vst_gene_id <- assay(vsd)
if (!identical(colnames(vst_gene_id), expected_samples) || anyDuplicated(rownames(vst_gene_id))) {
  stop("VST matrix failed sample-order or gene-ID uniqueness checks.")
}
write_tsv(as.data.frame(vst_gene_id) %>% tibble::rownames_to_column("gene_id"),
          file.path(out, "INPUT_VST_10SAMPLES.tsv"))
saveRDS(vsd, file.path(out, "VST_10SAMPLES_blindFALSE.rds"))

# 2. GENCODE v44 gene-ID mapping and deterministic duplicate-symbol collapse.
ann <- read_tsv(ann_file, show_col_types = FALSE)
mapping <- tibble(gene_id = rownames(vst_gene_id), mean_VST = rowMeans(vst_gene_id)) %>%
  left_join(ann %>% select(gene_id, gene_id_stripped, gene_name, gene_type), by = "gene_id") %>%
  mutate(
    clear_symbol = !is.na(gene_name) & gene_name != "",
    usable_for_activity = clear_symbol & gene_type == "protein_coding",
    duplicate_symbol = usable_for_activity & (duplicated(gene_name) | duplicated(gene_name, fromLast = TRUE))
  ) %>%
  group_by(gene_name) %>%
  arrange(desc(mean_VST), gene_id, .by_group = TRUE) %>%
  mutate(retained_yes_no = usable_for_activity & row_number() == 1L) %>%
  ungroup() %>%
  mutate(
    mapping_status = case_when(
      !clear_symbol ~ "UNMAPPED_SYMBOL",
      gene_type != "protein_coding" ~ "NON_PROTEIN_CODING_EXCLUDED",
      duplicate_symbol & retained_yes_no ~ "DUPLICATE_SYMBOL_RETAINED",
      duplicate_symbol & !retained_yes_no ~ "DUPLICATE_SYMBOL_COLLAPSED",
      retained_yes_no ~ "UNIQUE_SYMBOL_RETAINED",
      TRUE ~ "EXCLUDED"
    ),
    retention_rule = case_when(
      duplicate_symbol & retained_yes_no ~ "Highest mean VST; gene_id tie-breaker",
      duplicate_symbol & !retained_yes_no ~ "Collapsed to duplicate symbol entry with highest mean VST",
      retained_yes_no ~ "Unique protein-coding HGNC symbol",
      TRUE ~ "Not eligible for footprint activity"
    )
  )
write_tsv(mapping %>%
            select(gene_id, gene_id_stripped, gene_name, gene_type, mean_VST, mapping_status,
                   duplicate_symbol, retained_yes_no, retention_rule),
          file.path(out, "ACTIVITY_GENE_MAPPING_AUDIT.tsv"))

retained <- filter(mapping, retained_yes_no)
vst_symbol <- vst_gene_id[retained$gene_id, , drop = FALSE]
rownames(vst_symbol) <- retained$gene_name
if (anyDuplicated(rownames(vst_symbol))) stop("HGNC symbol collapse did not produce unique rows.")
write_tsv(as.data.frame(vst_symbol) %>% tibble::rownames_to_column("gene_name"),
          file.path(out, "INPUT_VST_HGNC_COLLAPSED.tsv"))

# 3. Frozen footprint resources and coverage.
progeny_net <- readRDS(progeny_file)
collectri_net <- readRDS(collectri_file)
if (!identical(names(progeny_net), c("source", "target", "weight", "p_value"))) stop("Unexpected PROGENy resource schema.")
if (!identical(names(collectri_net), c("source", "target", "mor"))) stop("Unexpected CollecTRI resource schema.")
coverage <- bind_rows(
  progeny_net %>% group_by(source) %>% summarise(resource = "PROGENy", resource_target_n = n_distinct(target),
                                                      matched_target_n = n_distinct(target[target %in% rownames(vst_symbol)]), .groups = "drop"),
  collectri_net %>% group_by(source) %>% summarise(resource = "CollecTRI", resource_target_n = n_distinct(target),
                                                       matched_target_n = n_distinct(target[target %in% rownames(vst_symbol)]), .groups = "drop")
) %>% mutate(coverage_fraction = matched_target_n / resource_target_n)
write_tsv(coverage, file.path(out, "ACTIVITY_RESOURCE_TARGET_COVERAGE.tsv"))

# decoupleR ULM returns one activity score per source and sample.
progeny_long <- run_ulm(vst_symbol, progeny_net, .source = "source", .target = "target", .mor = "weight", minsize = 5)
tf_long <- run_ulm(vst_symbol, collectri_net, .source = "source", .target = "target", .mor = "mor", minsize = 5)

activity_matrix <- function(x) {
  x %>% select(source, condition, score) %>%
    pivot_wider(names_from = condition, values_from = score) %>%
    arrange(source)
}
progeny_mat_df <- activity_matrix(progeny_long)
tf_mat_df <- activity_matrix(tf_long)
write_tsv(progeny_mat_df %>% rename(pathway = source), file.path(out, "PROGENY_ACTIVITY_MATRIX.tsv"))
write_tsv(tf_mat_df %>% rename(TF = source), file.path(out, "TF_ACTIVITY_MATRIX.tsv"))

progeny_mat <- as.matrix(progeny_mat_df[, expected_samples]); rownames(progeny_mat) <- progeny_mat_df$source
tf_mat <- as.matrix(tf_mat_df[, expected_samples]); rownames(tf_mat) <- tf_mat_df$source

# 4. Paired limma model, Edge/TI minus Core/C.
design_mat <- model.matrix(~ patient_id + region, data = meta)
if (!"regionti" %in% colnames(design_mat)) stop("regionti coefficient absent from paired design matrix.")

paired_limma <- function(mat, feature_name) {
  fit <- eBayes(lmFit(mat, design_mat))
  tt <- topTable(fit, coef = "regionti", number = Inf, sort.by = "none")
  se <- fit$stdev.unscaled[, "regionti"] * sqrt(fit$s2.post)
  core_ids <- rownames(meta)[meta$region == "c"]
  edge_ids <- rownames(meta)[meta$region == "ti"]
  out_tbl <- tibble(
    feature = rownames(tt),
    mean_Core = rowMeans(mat[, core_ids, drop = FALSE]),
    mean_Edge = rowMeans(mat[, edge_ids, drop = FALSE]),
    delta_Edge_minus_Core = tt$logFC,
    SE = se,
    CI95_low = tt$logFC - qt(0.975, df = fit$df.total) * se,
    CI95_high = tt$logFC + qt(0.975, df = fit$df.total) * se,
    moderated_t = tt$t,
    pvalue = tt$P.Value,
    BH_FDR = tt$adj.P.Val
  )
  deltas <- lapply(levels(meta$patient_id), function(pid) {
    c_id <- rownames(meta)[meta$patient_id == pid & meta$region == "c"]
    e_id <- rownames(meta)[meta$patient_id == pid & meta$region == "ti"]
    tibble(feature = rownames(mat), patient_id = pid, Core = mat[, c_id], Edge = mat[, e_id],
           delta_Edge_minus_Core = mat[, e_id] - mat[, c_id])
  }) %>% bind_rows()
  consistency <- deltas %>% group_by(feature) %>% summarise(
    direction_consistency_Edge_higher = sum(delta_Edge_minus_Core > 0),
    direction_consistency_Edge_lower = sum(delta_Edge_minus_Core < 0), .groups = "drop")
  list(result = out_tbl %>% left_join(consistency, by = "feature") %>%
         arrange(BH_FDR, desc(abs(moderated_t))), deltas = deltas)
}

progeny_fit <- paired_limma(progeny_mat, "pathway")
tf_fit <- paired_limma(tf_mat, "TF")
progeny_res <- progeny_fit$result %>% rename(pathway = feature)
tf_res <- tf_fit$result %>% rename(TF = feature)
progeny_deltas <- progeny_fit$deltas %>% rename(pathway = feature)
tf_deltas <- tf_fit$deltas %>% rename(TF = feature)
write_tsv(progeny_res, file.path(out, "PROGENY_EDGE_VS_CORE_LIMMA.tsv"))
write_tsv(progeny_deltas, file.path(out, "PROGENY_PATIENTWISE_DELTAS.tsv"))
write_tsv(tf_res, file.path(out, "TF_ACTIVITY_EDGE_VS_CORE_LIMMA.tsv"))
write_tsv(tf_deltas, file.path(out, "TF_ACTIVITY_PATIENTWISE_DELTAS.tsv"))

# 5. Freeze global results before contextual extraction.
top_description <- function(tbl, name_col, direction, n = 10L) {
  x <- if (direction == "up") filter(tbl, delta_Edge_minus_Core > 0) else filter(tbl, delta_Edge_minus_Core < 0)
  x <- x %>% arrange(BH_FDR, desc(abs(moderated_t))) %>% slice_head(n = n)
  paste(sprintf("- %s: delta=%.3f, t=%.3f, P=%.4g, FDR=%.4g, consistency=%d/5%s",
                x[[name_col]], x$delta_Edge_minus_Core, x$moderated_t, x$pvalue, x$BH_FDR,
                if (direction == "up") x$direction_consistency_Edge_higher else x$direction_consistency_Edge_lower,
                ifelse(x$BH_FDR < 0.05, "", " [NOT FDR SIGNIFICANT]")), collapse = "\n")
}
freeze_lines <- c(
  "# Global activity result freeze",
  "",
  paste0("Frozen before PDE2A contextual extraction: ", Sys.time()),
  "",
  "Sorting rule: BH FDR ascending; tie-breaker absolute moderated t descending.",
  "",
  "## Top Edge-increased PROGENy pathways",
  "", top_description(progeny_res, "pathway", "up", 14L),
  "", "## Top Edge-decreased PROGENy pathways",
  "", top_description(progeny_res, "pathway", "down", 14L),
  "", "## Top Edge-increased TF activities",
  "", top_description(tf_res, "TF", "up", 10L),
  "", "## Top Edge-decreased TF activities",
  "", top_description(tf_res, "TF", "down", 10L)
)
writeLines(freeze_lines, file.path(out, "GLOBAL_ACTIVITY_RESULT_FREEZE.md"))

# 6. Global figures. Development outputs: PNG 300 dpi plus source data and script.
fig1_data <- progeny_res
write_tsv(fig1_data, file.path(fig_dir, "Figure_Activity_01_PROGENy_Global_source_data.tsv"))
p1 <- ggplot(fig1_data, aes(x = delta_Edge_minus_Core, y = reorder(pathway, delta_Edge_minus_Core))) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50", linewidth = 0.45) +
  geom_errorbar(aes(xmin = CI95_low, xmax = CI95_high), orientation = "y", width = 0.16,
                linewidth = 0.5, color = "grey45") +
  geom_point(aes(color = BH_FDR < 0.05), size = 2.6) +
  scale_color_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#4D4D4D"), labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "Not FDR significant")) +
  labs(x = "Inferred activity difference (Edge - Core)", y = NULL, color = NULL,
       title = "Global PROGENy pathway activity") +
  theme_classic(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "Figure_Activity_01_PROGENy_Global.png"), p1, width = 7.0, height = 5.4, dpi = 300)

heat_data <- as.data.frame(progeny_mat) %>% tibble::rownames_to_column("pathway") %>%
  pivot_longer(-pathway, names_to = "sample_id", values_to = "activity") %>%
  group_by(pathway) %>% mutate(row_zscore = as.numeric(scale(activity))) %>% ungroup() %>%
  left_join(meta %>% tibble::rownames_to_column("row_id") %>% select(sample_id, patient_id, region), by = "sample_id")
write_tsv(heat_data, file.path(fig_dir, "Figure_Activity_02_PROGENy_PairedHeatmap_source_data.tsv"))
heat_mat <- t(apply(progeny_mat, 1, scale)); dimnames(heat_mat) <- dimnames(progeny_mat)
ann_col <- data.frame(Patient = meta$patient_id, Region = ifelse(meta$region == "c", "Core", "Edge"), row.names = rownames(meta))
png(file.path(fig_dir, "Figure_Activity_02_PROGENy_PairedHeatmap.png"), width = 2400, height = 1800, res = 300)
pheatmap(heat_mat, cluster_rows = TRUE, cluster_cols = FALSE, annotation_col = ann_col,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101), border_color = NA,
         main = "PROGENy inferred activity across paired samples", fontsize = 8)
dev.off()

fig3_data <- tf_res %>% arrange(desc(abs(moderated_t))) %>% slice_head(n = 20L)
write_tsv(fig3_data, file.path(fig_dir, "Figure_Activity_03_TF_Global_source_data.tsv"))
p3 <- ggplot(fig3_data, aes(x = delta_Edge_minus_Core, y = reorder(TF, delta_Edge_minus_Core))) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50", linewidth = 0.45) +
  geom_errorbar(aes(xmin = CI95_low, xmax = CI95_high), orientation = "y", width = 0.14,
                linewidth = 0.45, color = "grey50") +
  geom_point(aes(color = BH_FDR < 0.05), size = 2.4) +
  scale_color_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#4D4D4D"), labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "Not FDR significant")) +
  labs(x = "Inferred TF-activity difference (Edge - Core)", y = NULL, color = NULL,
       title = "Top 20 TF activities by absolute moderated t") +
  theme_classic(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "Figure_Activity_03_TF_Global.png"), p3, width = 7.2, height = 6.2, dpi = 300)

# 7. PDE2A context extraction after global tables and freeze are written.
relevant_pathways <- c("MAPK", "PI3K", "EGFR", "TGFb", "VEGF", "JAK-STAT", "WNT")
pde_path <- tibble(pathway_requested = relevant_pathways) %>%
  left_join(progeny_res, by = c("pathway_requested" = "pathway")) %>%
  mutate(present = !is.na(delta_Edge_minus_Core), interpretation = case_when(
    !present ~ "NOT PRESENT IN FROZEN PROGENy RESOURCE",
    BH_FDR < 0.05 ~ "FDR SIGNIFICANT",
    TRUE ~ "NOT FDR SIGNIFICANT"
  ))
write_tsv(pde_path, file.path(out, "PDE2A_RELEVANT_PROGENY_CONTEXT.tsv"))

tf_names <- tf_res$TF
requested_tf <- unique(c(
  "CREB1", "ATF1", "JUN", "FOS", "ELK1", "MYC",
  tf_names[grepl("^(ETS|ELK)", tf_names)],
  intersect(tf_names, c("ERG", "FLI1", "FEV", "JUNB", "JUND", "FOSB", "FOSL1", "FOSL2", "ATF2", "ATF3", "BATF", "BATF3", "AP1"))
))
pde_tf <- tibble(TF_requested = requested_tf) %>%
  left_join(tf_res, by = c("TF_requested" = "TF")) %>%
  mutate(present = !is.na(delta_Edge_minus_Core), interpretation = case_when(
    !present ~ "NOT COMPUTABLE WITH FROZEN CollecTRI/ULM",
    BH_FDR < 0.05 ~ "FDR SIGNIFICANT",
    TRUE ~ "NOT FDR SIGNIFICANT"
  ))
write_tsv(pde_tf, file.path(out, "PDE2A_RELEVANT_TF_CONTEXT.tsv"))

# Figure 4 is generated only if at least one prelisted relevant pathway or TF is FDR significant.
if (any(pde_path$BH_FDR < 0.05, na.rm = TRUE) || any(pde_tf$BH_FDR < 0.05, na.rm = TRUE)) {
  context_data <- bind_rows(
    pde_path %>% filter(present) %>% transmute(feature = pathway_requested, type = "PROGENy pathway", delta = delta_Edge_minus_Core, FDR = BH_FDR),
    pde_tf %>% filter(present, TF_requested %in% c("CREB1", "ATF1", "JUN", "FOS", "ELK1", "MYC", "AP1")) %>%
      transmute(feature = TF_requested, type = "TF activity", delta = delta_Edge_minus_Core, FDR = BH_FDR)
  )
  write_tsv(context_data, file.path(fig_dir, "Figure_Activity_04_PDE2A_MechanisticContext_source_data.tsv"))
  p4 <- ggplot(context_data, aes(x = delta, y = reorder(feature, delta), color = type, shape = FDR < 0.05)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey50", linewidth = 0.45) + geom_point(size = 2.8) +
    labs(x = "Transcriptome-inferred activity difference (Edge - Core)", y = NULL,
         title = "PDE2A-associated regulatory context (non-causal)", color = NULL, shape = "FDR < 0.05") +
    theme_classic(base_size = 10) + theme(legend.position = "bottom")
  ggsave(file.path(fig_dir, "Figure_Activity_04_PDE2A_MechanisticContext.png"), p4, width = 7.0, height = 5.4, dpi = 300)
}

# Machine-readable run metadata; audit prose is generated after result verification.
resource_meta <- read_tsv(resource_meta_file, show_col_types = FALSE)
run_meta <- tibble(
  field = c("run_date", "R_version", "DESeq2_version", "limma_version", "decoupleR_version", "OmnipathR_version",
            "dds_input", "dds_dimensions", "VST_blind", "activity_gene_n", "PROGENy_pathway_n", "TF_n",
            "design", "contrast", "GSE178537_used", "optional_second_messenger_scores"),
  value = c(as.character(Sys.time()), R.version.string, as.character(packageVersion("DESeq2")),
            as.character(packageVersion("limma")), as.character(packageVersion("decoupleR")),
            as.character(packageVersion("OmnipathR")), dds_file, paste(dim(vst_gene_id), collapse = "x"), "FALSE",
            nrow(vst_symbol), nrow(progeny_res), nrow(tf_res), "~ patient_id + region", "Edge_minus_Core",
            "FALSE", "SKIPPED: optional and not independent evidence")
)
write_tsv(run_meta, file.path(out, "CHECKPOINT10_RUN_METADATA.tsv"))

cat("Checkpoint 10 global activity inference completed.\n")
cat("PROGENy pathways:", nrow(progeny_res), "TFs:", nrow(tf_res), "\n")
