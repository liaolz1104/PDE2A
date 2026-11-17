#!/usr/bin/env Rscript

# Checkpoint 9: candidate-level replication in GSE178537.
# Frozen direction: VT (invasive front/Edge) minus IT (inner tumor/Core).
# This script stops before statistics if any provenance assertion fails.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: 09_gse178537_candidate_replication.R <project_root> <rld_matrix.gz> <official_geo_soft.txt> <corrected_metadata.csv>")
}

project_root <- normalizePath(args[[1]], mustWork = TRUE)
rld_file <- normalizePath(args[[2]], mustWork = TRUE)
soft_file <- normalizePath(args[[3]], mustWork = TRUE)
metadata_file <- normalizePath(args[[4]], mustWork = TRUE)
out_dir <- file.path(project_root, "results", "09_external_replication")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

panel_file <- file.path(out_dir, "FROZEN_CANDIDATE_PANEL.tsv")
if (!file.exists(panel_file)) stop("Frozen candidate panel is missing: ", panel_file)
panel <- read_tsv(panel_file, show_col_types = FALSE)
expected_genes <- c("ADCY4", "KDR", "CAMK2B", "GNG11", "PDE1A", "PDE2A", "NR4A1")
stopifnot(nrow(panel) == 7L, identical(panel$gene, expected_genes))

# Parse official GEO SOFT sample metadata.
soft_lines <- readLines(soft_file, warn = FALSE)
records <- list()
current <- NULL
for (line in soft_lines) {
  if (startsWith(line, "^SAMPLE = ")) {
    current <- sub("^\\^SAMPLE = ", "", line)
    records[[current]] <- list(gsm = current, title = NA_character_, source = NA_character_)
  } else if (!is.null(current) && startsWith(line, "!Sample_title = ")) {
    records[[current]]$title <- sub("^!Sample_title = ", "", line)
  } else if (!is.null(current) && startsWith(line, "!Sample_source_name_ch1 = ")) {
    records[[current]]$source <- sub("^!Sample_source_name_ch1 = ", "", line)
  }
}
official <- bind_rows(lapply(records, as.data.frame, stringsAsFactors = FALSE))
if (nrow(official) != 65L) stop("Expected 65 GEO sample records, found ", nrow(official))

official <- official %>%
  mutate(
    tag = sub("^[^-]+-", "", title),
    patient_id = sub("-.*$", "", title),
    region = case_when(
      tag == "VT" ~ "Edge",
      tag == "IT" ~ "Core",
      TRUE ~ NA_character_
    )
  )

meta <- read_csv(metadata_file, show_col_types = FALSE) %>%
  mutate(
    patient = as.character(patient),
    expected_region = case_when(tag == "VT" ~ "Edge (VT)", tag == "IT" ~ "Core (IT)", TRUE ~ region)
  )

paired_meta <- official %>% filter(tag %in% c("VT", "IT"))
pair_counts <- paired_meta %>% count(patient_id, tag) %>% pivot_wider(names_from = tag, values_from = n, values_fill = 0)
complete_pairs <- pair_counts %>% filter(VT == 1L, IT == 1L)
expected_pair_ids <- c("1", "6", "20", "22", "26", "44", "47", "54", "56", "58", "59", "61", "63", "64", "65")

assertions <- tibble(
  assertion_id = 1:6,
  assertion = c(
    "GSM5393476 / 1-VT maps to Edge",
    "GSM5393477 / 1-IT maps to Core",
    "All VT samples map to Edge",
    "All IT samples map to Core",
    "Exactly 15 patients each have one VT and one IT",
    "Contrast string is VT_minus_IT"
  ),
  pass = c(
    nrow(filter(official, gsm == "GSM5393476", title == "1-VT", source == "Primary Tumor (Invasive front)", region == "Edge")) == 1L,
    nrow(filter(official, gsm == "GSM5393477", title == "1-IT", source == "Primary Tumor (Inner site)", region == "Core")) == 1L,
    all(filter(official, tag == "VT")$source == "Primary Tumor (Invasive front)") && all(filter(official, tag == "VT")$region == "Edge"),
    all(filter(official, tag == "IT")$source == "Primary Tumor (Inner site)") && all(filter(official, tag == "IT")$region == "Core"),
    nrow(complete_pairs) == 15L && identical(sort(as.integer(complete_pairs$patient_id)), sort(as.integer(expected_pair_ids))),
    identical("VT_minus_IT", "VT_minus_IT")
  ),
  observed = c(
    paste(filter(official, gsm == "GSM5393476")$title, filter(official, gsm == "GSM5393476")$source, sep = "; "),
    paste(filter(official, gsm == "GSM5393477")$title, filter(official, gsm == "GSM5393477")$source, sep = "; "),
    paste(nrow(filter(official, tag == "VT")), "VT samples; all Edge"),
    paste(nrow(filter(official, tag == "IT")), "IT samples; all Core"),
    paste(nrow(complete_pairs), "complete pairs:", paste(sort(as.integer(complete_pairs$patient_id)), collapse = ",")),
    "VT_minus_IT"
  )
)
write_tsv(assertions, file.path(out_dir, "GSE178537_LABEL_ASSERTION_RESULTS.tsv"))
if (!all(assertions$pass)) stop("One or more GSE178537 label assertions failed; statistics were not run.")

# Export the official mapping used by this analysis.
official %>%
  select(gsm, title, patient_id, tag, source, region) %>%
  write_tsv(file.path(out_dir, "GSE178537_OFFICIAL_SAMPLE_METADATA.tsv"))

# Load the official RLD matrix. Gene rows are named ENSG..._SYMBOL; duplicate symbols are averaged.
rld <- read_tsv(rld_file, show_col_types = FALSE)
gene_col <- names(rld)[[1]]
rld <- rld %>%
  rename(gene_raw = all_of(gene_col)) %>%
  mutate(gene = sub("^ENSG[0-9.]+_", "", gene_raw)) %>%
  select(-gene_raw) %>%
  group_by(gene) %>%
  summarise(across(everything(), ~ mean(as.numeric(.x), na.rm = TRUE)), .groups = "drop")

missing_genes <- setdiff(panel$gene, rld$gene)
if (length(missing_genes) > 0L) stop("Candidate genes missing from RLD matrix: ", paste(missing_genes, collapse = ", "))

pair_ids <- as.character(sort(as.integer(complete_pairs$patient_id)))
required_columns <- as.vector(rbind(paste0(pair_ids, "-VT"), paste0(pair_ids, "-IT")))
if (!all(required_columns %in% names(rld))) stop("One or more paired sample columns are absent from RLD matrix.")

expr <- rld %>%
  filter(gene %in% panel$gene) %>%
  pivot_longer(cols = all_of(required_columns), names_to = "sample", values_to = "expression_RLD") %>%
  mutate(
    patient_id = sub("-.*$", "", sample),
    tag = sub("^[^-]+-", "", sample),
    region = if_else(tag == "VT", "Edge", "Core")
  ) %>%
  select(gene, patient_id, sample, tag, region, expression_RLD)

paired <- expr %>%
  select(gene, patient_id, tag, expression_RLD) %>%
  pivot_wider(names_from = tag, values_from = expression_RLD) %>%
  mutate(delta_VT_minus_IT = VT - IT) %>%
  arrange(match(gene, panel$gene), as.integer(patient_id))
write_tsv(paired, file.path(out_dir, "GSE178537_CANDIDATE_PAIRED_SOURCE_DATA.tsv"))

bootstrap_median_ci <- function(delta, seed, reps = 10000L) {
  set.seed(seed)
  boot <- replicate(reps, median(sample(delta, length(delta), replace = TRUE)))
  as.numeric(quantile(boot, probs = c(0.025, 0.975), names = FALSE, type = 7))
}

bootstrap_dz_ci <- function(delta, seed, reps = 10000L) {
  set.seed(seed)
  boot <- replicate(reps, {
    sampled <- sample(delta, length(delta), replace = TRUE)
    sampled_sd <- sd(sampled)
    if (is.finite(sampled_sd) && sampled_sd > 0) mean(sampled) / sampled_sd else NA_real_
  })
  as.numeric(quantile(boot, probs = c(0.025, 0.975), names = FALSE, na.rm = TRUE, type = 7))
}

stats_rows <- lapply(seq_along(panel$gene), function(i) {
  gene_i <- panel$gene[[i]]
  dat <- filter(paired, gene == gene_i)
  delta <- dat$delta_VT_minus_IT
  has_ties_or_zero <- any(delta == 0) || any(duplicated(abs(delta)))
  wt <- wilcox.test(dat$VT, dat$IT, paired = TRUE, alternative = "two.sided", exact = !has_ties_or_zero, correct = has_ties_or_zero)
  ci <- bootstrap_median_ci(delta, seed = 20260829L + i)
  dz_ci <- bootstrap_dz_ci(delta, seed = 20261829L + i)
  tibble(
    gene = gene_i,
    n_pairs = nrow(dat),
    mean_VT_Edge = mean(dat$VT),
    mean_IT_Core = mean(dat$IT),
    mean_delta_VT_minus_IT = mean(delta),
    median_delta_VT_minus_IT = median(delta),
    bootstrap_median_delta_CI_low = ci[[1]],
    bootstrap_median_delta_CI_high = ci[[2]],
    paired_standardized_effect_dz = ifelse(sd(delta) > 0, mean(delta) / sd(delta), NA_real_),
    bootstrap_dz_CI_low = dz_ci[[1]],
    bootstrap_dz_CI_high = dz_ci[[2]],
    direction_consistency_VT_lower = sum(dat$VT < dat$IT),
    direction_consistency_VT_higher = sum(dat$VT > dat$IT),
    primary_paired_test = ifelse(has_ties_or_zero, "Two-sided paired Wilcoxon signed-rank (normal approximation with continuity correction)", "Two-sided exact paired Wilcoxon signed-rank"),
    wilcoxon_statistic = unname(wt$statistic),
    primary_pvalue = wt$p.value,
    direction_replicated = median(delta) < 0,
    patient_majority_replication = sum(dat$VT < dat$IT) > nrow(dat) / 2
  )
})

replication <- bind_rows(stats_rows) %>%
  mutate(candidate_panel_BH_FDR = p.adjust(primary_pvalue, method = "BH")) %>%
  left_join(panel, by = "gene") %>%
  arrange(panel_order) %>%
  mutate(
    external_absolute_effect_rank_among_7 = rank(-abs(paired_standardized_effect_dz), ties.method = "min"),
    external_pvalue_rank_among_7 = rank(primary_pvalue, ties.method = "min"),
    contrast = "VT_minus_IT",
    negative_effect_definition = "Edge_lower_than_Core"
  )
write_tsv(replication, file.path(out_dir, "GSE178537_CANDIDATE_REPLICATION.tsv"))

concordance <- replication %>%
  transmute(
    gene,
    current_shrunken_LFC = current_shrunken_log2FC,
    current_padj,
    current_5patient_consistency = current_patient_consistency,
    current_leading_edge_count,
    current_module_count = current_functional_module_count,
    current_PPI_degree,
    GSE178537_paired_effect_dz = paired_standardized_effect_dz,
    GSE178537_median_delta_VT_minus_IT = median_delta_VT_minus_IT,
    GSE178537_paired_p = primary_pvalue,
    GSE178537_candidate_panel_BH_FDR = candidate_panel_BH_FDR,
    GSE178537_direction_consistency_VT_lower = direction_consistency_VT_lower,
    direction_replicated,
    current_effect_rank_among_7 = rank(-abs(current_shrunken_LFC), ties.method = "min"),
    current_leading_edge_rank_among_7 = rank(-current_leading_edge_count, ties.method = "min"),
    external_effect_rank_among_7 = external_absolute_effect_rank_among_7,
    external_p_rank_among_7 = external_pvalue_rank_among_7
  )
write_tsv(concordance, file.path(out_dir, "CROSS_COHORT_CANDIDATE_CONCORDANCE.tsv"))
write_tsv(filter(replication, gene == "PDE2A"), file.path(out_dir, "PDE2A_GSE178537_REPLICATION.tsv"))
write_tsv(replication, file.path(out_dir, "NATURAL_CANDIDATES_CROSS_COHORT_COMPARISON.tsv"))

# Figure 1: external effects, retaining frozen panel order.
fig1_data <- replication %>%
  select(panel_order, gene, paired_standardized_effect_dz, median_delta_VT_minus_IT,
         bootstrap_median_delta_CI_low, bootstrap_median_delta_CI_high,
         bootstrap_dz_CI_low, bootstrap_dz_CI_high,
         primary_pvalue, candidate_panel_BH_FDR, direction_consistency_VT_lower)
write_tsv(fig1_data, file.path(out_dir, "Figure_ExternalReplication_01_CandidateEffects_source_data.tsv"))
p1 <- ggplot(fig1_data, aes(x = paired_standardized_effect_dz, y = factor(gene, levels = rev(panel$gene)))) +
  geom_vline(xintercept = 0, color = "grey45", linewidth = 0.5, linetype = 2) +
  geom_errorbar(aes(xmin = bootstrap_dz_CI_low, xmax = bootstrap_dz_CI_high), orientation = "y",
                width = 0.16, linewidth = 0.55, color = "#355C7D") +
  geom_point(size = 2.8, color = "#355C7D") +
  labs(x = "Paired standardized effect (VT/Edge - IT/Core)", y = NULL,
       title = "GSE178537 candidate-level paired effects") +
  theme_classic(base_size = 10)
ggsave(file.path(out_dir, "Figure_ExternalReplication_01_CandidateEffects.png"), p1, width = 6.4, height = 4.2, dpi = 300)

# Figure 2: direction/effect table plot. Cross-cohort standardized effects are not forced
# because only the external RLD paired values are available in a directly reproducible form here.
fig2_data <- concordance %>%
  select(gene, current_shrunken_LFC, GSE178537_paired_effect_dz, direction_replicated,
         GSE178537_direction_consistency_VT_lower)
write_tsv(fig2_data, file.path(out_dir, "Figure_ExternalReplication_02_DirectionEffectTable_source_data.tsv"))
fig2_long <- fig2_data %>%
  select(gene, current_shrunken_LFC, GSE178537_paired_effect_dz) %>%
  pivot_longer(-gene, names_to = "cohort_metric", values_to = "effect") %>%
  mutate(cohort_metric = recode(cohort_metric,
    current_shrunken_LFC = "Current: shrunken log2FC",
    GSE178537_paired_effect_dz = "GSE178537: paired effect dz"))
p2 <- ggplot(fig2_long, aes(x = effect, y = factor(gene, levels = rev(panel$gene)), color = cohort_metric)) +
  geom_vline(xintercept = 0, color = "grey45", linewidth = 0.5, linetype = 2) +
  geom_point(size = 2.5, position = position_dodge(width = 0.45)) +
  scale_color_manual(values = c("Current: shrunken log2FC" = "#7C3AED", "GSE178537: paired effect dz" = "#008B8B")) +
  labs(x = "Within-cohort effect metric (not a shared absolute scale)", y = NULL, color = NULL,
       title = "Cross-cohort direction and effect summary") +
  theme_classic(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "Figure_ExternalReplication_02_DirectionEffectTable.png"), p2, width = 7.2, height = 4.5, dpi = 300)

# Figure 3: PDE2A paired RLD expression.
pde <- filter(paired, gene == "PDE2A")
pde_long <- pde %>%
  select(patient_id, IT, VT) %>%
  pivot_longer(c(IT, VT), names_to = "tag", values_to = "expression_RLD") %>%
  mutate(region = factor(if_else(tag == "IT", "Core", "Edge"), levels = c("Core", "Edge")))
write_tsv(pde_long, file.path(out_dir, "Figure_ExternalReplication_03_PDE2A_Paired_source_data.tsv"))
pde_stat <- filter(replication, gene == "PDE2A")
p3 <- ggplot(pde_long, aes(x = region, y = expression_RLD, group = patient_id)) +
  geom_line(color = "grey70", linewidth = 0.55) +
  geom_point(aes(color = region), size = 2.3) +
  scale_color_manual(values = c(Core = "#4C78A8", Edge = "#E45756"), guide = "none") +
  labs(x = NULL, y = "PDE2A RLD expression",
       title = "PDE2A independent paired validation",
       subtitle = sprintf("n = 15; median VT-IT = %.3f; two-sided paired P = %.4g",
                          pde_stat$median_delta_VT_minus_IT, pde_stat$primary_pvalue)) +
  theme_classic(base_size = 10)
ggsave(file.path(out_dir, "Figure_ExternalReplication_03_PDE2A_Paired.png"), p3, width = 4.6, height = 4.4, dpi = 300)

# Machine-readable run metadata.
run_meta <- tibble(
  field = c("accession", "input_type", "input_file", "candidate_panel_n", "paired_patient_n",
            "contrast", "primary_test", "multiple_testing_scope", "bootstrap_replicates", "bootstrap_seed_base"),
  value = c("GSE178537", "GEO processed regularized-log/RLD continuous expression", rld_file, "7", "15",
            "VT_minus_IT", "two-sided paired Wilcoxon signed-rank", "BH across frozen 7-gene panel", "10000", "20260829")
)
write_tsv(run_meta, file.path(out_dir, "GSE178537_CHECKPOINT9_RUN_METADATA.tsv"))

cat("Checkpoint 9 analysis completed successfully.\n")
