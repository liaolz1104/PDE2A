#!/usr/bin/env Rscript

# Exploratory gene-set and immune-context sensitivity analyses for the frozen
# five-patient paired Core/Edge RNA-seq cohort. This script does not rerun or
# modify the frozen DESeq2, GSEA, PPI, or candidate-prioritization results.

.libPaths(c(
  "D:/CodexHome/r-libs/pde2a-sensitivity",
  "C:/Users/liaol/AppData/Local/R/win-library/4.5",
  .libPaths()
))

suppressPackageStartupMessages({
  library(limma)
  library(PADOG)
  library(GSReg)
  library(GSVA)
  library(MCPcounter)
  library(tidyestimate)
  library(ggplot2)
  library(ragg)
})

set.seed(20261002)

# Keep the ASCII junction path used on Windows; normalizePath() resolves the
# junction to a non-ASCII target that this R locale cannot decode.
ROOT <- gsub("\\\\", "/", getwd())
OUT <- file.path(ROOT, "results", "12_gene_set_immune_sensitivity")
FIG <- file.path(OUT, "figures")
RES <- file.path(OUT, "resources")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(RES, recursive = TRUE, showWarnings = FALSE)

FONT_FAMILY <- "Arial"
FONT_PT <- 7
COL_CORE <- "#3C5488"
COL_EDGE <- "#E64B35"
COL_NEUTRAL <- "#8A8A8A"

read_tsv <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
             quote = "", comment.char = "", na.strings = c("NA", ""))
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

save_png <- function(filename, plot, width, height) {
  # ragg cannot open a device through this Windows directory junction, so it
  # renders to an ASCII-only temporary path and copies the finished PNG back.
  tmp_dir <- "D:/CodexHome/tmp/pde2a-gene-set-immune-figures"
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  tmp <- file.path(tmp_dir, paste0(sprintf("%08d", sample.int(1e8, 1)), ".png"))
  ggsave(tmp, plot, width = width, height = height, units = "in", dpi = 300,
         device = ragg::agg_png, bg = "white")
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(tmp, filename, overwrite = TRUE)
  unlink(tmp)
  if (!ok || !file.exists(filename) || file.info(filename)$size == 0) {
    stop("Failed to export PNG: ", filename)
  }
  invisible(filename)
}

clean_term <- function(x) {
  y <- gsub("^HALLMARK_", "", x)
  y <- gsub("_", " ", y, fixed = TRUE)
  y <- tools::toTitleCase(tolower(y))
  replacements <- c(
    "G2m Checkpoint" = "G2/M Checkpoint",
    "Tnfa Signaling Via Nfkb" = "TNF-alpha Signaling via NF-kappaB",
    "Interferon Alpha Response" = "Interferon-alpha Response",
    "Interferon Gamma Response" = "Interferon-gamma Response",
    "Tgf Beta Signaling" = "TGF-beta Signaling",
    "Il6 Jak Stat3 Signaling" = "IL-6/JAK/STAT3 Signaling",
    "Pi3k Akt Mtor Signaling" = "PI3K/AKT/mTOR Signaling",
    "Mtorc1 Signaling" = "mTORC1 Signaling",
    "Kras Signaling Dn" = "KRAS Signaling DN",
    "Kras Signaling Up" = "KRAS Signaling UP",
    "Myc Targets V1" = "MYC Targets V1",
    "Myc Targets V2" = "MYC Targets V2",
    "P53 Pathway" = "p53 Pathway",
    "Dna Repair" = "DNA Repair",
    "Uv Response Dn" = "UV Response DN",
    "Uv Response Up" = "UV Response UP"
  )
  hit <- match(y, names(replacements))
  y[!is.na(hit)] <- unname(replacements[hit[!is.na(hit)]])
  y
}

read_gmt <- function(path, universe, min_size = 10L, max_size = 500L) {
  zz <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  gs <- lapply(zz, function(z) unique(intersect(z[-c(1, 2)], universe)))
  names(gs) <- vapply(zz, function(z) z[1], character(1))
  gs[lengths(gs) >= min_size & lengths(gs) <= max_size]
}

paired_meta <- function(samples) {
  data.frame(
    sample = samples,
    patient = factor(sub("_(C|TI)$", "", samples)),
    region = factor(ifelse(grepl("_TI$", samples), "Edge", "Core"),
                    levels = c("Core", "Edge")),
    stringsAsFactors = FALSE
  )
}

paired_limma <- function(score_matrix, meta) {
  stopifnot(identical(colnames(score_matrix), meta$sample))
  design <- model.matrix(~ patient + region, data = meta)
  coef_name <- "regionEdge"
  fit <- eBayes(lmFit(score_matrix, design))
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  tt$SE <- ifelse(is.finite(tt$t) & tt$t != 0, abs(tt$logFC / tt$t), NA_real_)
  tt$feature <- rownames(tt)
  deltas <- do.call(rbind, lapply(seq_len(nrow(score_matrix)), function(i) {
    vals <- score_matrix[i, ]
    names(vals) <- meta$sample
    vapply(levels(meta$patient), function(p) {
      vals[meta$patient == p & meta$region == "Edge"] -
        vals[meta$patient == p & meta$region == "Core"]
    }, numeric(1))
  }))
  rownames(deltas) <- rownames(score_matrix)
  tt$mean_Core <- rowMeans(score_matrix[, meta$region == "Core", drop = FALSE])
  tt$mean_Edge <- rowMeans(score_matrix[, meta$region == "Edge", drop = FALSE])
  tt$delta_Edge_minus_Core <- tt$mean_Edge - tt$mean_Core
  tt$Edge_higher_count <- rowSums(deltas > 0, na.rm = TRUE)
  tt$Edge_lower_count <- rowSums(deltas < 0, na.rm = TRUE)
  tt <- tt[, c("feature", "mean_Core", "mean_Edge", "delta_Edge_minus_Core",
               "logFC", "SE", "t", "P.Value", "adj.P.Val", "B",
               "Edge_higher_count", "Edge_lower_count")]
  colnames(tt) <- c("feature", "mean_Core", "mean_Edge", "delta_Edge_minus_Core",
                    "limma_difference", "SE", "moderated_t", "pvalue", "BH_FDR",
                    "B", "Edge_higher_count", "Edge_lower_count")
  list(stats = tt, deltas = deltas)
}

plot_effect <- function(tab, feature_col, title, xlab, out_dir, n_show = 20L) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  z <- tab[order(tab$BH_FDR, -abs(tab$moderated_t)), , drop = FALSE]
  z <- head(z, n_show)
  z$display <- clean_term(z[[feature_col]])
  z$display <- factor(z$display, levels = rev(z$display))
  z$direction <- ifelse(z$delta_Edge_minus_Core >= 0, "Edge-increased", "Edge-decreased")
  z$lower <- z$delta_Edge_minus_Core - 1.96 * z$SE
  z$upper <- z$delta_Edge_minus_Core + 1.96 * z$SE
  write_tsv(z, file.path(out_dir, "source_data.tsv"))
  p <- ggplot(z, aes(x = delta_Edge_minus_Core, y = display, colour = direction)) +
    geom_vline(xintercept = 0, colour = "#B3B3B3", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, linewidth = 0.35) +
    geom_point(size = 1.8) +
    scale_colour_manual(values = c("Edge-decreased" = COL_CORE, "Edge-increased" = COL_EDGE)) +
    labs(title = title, x = xlab, y = NULL, colour = NULL) +
    theme_classic(base_family = FONT_FAMILY, base_size = FONT_PT) +
    theme(plot.title = element_text(face = "bold", size = FONT_PT),
          legend.position = "top", legend.justification = "left",
          legend.text = element_text(size = FONT_PT),
          axis.text = element_text(size = FONT_PT),
          axis.title = element_text(size = FONT_PT))
  save_png(file.path(out_dir, "figure.png"), p, width = 4.5, height = 4.8)
}

plot_signed_significance <- function(tab, title, out_dir, n_show = 20L) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  z <- tab[order(tab$FDR, tab$PValue), , drop = FALSE]
  z <- head(z, n_show)
  z$display <- clean_term(z$pathway)
  z$display <- factor(z$display, levels = rev(z$display))
  z$signed_log10_FDR <- ifelse(z$Direction == "Up", 1, -1) *
    -log10(pmax(z$FDR, .Machine$double.xmin))
  z$direction <- ifelse(z$Direction == "Up", "Edge-enriched", "Core-enriched")
  write_tsv(z, file.path(out_dir, "source_data.tsv"))
  p <- ggplot(z, aes(x = signed_log10_FDR, y = display, colour = direction)) +
    geom_vline(xintercept = 0, colour = "#B3B3B3", linewidth = 0.3) +
    geom_segment(aes(x = 0, xend = signed_log10_FDR, yend = display), linewidth = 0.5) +
    geom_point(size = 1.8) +
    scale_colour_manual(values = c("Core-enriched" = COL_CORE, "Edge-enriched" = COL_EDGE)) +
    labs(title = title, x = "Signed -log10(FDR)", y = NULL, colour = NULL) +
    theme_classic(base_family = FONT_FAMILY, base_size = FONT_PT) +
    theme(plot.title = element_text(face = "bold", size = FONT_PT),
          legend.position = "top", legend.justification = "left",
          legend.text = element_text(size = FONT_PT),
          axis.text = element_text(size = FONT_PT),
          axis.title = element_text(size = FONT_PT))
  save_png(file.path(out_dir, "figure.png"), p, width = 4.8, height = 4.8)
}

# Frozen VST matrix and paired metadata.
vst_df <- read_tsv(file.path(ROOT, "results", "10_regulatory_activity", "INPUT_VST_HGNC_COLLAPSED.tsv"))
vst <- as.matrix(vst_df[, -1, drop = FALSE])
storage.mode(vst) <- "double"
rownames(vst) <- vst_df[[1]]
meta <- paired_meta(colnames(vst))
stopifnot(nrow(meta) == 10L, length(unique(meta$patient)) == 5L,
          all(table(meta$patient, meta$region) == 1L))

gmt_paths <- c(
  Hallmark = file.path(ROOT, "results", "06_enrichment", "gene_sets", "MSigDB_Hallmark_2020.gmt"),
  Reactome = file.path(ROOT, "results", "06_enrichment", "gene_sets", "Reactome_2022.gmt"),
  GO_BP = file.path(ROOT, "results", "06_enrichment", "gene_sets", "GO_Biological_Process_2023.gmt"),
  KEGG = file.path(ROOT, "results", "06_enrichment", "gene_sets", "KEGG_2021_Human.gmt")
)
gene_sets <- lapply(gmt_paths, read_gmt, universe = rownames(vst))

resource_manifest <- data.frame(
  database = names(gene_sets),
  source_file = unname(gmt_paths),
  tested_gene_sets = vapply(gene_sets, length, integer(1)),
  min_size = 10L,
  max_size = 500L,
  stringsAsFactors = FALSE
)
write_tsv(resource_manifest, file.path(OUT, "GENE_SET_RESOURCE_MANIFEST.tsv"))

# CAMERA and fry: paired limma design, Edge - Core coefficient.
design <- model.matrix(~ patient + region, data = meta)
coef_idx <- which(colnames(design) == "regionEdge")
camera_path <- file.path(OUT, "CAMERA_EDGE_VS_CORE.tsv")
fry_path <- file.path(OUT, "FRY_EDGE_VS_CORE.tsv")
if (file.exists(camera_path) && file.exists(fry_path)) {
  camera_tab <- read_tsv(camera_path)
  fry_tab <- read_tsv(fry_path)
} else {
  camera_all <- list()
  fry_all <- list()
  for (db in names(gene_sets)) {
    idx <- ids2indices(gene_sets[[db]], rownames(vst))
    ca <- camera(vst, idx, design, contrast = coef_idx)
    ca$pathway <- rownames(ca)
    ca$database <- db
    camera_all[[db]] <- ca[, c("database", "pathway", "NGenes", "Direction", "PValue", "FDR")]
    fr <- fry(vst, idx, design, contrast = coef_idx, sort = "mixed")
    fr$pathway <- rownames(fr)
    fr$database <- db
    fry_all[[db]] <- fr[, c("database", "pathway", "NGenes", "Direction", "PValue", "FDR")]
  }
  camera_tab <- do.call(rbind, camera_all)
  fry_tab <- do.call(rbind, fry_all)
  write_tsv(camera_tab, camera_path)
  write_tsv(fry_tab, fry_path)
}
plot_signed_significance(camera_tab, "CAMERA competitive gene-set analysis",
                         file.path(FIG, "Figure_CAMERA_Global"))
plot_signed_significance(fry_tab, "fry self-contained gene-set analysis",
                         file.path(FIG, "Figure_FRY_Global"))

# PADOG: paired design with patient block. Ppadog is permutation P; BH is
# calculated within each frozen database. PADOG is non-directional, so pathway
# direction is descriptive and derived from mean member-gene Edge-Core delta.
group_padog <- ifelse(meta$region == "Edge", "d", "c")
padog_path <- file.path(OUT, "PADOG_PAIRED_EDGE_VS_CORE.tsv")
if (file.exists(padog_path)) {
  padog_tab <- read_tsv(padog_path)
} else {
  padog_all <- list()
  gene_delta <- rowMeans(vst[, meta$region == "Edge", drop = FALSE] -
                           vst[, meta$region == "Core", drop = FALSE])
  for (db in names(gene_sets)) {
    pd <- PADOG::padog(
      esetm = vst, group = group_padog, paired = TRUE,
      block = as.character(meta$patient), gslist = gene_sets[[db]],
      NI = 32, plots = FALSE, Nmin = 10, verbose = FALSE,
      parallel = FALSE, dseed = 20261002
    )
    pd$pathway <- rownames(pd)
    pd$database <- db
    pd$BH_FDR <- p.adjust(pd$Ppadog, method = "BH")
    mean_delta <- vapply(gene_sets[[db]], function(g) mean(gene_delta[g], na.rm = TRUE), numeric(1))
    lookup <- setNames(mean_delta, clean_term(names(mean_delta)))
    pd$mean_member_gene_delta <- unname(lookup[clean_term(pd$pathway)])
    pd$descriptive_direction <- ifelse(pd$mean_member_gene_delta >= 0,
                                       "Edge-increased", "Edge-decreased")
    padog_all[[db]] <- pd[, c("database", "pathway", "Size", "padog0", "Ppadog",
                               "BH_FDR", "mean_member_gene_delta", "descriptive_direction")]
  }
  padog_tab <- do.call(rbind, padog_all)
}

# Five paired subjects yield only 31 non-identity pair-label flips. PADOG
# replaces a zero exceedance count by 1/(NI*100), which is anti-conservative
# for an exact randomization test. Preserve the package value for audit, and
# add the conservative exact permutation P=(b+1)/(NI+1).
padog_tab$package_Ppadog <- padog_tab$Ppadog
padog_tab$package_BH_FDR <- padog_tab$BH_FDR
padog_tab$permutations_used <- 31L
padog_tab$package_zero_exceedance_floor <-
  abs(padog_tab$package_Ppadog - 1 / 31 / 100) < 1e-12
padog_tab$exceedance_count <- ifelse(
  padog_tab$package_zero_exceedance_floor,
  0L,
  as.integer(round(padog_tab$package_Ppadog * 31))
)
padog_tab$conservative_exact_pvalue <- (padog_tab$exceedance_count + 1) / 32
padog_tab$conservative_BH_FDR <- ave(
  padog_tab$conservative_exact_pvalue,
  padog_tab$database,
  FUN = function(z) p.adjust(z, method = "BH")
)
write_tsv(padog_tab, padog_path)

padog_plot <- padog_tab[order(padog_tab$conservative_BH_FDR,
                              padog_tab$conservative_exact_pvalue), , drop = FALSE]
padog_plot <- head(padog_plot, 20)
padog_plot$display <- factor(clean_term(padog_plot$pathway),
                             levels = rev(clean_term(padog_plot$pathway)))
write_tsv(padog_plot, file.path(FIG, "Figure_PADOG_Global", "source_data.tsv"))
p_padog <- ggplot(padog_plot,
                  aes(x = mean_member_gene_delta, y = display,
                      colour = descriptive_direction,
                      size = -log10(pmax(conservative_BH_FDR, .Machine$double.xmin)))) +
  geom_vline(xintercept = 0, colour = "#B3B3B3", linewidth = 0.3) +
  geom_point(alpha = 0.9) +
  scale_colour_manual(values = c("Edge-decreased" = COL_CORE, "Edge-increased" = COL_EDGE)) +
  scale_size_continuous(range = c(1.2, 3.2), name = "-log10(conservative BH FDR)") +
  labs(title = "Paired PADOG pathway perturbation",
       x = "Mean member-gene expression difference (Edge - Core)", y = NULL, colour = NULL) +
  theme_classic(base_family = FONT_FAMILY, base_size = FONT_PT) +
  theme(plot.title = element_text(face = "bold", size = FONT_PT),
        legend.position = "top", axis.text = element_text(size = FONT_PT),
        axis.title = element_text(size = FONT_PT), legend.text = element_text(size = FONT_PT),
        legend.title = element_text(size = FONT_PT))
save_png(file.path(FIG, "Figure_PADOG_Global", "figure.png"), p_padog,
         width = 5.0, height = 4.8)

# DIRAC: Hallmark-only exploratory analysis. The standard GSReg implementation
# does not include patient blocks; this is therefore not a paired confirmatory
# test. mu is within-phenotype ordering variability (lower = more conserved).
dirac_gs <- gene_sets$Hallmark
dirac_path <- file.path(OUT, "DIRAC_HALLMARK_EXPLORATORY.tsv")
if (file.exists(dirac_path)) {
  dirac_tab <- read_tsv(dirac_path)
} else {
  dirac <- GSReg::GSReg.GeneSets.DIRAC(
    geneexpres = vst,
    pathways = dirac_gs,
    phenotypes = factor(meta$region, levels = c("Core", "Edge")),
    Nperm = 1000,
    alpha = 1,
    minGeneNum = 10
  )
  dirac_tab <- data.frame(
    database = "Hallmark",
    pathway = names(dirac$pvalues),
    Core_variability_index = unname(dirac$mu1[names(dirac$pvalues)]),
    Edge_variability_index = unname(dirac$mu2[names(dirac$pvalues)]),
    delta_Edge_minus_Core = unname(dirac$mu2[names(dirac$pvalues)] - dirac$mu1[names(dirac$pvalues)]),
    zscore = unname(dirac$zscores[names(dirac$pvalues)]),
    pvalue = unname(dirac$pvalues),
    stringsAsFactors = FALSE
  )
  dirac_tab$BH_FDR <- p.adjust(dirac_tab$pvalue, method = "BH")
  dirac_tab$interpretation <- ifelse(dirac_tab$delta_Edge_minus_Core > 0,
                                     "Less rank-conserved at Edge",
                                     "More rank-conserved at Edge")
  write_tsv(dirac_tab, dirac_path)
}

dirac_plot <- head(dirac_tab[order(dirac_tab$BH_FDR, dirac_tab$pvalue), ], 20)
dirac_plot$display <- factor(clean_term(dirac_plot$pathway),
                             levels = rev(clean_term(dirac_plot$pathway)))
write_tsv(dirac_plot, file.path(FIG, "Figure_DIRAC_Hallmark", "source_data.tsv"))
p_dirac <- ggplot(dirac_plot, aes(x = delta_Edge_minus_Core, y = display,
                                  colour = interpretation)) +
  geom_vline(xintercept = 0, colour = "#B3B3B3", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = delta_Edge_minus_Core, yend = display), linewidth = 0.45) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c("Less rank-conserved at Edge" = COL_EDGE,
                                 "More rank-conserved at Edge" = COL_CORE)) +
  labs(title = "DIRAC Hallmark rank conservation (exploratory)",
       x = "DIRAC variability difference (Edge - Core)", y = NULL, colour = NULL) +
  theme_classic(base_family = FONT_FAMILY, base_size = FONT_PT) +
  theme(plot.title = element_text(face = "bold", size = FONT_PT),
        legend.position = "top", axis.text = element_text(size = FONT_PT),
        axis.title = element_text(size = FONT_PT), legend.text = element_text(size = FONT_PT))
save_png(file.path(FIG, "Figure_DIRAC_Hallmark", "figure.png"), p_dirac,
         width = 5.0, height = 4.8)

# Hallmark ssGSEA sample-level scores followed by the paired limma model.
ss_path <- file.path(OUT, "SSGSEA_HALLMARK_SAMPLE_SCORES.tsv")
if (file.exists(ss_path)) {
  ss_df <- read_tsv(ss_path)
  ss_mat <- as.matrix(ss_df[, -1, drop = FALSE])
  storage.mode(ss_mat) <- "double"
  rownames(ss_mat) <- ss_df[[1]]
} else {
  ss_param <- GSVA::ssgseaParam(vst, gene_sets$Hallmark, minSize = 10, maxSize = 500,
                                alpha = 0.25, normalize = TRUE)
  ss_mat <- GSVA::gsva(ss_param, verbose = FALSE)
  write_tsv(data.frame(pathway = rownames(ss_mat), ss_mat, check.names = FALSE), ss_path)
}
ss_fit <- paired_limma(ss_mat, meta)
ss_tab <- ss_fit$stats
colnames(ss_tab)[1] <- "pathway"
write_tsv(ss_tab, file.path(OUT, "SSGSEA_HALLMARK_PAIRED_LIMMA.tsv"))
plot_effect(ss_tab, "pathway", "Hallmark ssGSEA paired differences",
            "ssGSEA score difference (Edge - Core)",
            file.path(FIG, "Figure_ssGSEA_Hallmark"))

# MCP-counter relative cell-population abundance estimates from the frozen VST.
mcp_genes_url <- "https://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"
mcp_genes <- read.delim(mcp_genes_url, check.names = FALSE, stringsAsFactors = FALSE,
                        colClasses = "character")
write_tsv(mcp_genes, file.path(RES, "MCPcounter_genes.txt"))
mcp <- MCPcounter::MCPcounter.estimate(vst, featuresType = "HUGO_symbols", genes = mcp_genes)
mcp <- as.matrix(mcp)
write_tsv(data.frame(population = rownames(mcp), mcp, check.names = FALSE),
          file.path(OUT, "MPCOUNTER_SAMPLE_SCORES.tsv"))
mcp_fit <- paired_limma(mcp, meta)
mcp_tab <- mcp_fit$stats
colnames(mcp_tab)[1] <- "population"
write_tsv(mcp_tab, file.path(OUT, "MPCOUNTER_PAIRED_LIMMA.tsv"))
plot_effect(mcp_tab, "population", "MCP-counter relative abundance",
            "Score difference (Edge - Core)",
            file.path(FIG, "Figure_MCPcounter"), n_show = nrow(mcp_tab))

# ESTIMATE stromal/immune scores. Tumour purity is not reported for RNA-seq.
estimate_input <- data.frame(GeneSymbol = rownames(vst), vst, check.names = FALSE)
est <- tidyestimate::estimate_score(estimate_input, is_affymetrix = FALSE)
est_mat <- t(as.matrix(est[, c("stromal", "immune", "estimate")]))
colnames(est_mat) <- est$sample
est_mat <- est_mat[, meta$sample, drop = FALSE]
rownames(est_mat) <- c("Stromal Score", "Immune Score", "ESTIMATE Score")
write_tsv(data.frame(score = rownames(est_mat), est_mat, check.names = FALSE),
          file.path(OUT, "ESTIMATE_SAMPLE_SCORES.tsv"))
est_fit <- paired_limma(est_mat, meta)
est_tab <- est_fit$stats
colnames(est_tab)[1] <- "score"
write_tsv(est_tab, file.path(OUT, "ESTIMATE_PAIRED_LIMMA.tsv"))
plot_effect(est_tab, "score", "ESTIMATE transcriptomic scores",
            "Score difference (Edge - Core)",
            file.path(FIG, "Figure_ESTIMATE"), n_show = nrow(est_tab))

# Compute log2(TPM+1) from raw counts and featureCounts gene length for IPS.
fc_path <- file.path(ROOT, "results", "04_featurecounts", "featureCounts_gene_id.txt")
fc <- read.delim(fc_path, comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
count_cols <- grep("Aligned.sortedByCoord.out.bam$", colnames(fc), value = TRUE)
sample_ids <- sub("\\.Aligned.sortedByCoord.out.bam$", "", basename(count_cols))
names(count_cols) <- sample_ids
count_cols <- count_cols[meta$sample]
counts <- as.matrix(fc[, unname(count_cols), drop = FALSE])
storage.mode(counts) <- "double"
rownames(counts) <- fc$Geneid
colnames(counts) <- meta$sample
length_kb <- fc$Length / 1000
rate <- counts / length_kb
tpm <- sweep(rate, 2, colSums(rate), "/") * 1e6

anno <- read_tsv(file.path(ROOT, "results", "04_featurecounts", "gene_annotation.tsv"))
map_idx <- match(rownames(tpm), anno$gene_id)
symbols <- anno$gene_name[map_idx]
valid <- !is.na(symbols) & nzchar(symbols)
tpm <- tpm[valid, , drop = FALSE]
symbols <- symbols[valid]
mean_tpm <- rowMeans(tpm)
ord <- order(symbols, -mean_tpm)
tpm <- tpm[ord, , drop = FALSE]
symbols <- symbols[ord]
tpm <- tpm[!duplicated(symbols), , drop = FALSE]
rownames(tpm) <- symbols[!duplicated(symbols)]
log2_tpm <- log2(tpm + 1)
write_tsv(data.frame(gene = rownames(log2_tpm), log2_tpm, check.names = FALSE),
          file.path(OUT, "INPUT_LOG2_TPM_10SAMPLES.tsv"))
write_tsv(data.frame(gene = rownames(log2_tpm), log2_tpm, check.names = FALSE),
          file.path(OUT, "TIDE_INPUT_LOG2_TPM_PREPARED_NOT_RUN.tsv"))

ips_path <- file.path(RES, "IPS_genes.txt")
ips_genes <- read.delim(ips_path, check.names = FALSE, stringsAsFactors = FALSE)
ips_missing <- setdiff(unique(ips_genes$GENE), rownames(log2_tpm))
ips_present <- intersect(unique(ips_genes$GENE), rownames(log2_tpm))
write_tsv(data.frame(
  gene = unique(ips_genes$GENE),
  status = ifelse(unique(ips_genes$GENE) %in% ips_present, "PRESENT", "MISSING")
), file.path(OUT, "IPS_GENE_MAPPING_AUDIT.tsv"))

# Follow the official IPS.R algorithm: sample-wise z-scores across all genes,
# collapse aliases by NAME, apply signed weights, and map AZ to the 0-10 IPS.
ipsmap <- function(x) ifelse(x <= 0, 0, ifelse(x >= 3, 10, round(x * 10 / 3)))
ips_components <- matrix(NA_real_, nrow = 6, ncol = ncol(log2_tpm),
                         dimnames = list(c("MHC", "CP", "EC", "SC", "AZ", "IPS"),
                                         colnames(log2_tpm)))
for (i in seq_len(ncol(log2_tpm))) {
  z_all <- (log2_tpm[, i] - mean(log2_tpm[, i])) / sd(log2_tpm[, i])
  z1 <- z_all[ips_genes$GENE]
  weighted_by_name <- vapply(unique(ips_genes$NAME), function(nm) {
    ii <- which(ips_genes$NAME == nm)
    mean(z1[ii], na.rm = TRUE) * mean(ips_genes$WEIGHT[ii], na.rm = TRUE)
  }, numeric(1))
  mhc <- mean(weighted_by_name[1:10], na.rm = TRUE)
  cp <- mean(weighted_by_name[11:20], na.rm = TRUE)
  ec <- mean(weighted_by_name[21:24], na.rm = TRUE)
  sc <- mean(weighted_by_name[25:26], na.rm = TRUE)
  az <- sum(mhc, cp, ec, sc)
  ips_components[, i] <- c(mhc, cp, ec, sc, az, ipsmap(az))
}
write_tsv(data.frame(component = rownames(ips_components), ips_components, check.names = FALSE),
          file.path(OUT, "IPS_SAMPLE_SCORES.tsv"))
ips_fit <- paired_limma(ips_components, meta)
ips_tab <- ips_fit$stats
colnames(ips_tab)[1] <- "component"
write_tsv(ips_tab, file.path(OUT, "IPS_PAIRED_LIMMA.tsv"))
plot_effect(ips_tab, "component", "Immunophenoscore components",
            "Score difference (Edge - Core)",
            file.path(FIG, "Figure_IPS"), n_show = nrow(ips_tab))

# TIDE execution status. The official local TIDEpy dependency was not added
# because the active project guard rejected new Python dependencies. A valid
# log2(TPM+1) input is prepared above; no proxy score is substituted.
tide_status <- data.frame(
  method = "TIDE",
  status = "NOT_RUN_DEPENDENCY_GUARD",
  input_prepared = "TIDE_INPUT_LOG2_TPM_PREPARED_NOT_RUN.tsv",
  reason = "Official TIDEpy dependency was blocked by active Stop/SCOPE dependency policy; no custom TIDE-like proxy was used.",
  stringsAsFactors = FALSE
)
write_tsv(tide_status, file.path(OUT, "TIDE_EXECUTION_STATUS.tsv"))

# Reproducibility metadata and an evidence-boundary audit.
pkgs <- c("R", "limma", "PADOG", "GSReg", "GSVA", "MCPcounter", "tidyestimate", "ggplot2", "ragg")
versions <- c(as.character(getRversion()), vapply(pkgs[-1], function(p) as.character(packageVersion(p)), character(1)))
metadata <- data.frame(
  item = c("run_date", "contrast", "design", "patients", "samples", paste0("version_", pkgs)),
  value = c(as.character(Sys.Date()), "Edge_minus_Core", "~ patient + region", "5", "10", versions),
  stringsAsFactors = FALSE
)
write_tsv(metadata, file.path(OUT, "RUN_METADATA.tsv"))

audit <- c(
  "# Gene-set and immune sensitivity analysis audit",
  "",
  "- Frozen cohort: 5 paired patients, 10 samples; contrast Edge - Core.",
  "- Frozen input: Checkpoint 10 HGNC-collapsed VST matrix; raw counts and featureCounts gene length used only to derive TPM for IPS/TIDE input.",
  "- CAMERA and fry: paired limma design (~ patient + region), all four frozen GMT databases.",
  "- PADOG: paired=TRUE with patient block. Five pairs provide 31 non-identity label flips. PADOG's zero-exceedance replacement (1/[31*100]) was retained for audit but not treated as an exact P value; conservative exact P=(b+1)/32 and BH correction within database are reported. No PADOG result remains FDR significant under this correction.",
  "- DIRAC: Hallmark-only exploratory standard GSReg implementation; the method does not model patient blocks and is not a paired confirmatory test.",
  "- ssGSEA: Hallmark sample scores followed by paired limma.",
  "- MCP-counter and ESTIMATE: relative transcriptome-derived abundance/context scores, not measured cell fractions.",
  "- IPS: official IPS gene weights and official log2(TPM+1) scoring algorithm; no immunotherapy response endpoint exists in this cohort.",
  "- TIDE: official local execution not run because the active dependency guard rejected TIDEpy installation. A valid input matrix was prepared; no custom proxy was substituted.",
  "- No Venn filtering, machine learning, composite score, PPI rerun, DEG threshold change, or PDE2A-directed post-hoc pathway selection was performed.",
  "- All figures are standalone PNG files at 300 dpi with Arial 7 pt and per-figure source data.",
  "- These analyses are exploratory sensitivity/context analyses and do not alter the frozen PDE2A nomination or external-replication conclusions."
)
writeLines(audit, file.path(OUT, "QC_GENE_SET_IMMUNE_SENSITIVITY_AUDIT.md"), useBytes = TRUE)

# Copy the reproducible plotting/analysis script beside every standalone figure.
script_path <- file.path(ROOT, "scripts", "32_gene_set_immune_sensitivity.R")
for (d in list.dirs(FIG, recursive = FALSE, full.names = TRUE)) {
  file.copy(script_path, file.path(d, "plotting_script.R"), overwrite = TRUE)
}

cat("Completed gene-set and immune sensitivity analyses in:", OUT, "\n")
