#!/usr/bin/env Rscript

# Redraw the DEG heatmap with region-grouped columns and separated directions.
# Rows include adjusted-P ranks 1 through PDE2A. Columns are grouped as all Core
# followed by all Edge. Edge-down and Edge-up rows are clustered separately.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 37_redraw_clustered_heatmap_through_pde2a.R <project_root>")
root <- sub("/+$", "", args[[1]])

helper_path <- file.path(root, "scripts", "30_manuscript_redraw_helpers.R")
deseq_path <- file.path(root, "results", "05_deseq2", "DESeq2_Edge_vs_Core_all_genes.tsv")
vst_path <- file.path(root, "results", "04_featurecounts", "VST_all12_QC.rds")
index_path <- file.path(
  root, "results", "07_figures", "manuscript_redraw_v3_single_panel",
  "STANDALONE_FIGURE_INDEX.tsv"
)
required <- c(helper_path, deseq_path, vst_path, index_path)
if (any(!file.exists(required))) {
  stop("Missing required input: ", paste(required[!file.exists(required)], collapse = "; "))
}

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(pheatmap)
  source(helper_path)
})

# Register the actual Arial TrueType files rather than accepting font fallback.
font_dir <- if (.Platform$OS.type == "windows") "C:/Windows/Fonts" else
  file.path(root, "assets", "fonts", "Arial")
font_files <- file.path(font_dir, c("arial.ttf", "arialbd.ttf", "ariali.ttf", "arialbi.ttf"))
if (!all(file.exists(font_files))) stop("Missing required Arial TrueType files")
if (!"Arial" %in% unique(systemfonts::system_fonts()$family)) {
  systemfonts::register_font("Arial", plain = font_files[1], bold = font_files[2],
                            italic = font_files[3], bolditalic = font_files[4])
}

deseq <- read_tsv(deseq_path)
deseq$padj_num <- suppressWarnings(as.numeric(deseq$padj))
deseq$pvalue_num <- suppressWarnings(as.numeric(deseq$pvalue))
deseq$gene_id_stripped_local <- sub("\\..*$", "", deseq$gene_id)
rankable <- deseq[is.finite(deseq$padj_num), , drop = FALSE]
rankable <- rankable[order(rankable$padj_num, rankable$pvalue_num, rankable$gene_id), , drop = FALSE]
rankable$adjusted_p_rank <- seq_len(nrow(rankable))
pde_rows <- rankable[rankable$gene_name == "PDE2A", , drop = FALSE]
if (nrow(pde_rows) != 1L) stop("PDE2A must map to exactly one ranked DESeq2 row; found ", nrow(pde_rows))
pde_rank <- pde_rows$adjusted_p_rank[[1]]
ranked_genes <- rankable[seq_len(pde_rank), , drop = FALSE]

vsd <- readRDS(vst_path)
vst_mat <- SummarizedExperiment::assay(vsd)
if (is.null(rownames(vst_mat)) || anyDuplicated(rownames(vst_mat))) {
  stop("VST gene identifiers are missing or duplicated")
}

patients <- c("P08", "P09", "P34", "P36", "P38")
core_samples <- paste0(patients, "_C")
edge_samples <- paste0(patients, "_TI")
main_samples <- c(core_samples, edge_samples)
if (!all(main_samples %in% colnames(vst_mat))) stop("VST matrix lacks frozen main samples")

vst_id <- rownames(vst_mat)
match_exact <- match(ranked_genes$gene_id, vst_id)
match_stripped <- match(ranked_genes$gene_id_stripped_local, sub("\\..*$", "", vst_id))
vst_row <- ifelse(is.na(match_exact), match_stripped, match_exact)
if (any(is.na(vst_row))) {
  stop("VST rows missing for ranked genes: ",
       paste(ranked_genes$gene_name[is.na(vst_row)], collapse = ", "))
}

ranked_genes$direction <- ifelse(
  ranked_genes$log2FoldChange_shrunk < 0, "Edge-down",
  ifelse(ranked_genes$log2FoldChange_shrunk > 0, "Edge-up", NA_character_)
)
if (anyNA(ranked_genes$direction)) {
  stop("Every displayed gene must have a non-zero shrunken Edge-Core direction")
}
down_idx <- which(ranked_genes$direction == "Edge-down")
up_idx <- which(ranked_genes$direction == "Edge-up")
if (!length(down_idx) || !length(up_idx)) {
  stop("Both Edge-down and Edge-up blocks are required")
}

expr_ranked <- vst_mat[vst_row, main_samples, drop = FALSE]
rownames(expr_ranked) <- ranked_genes$gene_id
row_z_ranked <- t(scale(t(expr_ranked)))
row_z_ranked[!is.finite(row_z_ranked)] <- 0

cluster_row_order <- function(x) {
  if (nrow(x) <= 1L) return(seq_len(nrow(x)))
  row_cor <- cor(t(x), use = "pairwise.complete.obs")
  row_cor[!is.finite(row_cor)] <- 0
  row_cor <- pmin(pmax(row_cor, -1), 1)
  hclust(as.dist(1 - row_cor), method = "complete")
}

up_tree <- cluster_row_order(row_z_ranked[up_idx, , drop = FALSE])
down_tree <- cluster_row_order(row_z_ranked[down_idx, , drop = FALSE])
up_order <- up_idx[up_tree$order]
down_order <- down_idx[down_tree$order]
row_order <- c(up_order, down_order)
ranked_genes <- ranked_genes[row_order, , drop = FALSE]
expr <- expr_ranked[row_order, , drop = FALSE]
row_z <- row_z_ranked[row_order, , drop = FALSE]

gene_label <- ranked_genes$gene_name
gene_label[is.na(gene_label) | gene_label == ""] <-
  ranked_genes$gene_id_stripped_local[is.na(gene_label) | gene_label == ""]
if (anyDuplicated(gene_label)) stop("Displayed gene labels are duplicated")
rownames(row_z) <- gene_label

patient <- sub("_(C|TI)$", "", main_samples)
region <- ifelse(grepl("_C$", main_samples), "Core", "Edge")
column_labels <- paste(patient, region, sep = "\n")
annotation_col <- data.frame(
  Patient = factor(patient, levels = patients),
  Region = factor(region, levels = c("Core", "Edge")),
  row.names = main_samples,
  stringsAsFactors = FALSE
)
annotation_colors <- list(
  Patient = c(P08 = "#4E79A7", P09 = "#59A14F", P34 = "#F28E2B",
              P36 = "#B07AA1", P38 = "#E15759"),
  Region = c(Core = COL_CORE, Edge = COL_EDGE)
)

mat <- row_z
colnames(mat) <- main_samples
max_abs <- max(abs(mat), na.rm = TRUE)
heat_colours <- colorRampPalette(c("#3B5791", "#F7F7F7", "#B40426"))(101)
heat_breaks <- seq(-max_abs, max_abs, length.out = length(heat_colours) + 1L)

plot_obj <- pheatmap(
  mat,
  color = heat_colours,
  breaks = heat_breaks,
  scale = "none",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_col = length(core_samples),
  gaps_row = length(up_order),
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  annotation_names_col = TRUE,
  annotation_names_row = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  labels_col = column_labels,
  fontfamily = "Arial",
  fontsize = 7.0,
  fontsize_row = 7.0,
  fontsize_col = 7.0,
  angle_col = 45,
  border_color = NA,
  treeheight_row = 52,
  main = "Core-Edge differential expression",
  legend = TRUE,
  silent = TRUE
)

# Draw the two independent hclust trees in the matrix's own row viewport.
# Match pheatmap's 4-bigpt block gap exactly; never join unrelated block roots.
tree_segments <- function(tree, offset, block) {
  n <- length(tree$order)
  leaf_y <- offset + match(seq_len(n), tree$order) - 0.5
  node_y <- numeric(n - 1L)
  segments <- list()
  for (i in seq_len(n - 1L)) {
    children <- tree$merge[i, ]
    cy <- vapply(children, function(k) if (k < 0) leaf_y[-k] else node_y[k], numeric(1))
    ch <- vapply(children, function(k) if (k < 0) 0 else tree$height[k], numeric(1))
    node_y[i] <- mean(cy)
    segments[[i]] <- data.frame(block = block,
      x0 = c(ch, tree$height[i]), x1 = rep(tree$height[i], 3),
      y0 = c(cy, cy[1]), y1 = c(cy, cy[2]))
  }
  do.call(rbind, segments)
}
branches <- rbind(tree_segments(up_tree, 0, "Edge-up"),
                  tree_segments(down_tree, length(up_order), "Edge-down"))
max_height <- max(c(up_tree$height, down_tree$height))
tree_y <- function(y, block) {
  unit(1, "npc") - (unit(1, "npc") - unit(4, "bigpts")) * y / nrow(mat) -
    unit(ifelse(block == "Edge-down", 4, 0), "bigpts")
}
tree_grob <- segmentsGrob(
  x0 = unit(0.96 - 0.9 * branches$x0 / max_height, "npc"),
  x1 = unit(0.96 - 0.9 * branches$x1 / max_height, "npc"),
  y0 = tree_y(branches$y0, branches$block),
  y1 = tree_y(branches$y1, branches$block),
  gp = gpar(col = "#444444", lwd = 0.7))
matrix_cell <- plot_obj$gtable$layout[plot_obj$gtable$layout$name == "matrix", ]
added_tree_width <- unit(18, "mm") - plot_obj$gtable$widths[1]
plot_obj$gtable$widths[matrix_cell$l] <- plot_obj$gtable$widths[matrix_cell$l] - added_tree_width
plot_obj$gtable$widths[1] <- unit(18, "mm")
plot_obj$gtable <- gtable::gtable_add_grob(plot_obj$gtable, tree_grob,
  t = matrix_cell$t, l = 1, name = "independent_direction_dendrograms", clip = "on")

# Apply Arial to every nested text grob, including annotation legends.
set_arial <- function(g) {
  if (inherits(g, "text")) {
    g$gp$fontfamily <- "Arial"
    g$gp$fontsize <- 7
  }
  if (length(g$grobs)) g$grobs <- lapply(g$grobs, set_arial)
  if (length(g$children)) {
    for (k in seq_along(g$children)) g$children[[k]] <- set_arial(g$children[[k]])
  }
  g
}
plot_obj$gtable <- set_arial(plot_obj$gtable)

heat_source <- do.call(rbind, lapply(seq_along(gene_label), function(i) {
  data.frame(
    adjusted_p_rank = ranked_genes$adjusted_p_rank[[i]],
    cluster_order = i,
    direction_block = ranked_genes$direction[[i]],
    direction_block_order = ifelse(
      ranked_genes$direction[[i]] == "Edge-up",
      i,
      i - length(up_order)
    ),
    gene_id = ranked_genes$gene_id[[i]],
    gene_id_stripped = ranked_genes$gene_id_stripped_local[[i]],
    gene_label = gene_label[[i]],
    target_label = ifelse(gene_label[[i]] == "PDE2A", "PDE2A", "Other"),
    sample_id = main_samples,
    sample_label = paste(patient, region, sep = "_"),
    patient_id = patient,
    region = region,
    row_z = as.numeric(row_z[i, ]),
    stringsAsFactors = FALSE
  )
}))

out_root <- file.path(root, "results", "07_figures", "manuscript_redraw_v4_Arial7pt")
figure_name <- "Figure_06_CoreEdge_SplitDendrogramHeatmap"
plotting_script <- file.path(root, "scripts", "40_redraw_heatmap_split_dendrograms.R")
save_bundle(
  output_root = out_root,
  figure_name = figure_name,
  plot = plot_obj,
  width_mm = 180,
  height_mm = 250,
  sources = list(data = heat_source),
  plotting_script = plotting_script
)
write_tsv(branches, file.path(out_root, figure_name, "dendrogram_segments.tsv"))
write_tsv(data.frame(family = "Arial", style = c("Regular", "Bold", "Italic", "Bold Italic"),
                    format = "TrueType", file = font_files,
                    md5 = unname(tools::md5sum(font_files))),
          file.path(out_root, figure_name, "FONT_AUDIT.tsv"))
saveRDS(list(Edge_up = up_tree, Edge_down = down_tree),
        file.path(out_root, figure_name, "row_clustering_trees.rds"))

index <- read_tsv(index_path)
old_names <- c("Figure_06_TopDEG_to_PDE2A_CoreEdgeGroupedHeatmap", "Figure_06_TopDEG_to_PDE2A_ClusteredHeatmap",
               "Figure_06_TopDEG_to_PDE2A_Heatmap", "Figure_06_Top30DEG_Heatmap")
old_present <- old_names[old_names %in% index$figure]
if (length(old_present) == 1L) {
  index$figure[index$figure == old_present] <- figure_name
  index$replaces[index$figure == figure_name] <- "Figure_DEG_03_TopDEG_heatmap"
  index$interpretation_boundary[index$figure == figure_name] <-
  sprintf("Adjusted-P ranks 1 through PDE2A rank %d; all Core columns followed by all Edge columns; Edge-down and Edge-up rows separated and clustered within block", pde_rank)
}
if (sum(index$figure == figure_name) != 1L) stop("Figure index update failed")
index <- index[order(as.numeric(index$figure_order)), , drop = FALSE]
write_tsv(index, index_path)

for (legacy_name in old_present) {
  legacy_dir <- file.path(out_root, legacy_name)
  if (dir.exists(legacy_dir) && legacy_name != figure_name) {
    writeLines(c(
      "PRESERVED HEATMAP FIGURE",
      "",
      paste0("Use ../", figure_name, " for the current clustered redraw."),
      "The prior source data and PNG are retained for provenance."
    ), file.path(legacy_dir, "SUPERSEDED_BY_CLUSTERED_HEATMAP.txt"), useBytes = TRUE)
  }
}

note_path <- file.path(out_root, "FIGURE_06_CLUSTERED_HEATMAP_REVISION_NOTE.md")
writeLines(c(
  "# Figure 06 revision",
  "",
  paste0("- Current folder: `", figure_name, "`"),
  paste0("- Rows: adjusted-P rank 1 through PDE2A (rank ", pde_rank, ")."),
  "- Columns: all five Core samples are on the left, followed by all five Edge samples on the right.",
  "- Edge-down and Edge-up rows are shown as separate blocks with a row gap; each block is hierarchically clustered using 1 - Pearson correlation of row-z-scored frozen VST values.",
  "- Patient and Region annotations are shown; the extra Target legend was removed to prevent right-side stacking; PDE2A remains a visible row label.",
  "- No differential-expression calculation or threshold was changed.",
  "- Prior heatmap versions are preserved for provenance."
), note_path, useBytes = TRUE)

cat("Clustered heatmap redraw completed: ", figure_name,
    "; rows through PDE2A = ", pde_rank, "\n", sep = "")
