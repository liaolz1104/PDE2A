#!/usr/bin/env Rscript
# Reuse the reviewed paired-expression plot, removing its heading and placing
# the frozen whole-transcriptome DESeq2 BH-adjusted P inside the plotting area.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
root <- sub("/+$", "", args[[1]])
Sys.setenv(PDE2A_SKIP_BUNDLE_EXPORT = "1")
# Evaluate the plot-building portion only; avoid exporting unrelated PCA files.
input_script <- file.path(root, "scripts/42_polish_paired_figures_arial.R")
stopifnot(file.exists(input_script))
code <- readLines(input_script, warn = FALSE)
export_line <- which(grepl('^script <- file.path', code))
stopifnot(length(export_line) == 1L)
eval(parse(text = code[seq_len(export_line - 1L)]))
Sys.setenv(PDE2A_SKIP_BUNDLE_EXPORT = "0")
stopifnot(length(unique(expr$DESeq2_padj)) == 1L)
q <- q + labs(title = NULL, subtitle = NULL) +
  annotate("text", x = 1.68, y = 9.32,
    label = sprintf("FDR = %.3f", expr$DESeq2_padj[1]),
    family = "Arial", size = FONT_MM, colour = "#222222") +
  theme(plot.title = element_blank(), plot.subtitle = element_blank())
name <- "Figure_03_PDE2A_Paired_Arial_InPlotFDR"
save_bundle(out, name, q, 155, 140, list(data = expr),
            file.path(root, "scripts/43_pde2a_paired_inplot_fdr.R"))
file.copy(input_script, file.path(out, name, "plot_building_dependency.R"), overwrite = TRUE)
write_tsv(data.frame(statistic = "Whole-transcriptome paired DESeq2 BH-adjusted P",
  exact_value = expr$DESeq2_padj[1], display = sprintf("FDR = %.3f", expr$DESeq2_padj[1])),
  file.path(out, name, "STATISTIC_ANNOTATION.tsv"))
cat("Expression figure completed with in-plot FDR and no heading.\n")
