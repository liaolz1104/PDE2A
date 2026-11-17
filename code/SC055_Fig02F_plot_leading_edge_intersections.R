# Figure 2E: exact leading-edge overlaps from frozen GSEA, drawn by ComplexUpset.
# Run: Rscript scripts/figure2e_complex_upset.R D:/CodexHome/workspaces/pde2a-bulk
suppressPackageStartupMessages({
  library(ggplot2)
  library(ComplexUpset)
})
root <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else stop('Project root required')
out <- file.path(root, 'results', 'publication_Figure2EF_ComplexUpset_Pareto_20261016')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
read_tsv <- function(path) read.delim(file.path(root, 'results', path), check.names=FALSE,
                                      stringsAsFactors=FALSE, fileEncoding='UTF-8')
gsea <- do.call(rbind, lapply(c('Hallmark','GO_BP','KEGG','Reactome'), function(db) {
  d <- read_tsv(paste0('06_enrichment/GSEA_',db,'.tsv'))
  d$database <- db
  d
}))
modules <- read_tsv('07_candidate_convergence/PATHWAY_TO_FUNCTIONAL_MODULE.tsv')
annotation <- read_tsv('04_featurecounts/gene_annotation.tsv')
front <- read_tsv('11_candidate_prioritization/CANDIDATE_PARETO_FRONTS.tsv')
split_genes <- function(s) unique(trimws(strsplit(s,';',fixed=TRUE)[[1]]))
has_pde2a <- vapply(gsea$leadingEdge, function(s) !is.na(s) && 'PDE2A' %in% split_genes(s), logical(1))
paths <- gsea[!is.na(gsea$padj) & gsea$padj < .05 & gsea$NES < 0 & has_pde2a,]
stopifnot(nrow(paths)==11, !anyDuplicated(paste(paths$database,paths$pathway)))
key <- paste(paths$database, paths$pathway,sep='\r')
module_key <- paste(modules$database, modules$pathway,sep='\r')
paths$functional_module <- modules$functional_module[match(key,module_key)]
stopifnot(!anyNA(paths$functional_module))
module_order <- c('Cyclic nucleotide signaling','NO / cGMP signaling',
                  'GPCR / adenylate cyclase signaling','Growth-factor signaling',
                  'Hemostasis / platelet','Metabolism / endocrine')
stopifnot(setequal(paths$functional_module,module_order))
paths <- paths[order(match(paths$functional_module,module_order),paths$pathway),]
rownames(paths) <- NULL
paths$set_id <- sprintf('S%02d',seq_len(nrow(paths)))
sets <- lapply(paths$leadingEdge,split_genes)
names(sets) <- paths$set_id
all_genes <- sort(unique(unlist(sets,use.names=FALSE)))
stopifnot(length(all_genes)==215, all(all_genes %in% annotation$gene_name))
binary <- as.data.frame(sapply(sets,function(s) all_genes %in% s),check.names=FALSE)
binary <- cbind(gene=all_genes,binary)
stopifnot(!anyDuplicated(binary$gene), all(colSums(binary[paths$set_id])==lengths(sets)))
signature <- apply(binary[paths$set_id],1,function(row) paste(paths$set_id[as.logical(row)],collapse=';'))
stopifnot(!any(signature==''))
binary$exact_intersection_id <- signature
summary <- aggregate(gene ~ exact_intersection_id,binary,length)
names(summary)[2] <- 'exact_size'
summary$degree <- lengths(strsplit(summary$exact_intersection_id,';',fixed=TRUE))
summary <- summary[order(-summary$exact_size,summary$exact_intersection_id),]
summary$intersection_rank <- seq_len(nrow(summary))
summary$shown_top15 <- summary$intersection_rank<=15
stopifnot(nrow(summary)==39, sum(summary$exact_size)==length(all_genes))
shown <- summary[summary$shown_top15,]
covered <- sum(shown$exact_size)
stopifnot(covered==186)
summary$inclusive_size <- vapply(summary$exact_intersection_id,function(sig) {
  selected <- strsplit(sig,';',fixed=TRUE)[[1]]
  sum(rowSums(binary[selected])==length(selected))
},integer(1))
stopifnot(all(summary$inclusive_size>=summary$exact_size))
binary$exact_size <- summary$exact_size[match(binary$exact_intersection_id,summary$exact_intersection_id)]
binary$shown_top15 <- binary$exact_intersection_id %in% shown$exact_intersection_id
paths$set_size <- lengths(sets)
paths$display_order <- seq_len(nrow(paths))
write.table(paths,file.path(out,'Figure2E_frozen_pathways_complete_leading_edges.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
write.table(binary,file.path(out,'Figure2E_all_215_gene_membership.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
write.table(summary,file.path(out,'Figure2E_all_39_intersections.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
write.table(paths[,c('set_id','database','pathway','functional_module','NES','padj','set_size')],
            file.path(out,'Figure2E_set_annotations.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')

# User-selected top 15 is defined by complete exclusive-intersection sizes.
chosen <- lapply(shown$exact_intersection_id,function(s) strsplit(s,';',fixed=TRUE)[[1]])
names(chosen) <- shown$exact_intersection_id
plot_data <- binary[,paths$set_id,drop=FALSE]
stopifnot(all(vapply(chosen,function(s) {
  sum(rowSums(plot_data[,s,drop=FALSE])==length(s) &
        rowSums(plot_data[,setdiff(paths$set_id,s),drop=FALSE])==0)
},integer(1))==shown$exact_size))

# The ComplexUpset plot draws all standard UpSet elements; labels carry the
# three metadata fields beside each formal pathway name.
labels <- setNames(vapply(seq_len(nrow(paths)), function(i) {
  p <- paths[i,]
  name <- paste(strwrap(p$pathway,width=57),collapse='\n')
  sprintf('%s\n%s  |  NES %.2f  |  FDR %s',name,p$functional_module,
          p$NES,if (p$padj<.001) '<0.001' else sprintf('%.3f',p$padj))
},character(1)),paths$set_id)
labeller <- function(x) unname(labels[as.character(x)])
theme_set(theme_classic(base_family='Arial',base_size=8))
p <- upset(plot_data,rev(paths$set_id),
           intersections=chosen,sort_intersections=FALSE,sort_sets=FALSE,
           mode='exclusive_intersection',encode_sets=FALSE,keep_empty_groups=TRUE,
           name='Core-enriched gene sets',labeller=labeller,
           width_ratio=.18,height_ratio=1.05,stripes='white',
           matrix=intersection_matrix(geom=geom_point(size=1.35),
                                      segment=geom_segment(linewidth=.5),
                                      outline_color=list(active='black',inactive='#C4CDD3')),
           base_annotations=list('Number of genes'=intersection_size(
             counts=TRUE,mode='exclusive_intersection') +
             scale_y_continuous(expand=expansion(mult=c(0,.08))) +
             ylab('Number of genes')),
           set_sizes=upset_set_size(filter_intersections=FALSE)+
             ylab('Leading-edge genes (n)')) &
             theme(text=element_text(family='Arial',size=8),
                   plot.background=element_rect(fill='white',colour=NA),
                   panel.background=element_rect(fill='white',colour=NA),
                   axis.line=element_line(linewidth=.5,color='black'),
                   axis.text=element_text(color='black',size=7))

caption <- sprintf('Exact intersections shown: 15/39; union genes covered: %d/%d (%.1f%%). Set sizes use all genes. Inclusive intersections are larger where other sets overlap.',
                   covered,length(all_genes),100*covered/length(all_genes))
p <- p + patchwork::plot_annotation(title='E  Leading-edge gene-set intersections',
                                    subtitle='Core-enriched PDE2A-containing gene sets; module labels are descriptive',
                                    caption=caption,
                                    theme=theme(plot.title=element_text(family='Arial',size=9,face='bold'),
                                                plot.subtitle=element_text(family='Arial',size=8),
                                                plot.caption=element_text(family='Arial',size=7,hjust=0)))
svglite::svglite(file.path(out,'Figure2E_ComplexUpset.svg'),width=8.2,height=7.4,system_fonts=list(sans='Arial'))
print(p)
dev.off()
cairo_pdf(file.path(out,'Figure2E_ComplexUpset.pdf'),width=8.2,height=7.4,family='Arial')
print(p)
dev.off()
png_tmp <- 'D:/CodexHome/figure2e_complex_upset_tmp.png'
ragg::agg_png(png_tmp,width=8.2,height=7.4,units='in',res=300,background='white')
print(p)
dev.off()
stopifnot(file.exists(png_tmp), file.copy(png_tmp,file.path(out,'Figure2E_ComplexUpset.png'),overwrite=TRUE))
unlink(png_tmp)
cat(sprintf('union=%d intersections=%d top15=%d coverage=%.4f\n',length(all_genes),nrow(summary),covered,covered/length(all_genes)))
