rm(list=ls())
library(ggplot2)
library(pheatmap)

source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
ddir = '/home/whou10/scratch4/whou10/encode4/data/pseudobulk/snrna/pb/'
list.files(ddir)
for (species in c('Homo_sapiens', 'Mus_musculus')){
  pd = readRDS(paste0(ddir, species, '_pb.rds'))
  meta = readRDS(paste0(ddir, species, "_meta.rds"))
  str(meta)
  ## retain pseudobulks with > 50 cells
  num_cell <- as.numeric(sapply(colnames(pd), function(i) strsplit(i, ':')[[1]][3]))
  summary(num_cell)
  pd <- pd[, num_cell > 50]
  pd <- pd[rowMeans(pd > median(pd[pd > 0])) > 0.1, ] 
  summary(rowsds(pd))
  pd <- pd[rowsds(pd) > quantile(rowsds(pd), 0.2), , drop = F]
  
  a=pheatmap(pd, show_rownames = F, show_colnames = F, scale = 'row', breaks = seq(-1.5, 1.5, 3/100))
  
  coln = a$tree_col$labels[a$tree_col$order]
  dataset = sapply(coln, function(i) strsplit(i, ':')[[1]][1])
  tissue = meta[match(dataset, meta[,'dataset']), 'generaltissue']
  
  annotation_col = data.frame(celltype = sapply(coln, function(i) strsplit(i, ':')[[1]][2]), 
                              tissue = tissue,
                              stringsAsFactors = F)
  rownames(annotation_col) = coln
  
  coltb = read.table('/home/whou10/scratch4/whou10/encode4/circlize/color_scheme.csv', sep = ',')
  tissuecol = coltb[,2]
  names(tissuecol) = coltb[,1]
  tissuecol = tissuecol[names(tissuecol) %in% unique(tissue)]
  str(tissuecol)
  annotation_colors = list(tissue = tissuecol)
  
  pdf(paste0('/home/whou10/scratch4/whou10/encode4/pseudobulk_umap/snrna/pb/plot/', ifelse(species=='Mus_musculus', 'mouse', 'human'), '/gene_by_pseudobulk_log2tpm_scaleByRow3.pdf'), width = 28, height = 25)
  pheatmap(pd, show_rownames = F, show_colnames = F, scale = 'row', breaks = seq(-1.5, 1.5, 3/100),
           annotation_col = annotation_col,
           annotation_colors = annotation_colors)
  dev.off()
}
