library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
alltissue = list.files('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_no_downsample/res/')

numcell <- sapply(alltissue, function(tissue){
  ddir <- paste0('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_no_downsample/res/', tissue, '/figures/')
  tb = read.csv(paste0(ddir, 'cell_participation.csv'), row.names = 1)
  nrow(tb)
})

numtopic <- sapply(alltissue, function(tissue){
  ddir <- paste0('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_no_downsample/res/', tissue, '/figures/')
  tb = read.csv(paste0(ddir, 'cell_participation.csv'), row.names = 1)
  ncol(tb)
})

for (tissue in names(sort(numcell))){
  print(tissue)
  pdir <- '/home/whou10/scratch4/whou10/encode4/topyfic/compareversion/dendrogram/original/plot/confusionMatrix_clusterByTopics_celltype/hclu/'
  dir.create(pdir, showWarnings = F)
  rdir <- '/home/whou10/scratch4/whou10/encode4/topyfic/compareversion/dendrogram/original/res/'
  clu = readRDS(paste0(rdir, tissue, '_hclu.rds'))
  
  meta <- readRDS('/home/whou10/scratch4/whou10/encode4/data/celltype/proc/final.rds')
  ct = meta$celltype
  names(ct) = paste0(meta$rna_dataset, ':', meta$rna_library, ':', meta$rna_barcode)
  ct = ct[names(clu)]
  
  all_levels = c(paste0('cluster',unique(clu)), unique(ct))
  
  ct_factor = factor(ct, levels = all_levels)
  clu_factor = factor(paste0('cluster', clu), levels = all_levels)
  
  ## generate confusion matrix
  tb <- table(celltype = ct_factor, clustered_by_topics = clu_factor)
  
  ## scaled by num of cells in each celltype to rule out the case that a large celltype drive the color scheme
  tb <- tb/rowSums(tb) 
  flag <- is.na(tb)
  tb2 <- tb[rowMeans(flag)!=1, ]
  flag <- is.na(tb2)
  tb3 <- tb2[,colMeans(flag)!=1]
  tb3 <- tb3[,colMeans(tb3)>0]
  pdf(paste0(pdir, tissue,'.pdf'))
  print(pheatmap(tb3, cluster_rows = T, cluster_cols = T))
  dev.off()
}

