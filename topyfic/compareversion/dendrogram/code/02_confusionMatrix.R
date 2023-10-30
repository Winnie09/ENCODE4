library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(caret)
source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
tissuefd <- commandArgs(trailingOnly = T)[1]
pdir <- commandArgs(trailingOnly = T)[2]
rdir <- commandArgs(trailingOnly = T)[3]
alltissue = list.files(tissuefd)

for (tissue in alltissue){
  print(tissue)
  dir.create(pdir, showWarnings = F)
  clu = readRDS(paste0(rdir, tissue, '_hclu.rds'))
  
  meta <- readRDS('/home/whou10/scratch4/whou10/encode4/data/celltype/proc/final.rds')
  ct = meta$celltype
  names(ct) = paste0(meta$rna_dataset, ':', meta$rna_library, ':', meta$rna_barcode)
  ct = ct[names(clu)]
  
  all_levels = c(paste0('cluster',unique(clu)), unique(ct))
  
  ct_factor = factor(ct, levels = all_levels)
  clu_factor = factor(paste0('cluster', clu), levels = all_levels)
  names(clu_factor) = names(clu)
  
  m <- caret::confusionMatrix(data = clu_factor, reference = ct_factor,
                              dnn = c('clustered_by_topics', 'celltype'))  
  tb <- t(m$table)
  tb <- tb/rowSums(tb)
  flag <- is.na(tb)
  tb2 <- tb[rowMeans(flag)!=1, ]
  flag <- is.na(tb2)
  tb3 <- tb2[,colMeans(flag)!=1]
  tb3 <- tb3[,colMeans(tb3)>0]
  pdf(paste0(pdir, tissue,'.pdf'))
  print(pheatmap(tb3, cluster_rows = F, cluster_cols = F))
  dev.off()
}


