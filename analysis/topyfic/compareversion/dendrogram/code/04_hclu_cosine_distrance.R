rm(list=ls())
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggdendro)
library(plotly)
library(pheatmap)
source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
setwd('/home/whou10/scratch4/whou10/encode4/topyfic/')

ddir1 = 'tissue_hvg/res/'
rdir = 'compareversion/dendrogram/hvg/res/'
pdir = 'compareversion/dendrogram/hvg/plot/dendrogram_celltype/'

alltissue = list.files(ddir1)

## calculate numcell for the tissue, so that we can run from smallest tissue, 
## avoid the job being killed, at the beginning
numcell <- sapply(alltissue, function(tissue){
  ddir <- paste0(ddir1, tissue, '/figures/')
  if (file.exists(paste0(ddir, 'cell_participation.csv'))){
    tb = read.csv(paste0(ddir, 'cell_participation.csv'), row.names = 1)
    nrow(tb)
  } else {
    NA
  }
})
numcell <- numcell[!is.na(numcell)]

tissue = "Homo_sapiens:adrenal_gland"
res <- lapply(names(sort(numcell))[21:length(numcell)], function(tissue){
  print(tissue)
  ## load topics
  ddir <- paste0(ddir1, tissue, '/figures/')
  tb = as.matrix(read.csv(paste0(ddir, 'cell_participation.csv'), row.names = 1))
  if (nrow(tb) > 1e4) {
    set.seed(123)
    tb = tb[sample(1:nrow(tb), 1e4, replace = F), ]
  }
  
  ## ct
  meta <- readRDS('/home/whou10/scratch4/whou10/encode4/data/celltype/proc/final.rds')
  meta <- meta[!is.na(meta$rna_dataset),]
  rownames(meta) <- paste0(meta$rna_dataset, ':', meta$rna_library, ':', meta$rna_barcode)
  ct <- meta[rownames(tb),'celltype']
  names(ct) = rownames(tb)
  
  ## topic mean
  rs <- rowsum(tb,ct)
  rs <- rs/as.vector(table(ct)[rownames(rs)])
  rs <- t(rs[,order(apply(rs, 2, which.max))])
  pdf(paste0(pdir, 'topic_mean/', tissue, '.pdf'))
  pheatmap::pheatmap(rs, cluster_cols = F, cluster_rows = F, scale = 'none')
  dev.off()
  
  ## calculate distance 
  tb2 <- tb[rep(1:nrow(tb),nrow(tb)),]
  tb3 <- tb[rep(1:nrow(tb),each=nrow(tb)),]
  dis <- matrix(rowSums(tb2 * tb3) / (sqrt(rowSums(tb2^2)) * sqrt(rowSums(tb3^2))),nrow=nrow(tb))
  dimnames(dis) <- list(rownames(tb), rownames(tb))
  saveRDS(dis, paste0(pdir, 'cosine_dist_of_cell_contribution/', tissue, '.rds'))
  
  ## clustering
  # hclu1 <- fastcluster::hclust(as.dist(1-dis))
  # clu1 <- cutree(hclu1, k = length(unique(ct)))
  hclu <- flashClust::flashClust(as.dist(1-dis))
  clu <- cutree(hclu, k = length(unique(ct)))
  saveRDS(clu, paste0(rdir, 'hclu_cosine/', tissue, '.rds' ))
  
  ## confusion
  mat = table(clu[names(ct)], ct)
  mat = t(t(mat)/rowSums(t(mat)))
  mat = mat[, order(apply(mat, 2, which.max))]
  pdf(paste0(pdir, 'hclu_cosine/', tissue, '.pdf'))
  pheatmap(mat, cluster_cols = F, cluster_rows = F, scale = 'none')
  dev.off()
  return(0)
})

  