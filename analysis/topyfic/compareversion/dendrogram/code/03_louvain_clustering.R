library(ggplot2)
library(dplyr)
library(RColorBrewer)
#library(ggdendro)
#library(plotly)
library(pheatmap)
suppressMessages(library(scran))
suppressMessages(library(igraph))
source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
setwd('/home/whou10/scratch4/whou10/encode4/topyfic/')


# ddir1 <- 'tissue_downsample/res/'
# rdir <- 'compareversion/dendrogram/downsample/res/'
# pdir <- 'compareversion/dendrogram/downsample/plot/confusionMatrix_clusterByTopics_celltype/'

ddir1 = 'tissue_hvg/res/'
rdir = 'compareversion/dendrogram/hvg/res/'
pdir = 'compareversion/dendrogram/hvg/plot/confusionMatrix_clusterByTopics_celltype/'

ddir1 <- commandArgs(trailingOnly = T)[1]
rdir = commandArgs(trailingOnly = T)[2]
pdir = commandArgs(trailingOnly = T)[3]
dir.create(pdir, showWarnings = F, recursive = T)

alltissue = list.files(ddir1)

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

for (tissue in names(sort(numcell))){
  tryCatch({
    # The code that might throw an error
    ## ==============================
    ## cluster cells based on topics
    ## ===============================
    print(tissue)
    ddir <- paste0(ddir1, tissue, '/figures/')
    tb3 = read.csv(paste0(ddir, 'cell_participation.csv'), row.names = 1)
    meta <- readRDS('/home/whou10/scratch4/whou10/encode4/data/celltype/proc/final.rds')
    rn = paste0(meta$rna_dataset, ':', meta$rna_library, ':', meta$rna_barcode)
    meta$label <- rn
    meta <- meta[meta$label %in% rownames(tb3), ]
    
    # Define the colors
    v = unique(meta$celltype)
    color <- colorRampPalette(brewer.pal(n = 8, name = "Set1"))(length(v))
    names(color) = v
    
    ## Louvain clustering
    k = length(unique(meta$celltype))
    graph = buildSNNGraph(tb3, transposed=T,k=k,d=NA)
    res = cluster_louvain(graph)$membership
    if (max(res) <= k){
      hclu <- flashClust::flashClust(dist(tb3))
      clu <- cutree(hclu,k)
    } else {
      cc <- aggregate(tb3,list(res),mean)
      cc <- as.matrix(cc[,-1])
      hclu <- hclust(dist(cc))
      clu <- cutree(hclu,k)
      clu <- clu[res]      
    }  
    names(clu) = rownames(tb3)
    saveRDS(clu, paste0(rdir, tissue, '_louvainclu.rds'))
    
    ## ================
    ## Confusion Matrix
    ## ================
    ## prepare cell types
    meta <- readRDS('/home/whou10/scratch4/whou10/encode4/data/celltype/proc/final.rds')
    ct = meta$celltype
    names(ct) = paste0(meta$rna_dataset, ':', meta$rna_library, ':', meta$rna_barcode)
    ct = ct[names(clu)]
    all_levels = c(paste0('cluster',unique(clu)), unique(ct))
    
    ct_factor = factor(ct, levels = all_levels)
    clu_factor = factor(paste0('cluster', clu), levels = all_levels)
    names(clu_factor) = names(clu)
    
    ## generate confusion matrix
    tb <- table(celltype = ct_factor, clustered_by_topics = clu_factor)
    
    ## scaled by num of cells in each celltype to rule out the case that a large celltype drive the color scheme
    tb <- tb/rowSums(tb)
    flag <- is.na(tb)
    tb2 <- tb[rowMeans(flag)!=1, ]
    flag <- is.na(tb2)
    tb3 <- tb2[,colMeans(flag)!=1]
    tb3 <- tb3[,colMeans(tb3)>0]
    pdf(paste0(getwd(), '/', pdir, 'louvain/',tissue,'.pdf'))
    print(pheatmap(tb3, cluster_rows = T, cluster_cols = T))
    dev.off()
  }, error = function(e) {
    # Optional: handle the error, e.g., print a message
    message("An error occurred", ": ", e$message)
    # The loop will continue with the next iteration
  })
}
