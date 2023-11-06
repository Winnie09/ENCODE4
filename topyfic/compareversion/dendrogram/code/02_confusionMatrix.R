library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
setwd('/home/whou10/scratch4/whou10/encode4/topyfic/')

ddir <- commandArgs(trailingOnly = T)[1]
rdir <- commandArgs(trailingOnly = T)[2]
pdir <- commandArgs(trailingOnly = T)[3]

alltissue = list.files(ddir)

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

for (tissue in sort(names(numcell))){
  tryCatch({
    # Code that might throw an error
    print(tissue)
    dir.create(pdir, showWarnings = F, recursive = T)
    clu = readRDS(paste0(rdir, tissue, '_hclu.rds'))
    
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
    pdf(paste0(pdir, 'hclu/', tissue,'.pdf'))
    print(pheatmap(tb3, cluster_rows = F, cluster_cols = F))
    dev.off()
  }, error = function(e) {
    # Optional: handle the error, e.g., print a message
    message("An error occurred ", ": ", e$message)
    # The loop will continue with the next iteration
  })
}
