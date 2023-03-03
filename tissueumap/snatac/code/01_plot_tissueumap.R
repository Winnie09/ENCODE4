## each tissue has many samples, each sample has multiple cells. 
## for each tissue, we obtainedintegrative umap
## plot sc level, color by 1. cell type, 2. study 
library(Seurat)
library(gridExtra)
ct <- readRDS('/home/whou10/data/zji/encode/data/celltype/proc.rds')
ct$cn <- paste0(ct$atac_dataset,':',ct$atac_barcode)
af <- c('mouse-adrenal_gland.rds',
        'mouse-cerebral.rds',
        'mouse-gastrocnemius.rds',
        'mouse-heart.rds',
        'human-colon.rds')
for (f in list.files('/home/whou10/data/zji/encode/data/snatac/peak/tissueseurat_noint/')) { ## direct PCA without sample integreation
  print(f)
  d <- readRDS(paste0('/home/whou10/data/zji/encode/data/snatac/peak/tissueseurat_noint/',f))
  d@meta.data$celltype <- ct[match(colnames(d),ct$cn),'cell_type_name']
  d@meta.data$sample <- sub(':.*','',rownames(d@meta.data))
  d1 <- DimPlot(d,group.by='celltype')
  d2 <- DimPlot(d,group.by='sample')
  if (f %in% af){
    if (length(unique(d@meta.data$celltype)) > 1.5*length(unique(d@meta.data$sample))){
      lay <- matrix(c(rep(1,6), rep(2,4)), nrow = 1)
    } else {
      lay <- matrix(c(rep(1,4), rep(2,6)), nrow = 1)
    }
  } else {
    lay <- matrix(c(rep(1,2), rep(2,2)), nrow = 1)
  }
  pdf(paste0('/home/whou10/scratch4/whou10/encode4/tissueumap/snatac/plot/', sub('.rds','.pdf',f)),width=15,height=5)
  grid.arrange(d1,d2,nrow=1, layout_matrix = lay)
  dev.off()
}




