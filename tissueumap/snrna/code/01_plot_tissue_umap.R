## each tissue has many samples, each sample has multiple cells. 
## for each tissue, we obtainedintegrative umap
## plot sc level, color by 1. cell type, 2. study 
library(Seurat)
library(gridExtra)
ct <- readRDS('/home/whou10/data/zji/encode/data/celltype/proc.rds')
ct$cn <- paste0(ct$rna_dataset,':',ct$rna_barcode)
af <- c('mouse-adrenal_gland.rds',
        'mouse-cerebral.rds',
        'mouse-gastrocnemius.rds',
        'mouse-heart.rds',
        'human-colon.rds')

for (f in list.files('/home/whou10/data/zji/encode/data/scrna/tissueseurat_int/')) {
  print(f)
  d <- readRDS(paste0('/home/whou10/data/zji/encode/data/scrna/tissueseurat_int/',f))
  d@meta.data$celltype <- ct[match(colnames(d),ct$cn),'cell_type_name']
  d@meta.data$sample <- sub(':.*','',rownames(d@meta.data))
  d1 <- DimPlot(d,group.by='celltype')
  d2 <- DimPlot(d,group.by='sample')
  if (f %in% af){
    if (length(d@meta.data$celltype) > length(d@meta.data$sample)){
      lay <- matrix(c(rep(1,5), rep(2,4)), nrow = 1)
    } else {
      lay <- matrix(c(rep(1,4), rep(2,5)), nrow = 1)
    }
  } else {
    lay <- matrix(c(rep(1,2), rep(2,2)), nrow = 1)
  }
  pdf(paste0('/home/whou10/scratch4/whou10/encode4/tissueumap/snrna/plot/', sub('.rds','.pdf',f)),width=15,height=4)
  grid.arrange(d1,d2,nrow=1, layout_matrix = lay)
  dev.off()
}



