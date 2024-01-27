suppressMessages(library(SingleCellExperiment))
suppressMessages(library(zellkonverter))
library(data.table)
library(Matrix)
dc <- readRDS('/home/whou10/scratch4/whou10/encode4/data/pseudobulk/snrna/pb/Homo_sapiens_pb.rds')
dc <- dc[rowMeans(dc > 2) > 0.2,]
meta <- readRDS('/home/whou10/scratch4/whou10/encode4/data/pseudobulk/snrna/pb/Homo_sapiens_meta.rds')
ctc = sapply(colnames(dc), function(i) strsplit(i, ':')[[1]][2])
sce <- SingleCellExperiment(list(count=as.matrix(dc)),colData=DataFrame(celltype=ctc))
writeH5AD(sce, file = paste0('/home/whou10/scratch4/whou10/encode4/data/scrna/anndata/pseudobulk/h5ad_human.h5ad'))

dc <- readRDS('/home/whou10/scratch4/whou10/encode4/data/pseudobulk/snrna/pb/Mus_musculus_pb.rds')
dc <- dc[rowMeans(dc > 2) > 0.2,]
meta <- readRDS('/home/whou10/scratch4/whou10/encode4/data/pseudobulk/snrna/pb/Mus_musculus_meta.rds')
ctc = sapply(colnames(dc), function(i) strsplit(i, ':')[[1]][2])
sce <- SingleCellExperiment(list(count=as.matrix(dc)),colData=DataFrame(celltype=ctc))
writeH5AD(sce, file = paste0('/home/whou10/scratch4/whou10/encode4/data/scrna/anndata/pseudobulk/h5ad_mouse.h5ad'))

tech <- unique(meta[,5])
for (i in tech){
  print(i)
  dc.tmp <- dc[, meta[,5]==i, drop = F]
  dc.tmp <- dc.tmp[rowMeans(dc.tmp > 2) > 0.2,]
  print(str(dc.tmp))
  ctc = sapply(colnames(dc.tmp), function(j) strsplit(j, ':')[[1]][2])
  sp <- sub(':.*','',colnames(dc.tmp))
  sce <- SingleCellExperiment(list(count=as.matrix(dc.tmp)),colData=DataFrame(celltype=ctc))
  writeH5AD(sce, file = paste0('/home/whou10/scratch4/whou10/encode4/data/scrna/anndata/pseudobulk/','h5ad_mouse_', sub(' .*', '',i), '.h5ad'))
}


