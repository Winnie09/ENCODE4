library(Seurat)
d1 <- readRDS('adrenal_Parse_10x_to_integrate.rds')
sapply(d1,function(i) table(i@meta.data$technology))
sapply(d1,function(i) table(i@meta.data$depth1))

d2 <- readRDS('cortex_Parse_10x_to_integrate.rds')
sapply(d2,function(i) table(i@meta.data$technology))
sapply(d2,function(i) table(i@meta.data$depth1))

d3 <- readRDS('gastrocnemius_Parse_10x_to_integrate.rds')
sapply(d3,function(i) table(i@meta.data$technology))
sapply(d3,function(i) table(i@meta.data$depth1))

d4 <- readRDS('heart_Parse_10x_to_integrate.rds')
sapply(d4,function(i) table(i@meta.data$technology))
sapply(d4,function(i) table(i@meta.data$depth1))

d5 <- readRDS('hippocampus_Parse_10x_to_integrate.rds')
sapply(d5,function(i) table(i@meta.data$technology))
sapply(d5,function(i) table(i@meta.data$depth1))

# d6 <- readRDS('psoas_integrated.rds')
# table(d6@meta.data$technology)
# d7 <- readRDS('ureter.rds')
# table(d7@meta.data$technology)

d <- list(d1[[2]],d2[[2]],d3[[2]],d4[[2]],d5[[2]])
library(Matrix)
cm <- sapply(c(1,2,4,5),function(i) {
  ct <- apply(d[[i]]@meta.data[,c('tissue','celltypes')],1,paste0,collapse=':')
  pb <- rowsum(t(as.matrix(d[[i]]@assays$RNA@counts)),ct)
  log2(t(pb/rowSums(pb)*1e6)+1)
},simplify = F)
cm <- do.call(cbind,cm)
saveRDS(cm,file='/home/zji4/data-hji7/zji/encode/data/ambient/parseshallow.rds')

d <- list(d1[[1]],d2[[1]],d3[[1]],d4[[1]],d5[[1]])
library(Matrix)
cm <- sapply(c(1,2,4,5),function(i) {
  ct <- apply(d[[i]]@meta.data[,c('tissue','celltypes')],1,paste0,collapse=':')
  pb <- rowsum(t(as.matrix(d[[i]]@assays$RNA@counts)),ct)
  log2(t(pb/rowSums(pb)*1e6)+1)
},simplify = F)
cm <- do.call(cbind,cm)
saveRDS(cm,file='/home/zji4/data-hji7/zji/encode/data/ambient/10x.rds')

