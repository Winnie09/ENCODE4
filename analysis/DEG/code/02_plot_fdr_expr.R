## plot fdr < 0.05 DEG's expr across general tissues, do clustering 
## 

rdir = '/home/whou10/scratch4/whou10/encode4/DEG/res/'
fdr = readRDS(paste0(rdir, 'fdr_macrophage.rds'))
rownames(fdr) = sub(':.*', '', rownames(fdr))
m = readRDS(paste0(rdir, 'macrophage_expr.rds'))
rownames(m) = sub(':.*', '', rownames(m))
allct = readRDS('/home/whou10/scratch4/whou10/encode4/data/celltype/proc/final.rds')
allct = allct[allct$species == 'Homo sapiens', ]
tmp = allct[allct$celltype=='macrophage', ]
v = sub(':.*','',colnames(m))
gts = tmp[tmp$rna_dataset %in% v, 'generaltissue']
ts = tmp[tmp$rna_dataset %in% v, 'tissue']



pdf(paste0(rdir, 'num_cells_generaltissue.pdf'), width = 10, height = 6)
par(mar=c(15,5,5,5))
barplot(table(gts), las = 2)
dev.off()

pdf(paste0(rdir, 'num_cells_tissue.pdf'), width = 10, height = 6)
par(mar=c(15,5,5,5))
barplot(table(ts), las = 2)
dev.off()



str(fdr)


library(ggplot2)
numsig = sapply(1:ncol(fdr), function(i) sum(fdr[,i] < 0.05))

names(numsig) = colnames(fdr)
pdf(paste0(rdir, 'num_siggene_fdr.pdf'), width = 15, height = 5)
barplot(numsig, names.arg = names(numsig))
dev.off()


res <- sapply(unique(gts), function(i) {
  print(i)
  rowMeans(m[, gts == i])
})

str(res)
i = 1
siggene = unique(as.vector(sapply(1:ncol(fdr), function(i){
  v = fdr[,i]
  rownames(fdr)[order(v)[1:50]]
})))


str(siggene)
m2 = res[siggene,] ###
s = sapply(1:nrow(m2), function(i) which.max(m2[i,])[1])
str(s)
o = order(s)
library(pheatmap)
pdf(paste0(rdir, 'expr_hm.pdf'), height = 45, width = 25)
pheatmap::pheatmap(m2[o,], cluster.rows = F, cluster.cols = T, show_rownames = T)
dev.off()


source('/home/whou10/scratch16/whou10/trajectory_variability/function/01_function.R')

