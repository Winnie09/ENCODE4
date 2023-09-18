# tissue = 'Homo_sapiens:adrenal_gland'
tissue = commandArgs(trailingOnly = T)[[1]][1]

library(grid)
library(gridExtra)
library(pheatmap)
source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
pdir = '/home/whou10/scratch4/whou10/encode4/topyfic/compareversion/plot/'
w.hvg.human = read.csv('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_hvg/res/Homo_sapiens:adrenal_gland/figures/gene_weights.csv', row.names = 1)
w.ori.human = read.csv('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_no_downsample/res/Homo_sapiens:adrenal_gland/figures/gene_weights.csv', row.names = 1)
int = intersect(rownames(w.hvg.human), rownames(w.ori.human))
cormat = matrix(0, nrow= ncol(w.hvg.human), ncol = ncol(w.ori.human))
dimnames(cormat) = list(paste0('hvg:', colnames(w.hvg.human)), paste0('ori:', colnames(w.ori.human)))
for (i in 1:ncol(w.hvg.human)){
  for (j in 1:ncol(w.ori.human)){
    cormat[i,j] = cor(w.hvg.human[int,i], w.ori.human[int,j])
  }
}
p1 = pheatmap::pheatmap(t(cormat), scale = 'none', na_col = 'grey',
                   cluster_rows = F,
                   cluster_cols = F)

####
w.dsp.human = read.csv('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_downsample/res/Homo_sapiens:adrenal_gland/figures/gene_weights.csv', row.names = 1)
int = intersect(rownames(w.dsp.human), rownames(w.ori.human))
cormat2 = matrix(0, nrow= ncol(w.dsp.human), ncol = ncol(w.ori.human))
dimnames(cormat2) = list(paste0('dsp:', colnames(w.dsp.human)), paste0('ori:', colnames(w.ori.human)))
for (i in 1:ncol(w.dsp.human)){
  for (j in 1:ncol(w.ori.human)){
    cormat2[i,j] = cor(w.dsp.human[int,i], w.ori.human[int,j])
  }
}
p2 = pheatmap::pheatmap(t(cormat2), scale = 'none', na_col = 'grey',
                        cluster_rows = F,
                        cluster_cols = F)
pdf(paste0(pdir, tissue, '_corhm.pdf'), width = 12, height = 6)
grid.arrange(p1$gtable, p2$gtable, ncol=2, widths=c(5,3))
dev.off()
