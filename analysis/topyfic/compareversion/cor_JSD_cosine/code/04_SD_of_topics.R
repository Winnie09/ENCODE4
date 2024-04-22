library(grid)
library(gridExtra)
library(pheatmap)
library(philentropy)
source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
pdir = '/home/whou10/scratch4/whou10/encode4/topyfic/compareversion/plot/'
for (tissue in list.files('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_no_downsample/res/')){
  print(tissue)
  fn1 = paste0('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_hvg/res/', tissue, '/figures/gene_weights.csv')
  fn2 = paste0('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_no_downsample/res/', tissue, '/figures/gene_weights.csv')
  fn3 = paste0('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_downsample/res/', tissue, '/figures/gene_weights.csv')
  if (!(file.exists(fn1) & file.exists(fn2) & file.exists(fn3))){
    next
  }
  w.hvg.human = read.csv(fn1, row.names = 1)
  w.ori.human = read.csv(fn2, row.names = 1)
  int = intersect(rownames(w.hvg.human), rownames(w.ori.human))
  cormat = matrix(0, nrow= ncol(w.hvg.human), ncol = ncol(w.ori.human))
  dimnames(cormat) = list(paste0('hvg:', colnames(w.hvg.human)), paste0('ori:', colnames(w.ori.human)))
  sdv1 = rep(0, ncol(w.ori.human))
  for (i in 1:ncol(w.ori.human)){
    sdv1[i] = sd(w.ori.human[,i])
  }
  sdv2 = rep(0, ncol(w.ori.human))
  for (i in 1:ncol(w.ori.human)){
    sdv2[i] = sd(w.ori.human[int,i])
  }
  w.dsp.human = read.csv(fn3, row.names = 1)
  int2 = intersect(rownames(w.dsp.human), rownames(w.ori.human))
  sdv3 = rep(0, ncol(w.ori.human))
  for (i in 1:ncol(w.ori.human)){
    sdv3[i] = sd(w.ori.human[int2,i])
  }
  pdf(paste0(pdir, tissue, '_original_version_gene_weight_sd.pdf'), width = 6.8, height = 2.2)
  par(mfrow = c(1,3))
  hist(sdv1, breaks = 100, main = 'Original version', xlab = 'SD of topics (all genes)')
  hist(sdv2, breaks = 100, main = 'Original version', xlab = 'SD of topics (intersect w/ hvg genes)')
  hist(sdv3, breaks = 100, main = 'Original version', xlab = 'SD of topics (intersect w/ dsp genes)')
  dev.off()
}



