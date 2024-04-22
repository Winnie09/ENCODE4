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
  for (i in 1:ncol(w.hvg.human)){
    for (j in 1:ncol(w.ori.human)){
      topic1 <- w.hvg.human[int,i]
      topic2 <- w.ori.human[int,j]
      
      # Calculate cosine similarity
      similarity <- sum(topic1 * topic2) / (sqrt(sum(topic1^2)) * sqrt(sum(topic2^2)))
      
      # Store the similarity score in the matrix
      cormat[i, j] <- similarity
    
    }
  }
  p1 = pheatmap::pheatmap(t(cormat), scale = 'none', na_col = 'grey',
                          cluster_rows = F,
                          cluster_cols = F)
  
  ###
  w.dsp.human = read.csv(fn3, row.names = 1)
  int = intersect(rownames(w.dsp.human), rownames(w.ori.human))
  cormat2 = matrix(0, nrow= ncol(w.dsp.human), ncol = ncol(w.ori.human))
  dimnames(cormat2) = list(paste0('dsp:', colnames(w.dsp.human)), paste0('ori:', colnames(w.ori.human)))
  for (i in 1:ncol(w.dsp.human)){
    for (j in 1:ncol(w.ori.human)){
      topic1 <- w.dsp.human[int,i]
      topic2 <- w.ori.human[int,j]
      
      # Calculate cosine similarity
      similarity <- sum(topic1 * topic2) / (sqrt(sum(topic1^2)) * sqrt(sum(topic2^2)))
      
      # Store the similarity score in the matrix
      cormat2[i, j] <- similarity
      
    }
  }

  p2 = pheatmap::pheatmap(t(cormat2), scale = 'none', na_col = 'grey',
                          cluster_rows = F,
                          cluster_cols = F)
  pdf(paste0(pdir, tissue, '_cosine_hm.pdf'), width = 12, height = 6)
  grid.arrange(p1$gtable, p2$gtable, ncol=2, widths=c(5,3))
  dev.off()
}
dev.off()




