at = list.files('/home/whou10/scratch4/whou10/encode4/topyfic/tissue/res/')

for (t in at){
  print(t)
  setwd(paste0('/home/whou10/scratch4/whou10/encode4/topyfic/tissue/res/', t, '/figures'))
  a = read.csv('gene_weights.csv', row.names = 1)
  
  g = sapply(1:ncol(a), function(j){
    v = a[,j]
    names(v) = rownames(a)
    names(sort(v, decreasing = T)[1:1e2])
  })
  unig = unique(as.vector(g))
  
  g = g[order(rowSums(a[g,]))]
  
  library(pheatmap)
  pdf('gene_by_topic_hm.pdf', height = 18, width = 15)
  pheatmap(a[unig,], cluster_rows = F, cluster_cols = F, show_rownames = T, show_colnames = T, scale = 'none', fontsize_row = 5)
  dev.off()
}
  
