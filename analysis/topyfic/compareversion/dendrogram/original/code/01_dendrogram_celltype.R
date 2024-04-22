library(ggplot2)
library(dplyr)
library(RColorBrewer)
# library(tidyverse)
library(ggdendro)
library(plotly)
source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
alltissue = list.files('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_no_downsample/res/')

numcell <- sapply(alltissue, function(tissue){
  ddir <- paste0('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_no_downsample/res/', tissue, '/figures/')
  tb = read.csv(paste0(ddir, 'cell_participation.csv'), row.names = 1)
  nrow(tb)
})

numtopic <- sapply(alltissue, function(tissue){
  ddir <- paste0('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_no_downsample/res/', tissue, '/figures/')
  tb = read.csv(paste0(ddir, 'cell_participation.csv'), row.names = 1)
  ncol(tb)
})

for (tissue in names(sort(numcell))){
  ddir <- paste0('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_no_downsample/res/', tissue, '/figures/')
  pdir <- '/home/whou10/scratch4/whou10/encode4/topyfic/compareversion/dendrogram/original/plot/'
  rdir <- '/home/whou10/scratch4/whou10/encode4/topyfic/compareversion/dendrogram/original/res/'
  tb3 = read.csv(paste0(ddir, 'cell_participation.csv'), row.names = 1)
  
  # pc = PCA(genebycellmat = t(tb3), maxVariableGenes = ncol(tb3))
  d = as.matrix(dist(tb3))
  hclu <- flashClust::flashClust(as.dist(d))
  dend <- as.dendrogram(hclu)

  meta <- readRDS('/home/whou10/scratch4/whou10/encode4/data/celltype/proc/final.rds')
  rn = paste0(meta$rna_dataset, ':', meta$rna_library, ':', meta$rna_barcode)
  meta$label <- rn
  meta <- meta[meta$label %in% hclu$labels, ]
  str(meta)
  
  # Define the colors
  v = unique(meta$celltype)
  color <- colorRampPalette(brewer.pal(n = 8, name = "Set1"))(length(v))
  names(color) = v
  
  # Cluster
  clusters <- cutree(hclu, k = length(v))
  saveRDS(clusters, paste0(rdir, tissue, '_hclu.rds'))

  # extract dendrogram segment data
  dendrogram = dend
  
  dendrogram_data <- dendro_data(dendrogram)
  dendrogram_segments <- dendrogram_data$segments # contains all dendrogram segment data

  # Filter for terminal dendrogram ends
  filtered_dendrogram <- filter(dendrogram_segments, yend == 0)

  # Join with dendrogram_data$labels on 'x'
  joined_dendrogram <- left_join(filtered_dendrogram, dendrogram_data$labels, by = "x")

  # Join with meta on 'label' which is cell names
  dendrogram_ends <- left_join(joined_dendrogram, meta, by = "label")

  # Plot
  p <- ggplot() +
   geom_segment(data = dendrogram_segments,
   aes(x=x, y=y, xend=xend, yend=yend)) +
   geom_segment(data = dendrogram_ends,
   aes(x=x, y=y.x, xend=xend, yend=yend, color = celltype, text = celltype)) + # test aes is for plotly
   scale_color_manual(values = color) +
   scale_y_reverse() +
   coord_flip() + theme_bw() + theme(legend.position = "right") + ylab("Distance")  # flipped x and y coordinates for aesthetic reasons

  png(paste0(pdir, 'dendrogram_celltype/',tissue, '_dendrogram_celltype.png'), res = 300, width = 3500, height = 3000)
  print(p)
  dev.off()
}

