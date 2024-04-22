rm(list=ls())

##### rockfish
## read in pseudobulk
source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
pdir <- '/home/whou10/scratch4/whou10/encode4/pseudobulk_umap/disco/plot/'
rdir <- '/home/whou10/scratch4/whou10/encode4/pseudobulk_umap/disco/result/'
pd <- readRDS('/home/whou10/scratch4/whou10/encode4/data/pseudobulk/disco/pb/pb.rds')  
str(pd)

## retain pseudobulks with > 50 cells
num_cell <- as.numeric(sub('.*:','',colnames(pd)))
summary(num_cell)
pd <- pd[, num_cell > 50]
pd <- pd[rowMeans(pd > mean(apply(pd, 1, median))) > 0.1, ] ## [1:26561, 1:350]
str(pd)


## PCA and UMAP
pc <- PCA(genebycellmat = pd, save.pca = TRUE, plot.statistics=FALSE, plot.dir = pdir, result.dir = rdir, PC_for_columns = TRUE, findVariableGenes = TRUE, findVariableGenesMethod = 'lm', maxVariableGenes = NULL, numPC = NULL, smoothFittingMethod = 'loess') ## 17 pc selected
## 6 pc selected
u <- UMAP(samplebyfeature_mat = pc, save.umap = TRUE, result.dir = rdir)


## plot UMAP on var
library(ggplot2)
library(RColorBrewer)

u <- readRDS(paste0(rdir, '/umap/umap.rds'))
colnames(u) <- c('UMAP1', 'UMAP2')
u <-
  data.frame(
    u,
    tissue = sub(':.*', '', rownames(u)),
    celltype = sapply(rownames(u), function(i) strsplit(i,':')[[1]][2]),
    num_cell = as.numeric(sapply(rownames(u), function(i) strsplit(i,':')[[1]][3])),
    sex = sapply(rownames(u), function(i)
      ifelse(
        grepl('male', i), 'Male', 'Female'
      ))
  )
u[,3] <- sapply(u[,3], tools::toTitleCase)
u[,4] <- sapply(u[,4], tools::toTitleCase)
for (s in c(' Cells$', ' Cell$', 'cell$', 'cells$', ' $', 's$', ' [0-9]$', '[0-9]$')){
  sel <- grepl(s,u[,4])
  u[sel,4] <- sub(s, '', u[sel,4])  
}

u = cbind(u, tissue_short = sub(' Tissue.*', '', sub('.*Sapiens ', '', u[, 3])))
for (var in c('celltype', 'tissue_short')){
  tab <- table(u[,var])
  u[,var] <- paste0(u[,var], '(', tab[u[,var]],')')  
}

source('/home/whou10/scratch16/whou10/resource/ggplot_theme.R')
theme_set(.new_theme)

for (colid in 3:ncol(u)) {
  assign('var', colnames(u)[colid])
  print(var)
  if (var == 'num_cell') next
  pdf(
    paste0(pdir, 'umap_colored_by_', var, '.pdf'),
    width = 6,
    height = 2
  )
  print(
    ggplot(data = u) +
      geom_point(aes_string(
        x = colnames(u)[1],
        y = colnames(u)[2],
        color = var
      ), size = 0.1) +
      scale_color_manual(values = colorRampPalette(rev(
        brewer.pal(8, 'Set1')
      ))(length(unique(
        u[, var]
      ))))
  )
  dev.off()
}

set.seed(123)
clu = kmeans(u[,1:2], 5)$cluster
var = 'cluster'
pdf(
  paste0(pdir, 'umap_colored_by_cluster.pdf'),
  width = 2.2,
  height = 2
)

ggplot(data = cbind(u,cluster = as.factor(clu))) +
  geom_point(aes_string(
    x = colnames(u)[1],
    y = colnames(u)[2],
    color = var), size = 0.1) +
  scale_color_manual(values = brewer.pal(8, 'Set1'))
dev.off()

library(ggrepel)
myvar <- c('celltype', 'tissue_short')
for (id in 1:length(myvar)){
  ## assign variable
  assign('var', myvar[id])
  ## define color lable
  v <- paste0(u[,var], ';', clu)
  col.label = rep(NA, nrow(u))
  names(col.label) <- names(v) <- rownames(u)
  for (i  in unique(v)){
    tmp = which(v == i, useNames = T)
    col.label[names(tmp)[1]] <- u[names(tmp)[1], var]
  }
  col.label = col.label[!is.na(col.label)]
  df <- cbind(u, col.label = ifelse(rownames(u) %in% names(col.label), u[,var], ""))
  pdf(paste0(pdir, 'umap_ggrepel_by_', var, '.pdf'), width = ifelse(id ==1, 8, 4), height = ifelse(id ==1, 8, 4))
  print(
    ggplot(
      data = df,
      aes_string(
        x = colnames(u)[1],
        y = colnames(u)[2],
        color = var,
        label = 'col.label'
      )
    ) +
      geom_point(size = 0.2) +
      scale_color_manual(values = colorRampPalette(rev(
        brewer.pal(8, 'Set1')
      ))(length(unique(
        u[, var]
      )))) +
      geom_text_repel(
        size = 2,
        max.overlaps = Inf,
        segment.size = 0.1,
        # min.segment.length = 2,
        # nudge_x = .15,
        # box.padding = 0.5,
        # nudge_y = 1,
        # segment.curvature = -1e-20,
        # segment.ncp = 3,
        # segment.angle = 20
      ) +
      # scale_x_continuous(expand = expansion(mult = 0.5)) +
      theme(legend.position = 'none')
  )
  dev.off()    
}

for (id in 1:length(myvar)){
  assign('var', myvar[id])
  assign('colvar', 'tissue_short')
  pdf(paste0(pdir, 'umap_facet_', var, '.pdf'), width = 11, height = 9)
  print(ggplot(data = u, aes_string(x = colnames(u)[1], y = colnames(u)[2], color = colvar)) +
          geom_point( size = 0.2) + 
          scale_color_manual(values = colorRampPalette(rev(brewer.pal(8, 'Set1')))(length(unique(u[, colvar])))) + 
          facet_wrap(~celltype) + 
          theme(legend.position = 'bottom'))
  dev.off()
}

## hist of num_cell in pseudobulk samples
pdf(paste0(pdir, 'hist_pseudobulk_num_cell.pdf'), width = 1.8, height = 1.8)
df <- data.frame(celltype = u[,'celltype'], num_cell = u[,'num_cell'])
ggplot(df, aes(x=num_cell)) + 
  geom_histogram(colour="black", fill="white", binwidth = 50) +
  xlab('Number of cells') + ylab('Count of pseudobulk') 
dev.off()

