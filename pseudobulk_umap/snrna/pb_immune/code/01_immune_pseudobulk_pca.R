rm(list=ls())
## ------------------------------------------
## do snRNA pseudobuk PCA using immune cells
## ------------------------------------------
## read in pseudobulk
source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
pdir <- '/home/whou10/scratch4/whou10/encode4/pseudobulk_umap/snrna/pb_immune/plot/'
rdir <- '/home/whou10/scratch4/whou10/encode4/pseudobulk_umap/snrna/pb_immune/result/'
pd <- readRDS('/home/whou10/scratch4/whou10/encode4/data/pseudobulk/snrna/pb/human.rds')  
colnames(pd) <- sapply(colnames(pd), tools::toTitleCase)
meta <- readRDS(paste0(rdir, 'umap_and_cell_meta.rds'))

## ---------------------------------------
## subset immune cells for downstream PCA
## ---------------------------------------
selcell <- rownames(meta)[grep(('(Immune)|(B\\()|(Lymphatic)|(Lymphoid)|(T\\()|(Macrophage)|(Myeloid)') , meta[,4])]
pd <- pd[, selcell]
str(pd)

## retain genes with expression greater than median in >10% pseudobulks
pd <- pd[rowMeans(pd > median(pd[pd > 0])) > 0.1, ] ##  [1:15681, 1:65]
str(pd)

## PCA and UMAP
pc <- PCA(genebycellmat = pd, save.pca = TRUE, plot.statistics=FALSE, 
          plot.dir = pdir, result.dir = rdir, PC_for_columns = TRUE, 
          findVariableGenes = TRUE, findVariableGenesMethod = 'lm', 
          maxVariableGenes = 2e3, numPC = 10, 
          smoothFittingMethod = 'loess') ## 
u <- UMAP(samplebyfeature_mat = pc, save.umap = TRUE, result.dir = rdir)

## plot PC on var
library(ggplot2)
library(RColorBrewer)
# u <- readRDS(paste0(rdir, '/umap/umap.rds'))
u <- readRDS(paste0(rdir, '/pr.rds'))
u <- u[, 1:2]
colnames(u) <- c('Principal component 1', 'Principal component 2')

prfull <- readRDS(paste0(rdir, '/prfull.rds'))
str(prfull)
pcvar = (prfull[[1]])^2
pcvar <- pcvar/sum(pcvar)
sum(pcvar)
rownames(u) <- sapply(rownames(u), tools::toTitleCase)
str(meta)
u <-
  data.frame(
    u,
    meta[rownames(u),3:ncol(meta)]
  )

for (var in c('celltype', 'tissue_short')){
  u[,var] <- sub('\\(.*','', u[,var])
  tab <- table(u[,var])
  u[,var] <- paste0(u[,var], '(', tab[u[,var]],')')  
}

str(u)

source('/home/whou10/scratch16/whou10/resource/ggplot_theme.R')
theme_set(.new_theme)

for (colid in 3:ncol(u)) {
  assign('var', colnames(u)[colid])
  print(var)
  if (var == 'num_cell') next
  pdf(
    paste0(pdir, 'pca_colored_by_', var, '.pdf'),
    width = 2.8 * length(unique(u[,colid]))/9,
    height = 2
  )
  print(
    ggplot(data = u) +
      geom_point(aes_string(
        x = colnames(u)[1],
        y = colnames(u)[2],
        color = var
      ), size = 0.1) +
      xlab(paste0('Principal component 1 (', round(pcvar[1],4)*100, '%)')) +
      ylab(paste0('Principal component 2 (', round(pcvar[2],4)*100, '%)')) +
      # scale_color_brewer(palette = 'Set1')
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
  paste0(pdir, 'pca_colored_by_cluster.pdf'),
  width = 2.2,
  height = 2
)
ggplot(data = cbind(u,cluster = as.factor(clu))) +
  geom_point(aes_string(
    x = colnames(u)[1],
    y = colnames(u)[2],
    color = var), size = 0.1) +
  scale_color_manual(values = brewer.pal(8, 'Set1'))+
  xlab(paste0('Principal component 1 (', round(pcvar[1],4)*100, '%)')) +
  ylab(paste0('Principal component 2 (', round(pcvar[2],4)*100, '%)'))
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
  pdf(paste0(pdir, 'pca_ggrepel_by_', var, '.pdf'), width = ifelse(id ==1, 3.5, 3.5), height = ifelse(id ==1, 3.5, 3.5))
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
      geom_point(size = 0.6) +
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
      theme(legend.position = 'none') + 
      xlab(paste0('Principal component 1 (', round(pcvar[1],4)*100, '%)')) +
      ylab(paste0('Principal component 2 (', round(pcvar[2],4)*100, '%)'))
  )
  dev.off()    
}

pdf(paste0(pdir, 'pca_facet_celltype.pdf'), width = 5, height = 4.5)
ggplot(data = u, aes_string(x = colnames(u)[1], y = colnames(u)[2], color = 'tissue_short')) +
  geom_point( size = 1) + 
  scale_color_manual(values = colorRampPalette(rev(brewer.pal(8, 'Set1')))(length(unique(u[, 'tissue_short'])))) + 
  facet_wrap(~celltype) + 
  theme(legend.position = 'bottom')+
  xlab(paste0('Principal component 1 (', round(pcvar[1],4)*100, '%)')) +
  ylab(paste0('Principal component 2 (', round(pcvar[2],4)*100, '%)'))
dev.off()


pdf(paste0(pdir, 'pca_facet_tissue_short.pdf'), width = 5, height = 4.5)
ggplot(data = u, aes_string(x = colnames(u)[1], y = colnames(u)[2], color = 'tissue_short')) +
  geom_point( size = 1) + 
  scale_color_manual(values = colorRampPalette(rev(brewer.pal(8, 'Set1')))(length(unique(u[, 'tissue_short'])))) + 
  facet_wrap(~tissue_short) + 
  theme(legend.position = 'none')+
  xlab(paste0('Principal component 1 (', round(pcvar[1],4)*100, '%)')) +
  ylab(paste0('Principal component 2 (', round(pcvar[2],4)*100, '%)'))
dev.off()

