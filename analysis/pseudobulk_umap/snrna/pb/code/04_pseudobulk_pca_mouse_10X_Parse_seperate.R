rm(list=ls())
## --------------------------------------
## redo RNA pseudobuk using more samples
## --------------------------------------
## read in pseudobulk
source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
tech = 'Parse' # '10X'
pdir <- paste0('/home/whou10/scratch4/whou10/encode4/pseudobulk_umap/snrna/pb/plot/mouse_', tech, '/')
rdir <- paste0('/home/whou10/scratch4/whou10/encode4/pseudobulk_umap/snrna/pb/result/mouse_', tech, '/')
dir.create(pdir)
dir.create(rdir)
pd <- readRDS('/home/whou10/scratch4/whou10/encode4/data/pseudobulk/snrna/pb/Mus_musculus_pb.rds')  
meta <- readRDS('/home/whou10/scratch4/whou10/encode4/data/pseudobulk/snrna/pb/Mus_musculus_meta.rds')  
list.files('/home/whou10/scratch4/whou10/encode4/data/pseudobulk/snrna/pb/')
str(pd)

table(meta$technology)
## select technology: 10X or Parse
str(meta)
pd <- pd[, grepl(tech, meta$technology)]
str(pd)



## retain pseudobulks with > 50 cells
num_cell <- as.numeric(sapply(colnames(pd), function(i) strsplit(i, ':')[[1]][3]))
summary(num_cell)
hist(log10(num_cell))
str(pd)
pd <- pd[, num_cell > 50]
str(pd)
## pd <- pd[rowMeans(pd > mean(apply(pd, 1, median))) > 0.1, ] ## [1:26561, 1:350]
pd <- pd[rowMeans(pd > median(pd[pd > 0])) > 0.1, ] ## [1:17849, 1:323] 
str(pd)

## PCA and UMAP
pc <- PCA(genebycellmat = pd, save.pca = TRUE, plot.statistics=FALSE, 
          plot.dir = pdir, result.dir = rdir, PC_for_columns = TRUE, 
          findVariableGenes = TRUE, findVariableGenesMethod = 'lm', 
          maxVariableGenes = 2e3, numPC = 50, 
          smoothFittingMethod = 'loess') ## 17 pc selected
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
pdf(paste0(pdir, 'PC_variance_explained.pdf'), width = 4, height = 4)
plot(x = seq(1, length(pcvar)), y = cumsum(pcvar), xlab = 'Top PCs', ylab = 'Variance explained', pch = 20)
abline(v = 50, col = 'red')
dev.off()


u <-
  data.frame(
    u,
    generaltissue = tools::toTitleCase(meta[match(sub(':.*','',rownames(u)), meta$dataset),'generaltissue']),
    tissue = tools::toTitleCase(meta[match(sub(':.*','',rownames(u)), meta$dataset),'tissue']),
    celltype = tools::toTitleCase(sapply(rownames(u), function(i) strsplit(i, ':')[[1]][2])),
    num_cell = as.numeric(sapply(rownames(u), function(i) strsplit(i,':')[[1]][3])),
    sex = tools::toTitleCase(meta[match(sub(':.*','',rownames(u)), meta$dataset),'sex']),
    technology = tools::toTitleCase(meta[match(sub(':.*','',rownames(u)), meta$dataset),'technology'])
  )

for (s in c(' Cells$', ' Cell$', ' $', 's$', ' [0-9]$', '[0-9]$')){
  sel <- grep(s,u[,5])
  sel.rm <- grep('Type', u[,5])
  sel <- setdiff(sel, sel.rm)
  u[sel,5] <- sub(s, '', u[sel,5])  
}
u[u[,5]=='t', 5] <- 'T'
u[u[,5]=='b', 5] <- 'B'
table(u[,5])
length(unique(u[,5])) ## 77
saveRDS(u, paste0(rdir, 'umap_and_cell_meta.rds'))


for (var in c('celltype', 'tissue', 'generaltissue')){
  tab <- table(u[,var])
  u[,var] <- paste0(u[,var], '(', tab[u[,var]],')')  
}

source('/home/whou10/scratch16/whou10/resource/ggplot_theme.R')
theme_set(.new_theme)

for (colid in 3:ncol(u)) {
  assign('var', colnames(u)[colid])
  print(var)
  
  if (var == 'num_cell') next
  width = 6 + length(unique(u[,colid]))/20
  height = 2.5
  print(width)
  print(height)
  pdf(
    paste0(pdir, 'pca_colored_by_', var, '.pdf'),
    width = 6 + length(unique(u[,colid]))/20,
    height = 2.5
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
myvar <- c('celltype', 'tissue')
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
  pdf(paste0(pdir, 'pca_ggrepel_by_', var, '.pdf'), width = ifelse(id ==1, 8, 4), height = ifelse(id ==1, 8, 4))
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
      theme(legend.position = 'none') + 
      xlab(paste0('Principal component 1 (', round(pcvar[1],4)*100, '%)')) +
      ylab(paste0('Principal component 2 (', round(pcvar[2],4)*100, '%)'))
  )
  dev.off()    
}

pdf(paste0(pdir, 'pca_facet_celltype.pdf'), width = 13, height = 11)
ggplot(data = u, aes_string(x = colnames(u)[1], y = colnames(u)[2], color = 'tissue')) +
  geom_point( size = 0.2) + 
  scale_color_manual(values = colorRampPalette(rev(brewer.pal(8, 'Set1')))(length(unique(u[, 'tissue'])))) + 
  facet_wrap(~celltype) + 
  theme(legend.position = 'bottom')+
  xlab(paste0('Principal component 1 (', round(pcvar[1],4)*100, '%)')) +
  ylab(paste0('Principal component 2 (', round(pcvar[2],4)*100, '%)'))
dev.off()


pdf(paste0(pdir, 'pca_facet_tissue.pdf'), width = 9, height = 11)
ggplot(data = u, aes_string(x = colnames(u)[1], y = colnames(u)[2], color = 'celltype')) +
  geom_point( size = 0.2) + 
  scale_color_manual(values = colorRampPalette(rev(brewer.pal(8, 'Set1')))(length(unique(u[, 'celltype'])))) + 
  facet_wrap(~tissue) + 
  theme(legend.position = 'bottom')+
  xlab(paste0('Principal component 1 (', round(pcvar[1],4)*100, '%)')) +
  ylab(paste0('Principal component 2 (', round(pcvar[2],4)*100, '%)'))
dev.off()

