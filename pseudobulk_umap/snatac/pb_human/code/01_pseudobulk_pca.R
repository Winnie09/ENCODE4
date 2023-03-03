rm(list=ls())
## --------------------------------------
## redo snATAC pseudobuk using more samples
## --------------------------------------
## read in pseudobulk
source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
pdir <- '/home/whou10/scratch4/whou10/encode4/pseudobulk_umap/snatac/pb_human/plot/'
rdir <- '/home/whou10/scratch4/whou10/encode4/pseudobulk_umap/snatac/pb_human/result/'
dir.create(pdir, recursive = T, showWarnings = F)
dir.create(rdir, recursive = T, showWarnings = F)
pd <- readRDS('/home/whou10/scratch4/whou10/encode4/data/pseudobulk/snatac/promoter/pb/human.rds')
str(pd) ## [1:58721, 1:354]
colnames(pd)[1]

## retain pseudobulks with > 50 cells
num_cell <- as.numeric(sapply(colnames(pd), function(i) strsplit(i, ':')[[1]][3]))
summary(num_cell)
pd <- pd[, num_cell > 50]
str(pd) ##  [1:58721, 1:220] 
## pd <- pd[rowMeans(pd > mean(apply(pd, 1, median))) > 0.1, ] ## [1:26561, 1:350]
summary(colMeans(pd))

## retain genes with expression greater than median in >10% pseudobulks
pd <- pd[rowMeans(pd > median(pd[pd > 0], na.rm = T), na.rm = T) > 0.1, ] 
str(pd) ## [1:36016, 1:220]

## PCA and UMAP
pc <- PCA(genebycellmat = pd, save.pca = TRUE, plot.statistics=FALSE, 
          plot.dir = pdir, result.dir = rdir, PC_for_columns = TRUE, 
          findVariableGenes = TRUE, findVariableGenesMethod = 'lm', 
          maxVariableGenes = 2e3, numPC = 20, 
          smoothFittingMethod = 'loess') ## 17 pc selected
u <- UMAP(samplebyfeature_mat = pc, save.umap = TRUE, result.dir = rdir)

## plot PC on var
library(ggplot2)
library(RColorBrewer)
# u <- readRDS(paste0(rdir, '/umap/umap.rds'))
u <- readRDS(paste0(rdir, '/pr.rds'))
str(u)

u <- u[, 1:2]
colnames(u) <- c('Principal component 1', 'Principal component 2')


prfull <- readRDS(paste0(rdir, '/prfull.rds'))
str(prfull)
pcvar = (prfull[[1]])^2
pcvar <- pcvar/sum(pcvar)
sum(pcvar)
head(pcvar, 20)
rownames(u) <- sapply(rownames(u), tools::toTitleCase)
u <-
  data.frame(
    u,
    tissue = as.vector(sapply(sub(':.*', '', rownames(u)), tools::toTitleCase)),
    celltype = sapply(rownames(u), function(i) tools::toTitleCase(strsplit(i, ':')[[1]][2])),
    num_cell = as.numeric(sapply(rownames(u), function(i) strsplit(i,':')[[1]][3])),
    sex = sapply(rownames(u), function(i)
      ifelse(grepl('Male', i, ignore.case = FALSE), 'Male', 'Female')
    ),
    lab = sapply(rownames(u), function(i) strsplit(i,':')[[1]][5])
  )
str(u)
u[,3] <- sapply(u[,3], tools::toTitleCase)
u[,4] <- sapply(u[,4], tools::toTitleCase)
u[,7] <- sapply(u[,7], tools::toTitleCase)
length(unique(u[,3]))
length(unique(u[,4]))
length(unique(u[,6]))
length(unique(u[,7]))
table(u[,7])
table(u[,4])

for (s in c(' Cells$', ' Cell$', ' $', 's$', ' [0-9]$', '[0-9]$', 's$')){
  sel <- grep(s,u[,4])
  sel.rm <- grep('Type', u[,4])
  sel <- setdiff(sel, sel.rm)
  u[sel,4] <- sub(s, '', u[sel,4])  
}
u[u[,4]=='b', 4] <- 'B'
u[u[,4]=='t', 4] <- 'T'
u[u[,4]=='Ce', 4] <- 'CE'
u[u[,4]=='Cf', 4] <- 'CF'
u[u[,4]=='Cm', 4] <- 'CM'
u[u[,4]=='Mtj', 4] <- 'MTJ'
u[u[,4]=='Nmj', 4] <- 'NMF'
u[u[,4]=='Opc', 4] <- 'OPC'
u[u[,4]=='Type2a', 4] <- 'Type2A'
u[u[,4]=='Type2b', 4] <- 'Type2B'
u[u[,4]=='Type2x', 4] <- 'Type2X'
u[u[,4]=='Vlmc', 4] <- 'VLMC'
u[u[,4]=='Lymphatic_ecs', 4] <- 'Lymphatic_ECs'
u[u[,4]=='Fap', 4] <- 'FAP'
u[u[,4]=='Lymphatic_ec', 4] <- 'Lymphatic_EC'
u[u[,4]=='Lymphatic_endothelial', 4] <- 'Lymphatic Endothelial'
u[u[,4]=='Monocytes_macrophage', 4] <- 'Mono_macrophage'

table(u[,4])

u = cbind(u, tissue_short = sub(' Tissue.*', '', sub('.*Sapiens ', '', u[, 3])))
table(u[,8])
saveRDS(u, paste0(rdir, 'umap_and_cell_meta.rds'))

for (var in c('celltype', 'tissue_short', 'lab', 'sex')){
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
    width = 2.8 * length(unique(u[,colid]))/7,
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
myvar <- c('celltype', 'tissue_short', 'lab', 'sex')
for (id in 1:length(myvar)){
  ## assign variable
  assign('var', myvar[id])
  ## define color lable
  v <- paste0(u[,var], ';', clu)
  col.label = rep(NA, nrow(u))
  names(col.label) <- names(v) <- rownames(u)
  for (i in unique(v)){
    tmp = which(v == i, useNames = T)
    col.label[names(tmp)[1]] <- u[names(tmp)[1], var]
  }
  col.label = col.label[!is.na(col.label)]
  df <- cbind(u, col.label = ifelse(rownames(u) %in% names(col.label), u[,var], ""))
  
  
  pdf(paste0(pdir, 'pca_ggrepel_by_', var, '.pdf'), width = ifelse(id ==1, 8, 4), height = ifelse(id ==1, 8, 4))
  
  p <- ggplot(
    data = df,
    aes_string(
      x = colnames(u)[1],
      y = colnames(u)[2],
      color = var,
      label = 'col.label'
    )
  ) +
    geom_point(size = 0.2) +
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
  
  if (length(unique(u[, var])) > 5) {
    p <- p + scale_color_manual(values = colorRampPalette(rev(
      brewer.pal(8, 'Set1')
    ))(length(unique(
      u[, var]
    ))))   
    print(p)
  } else {
    print(p + scale_color_brewer(palette = 'Set1') )
  }
  dev.off()    
}

assign('colvar', 'tissue_short')
pdf(paste0(pdir, 'pca_facet_celltype.pdf'), width = 11, height = 9)
ggplot(data = u, aes_string(x = colnames(u)[1], y = colnames(u)[2], color = colvar)) +
  geom_point( size = 0.2) + 
  scale_color_manual(values = colorRampPalette(rev(brewer.pal(8, 'Set1')))(length(unique(u[, colvar])))) + 
  facet_wrap(~celltype) + 
  theme(legend.position = 'bottom')+
  xlab(paste0('Principal component 1 (', round(pcvar[1],4)*100, '%)')) +
  ylab(paste0('Principal component 2 (', round(pcvar[2],4)*100, '%)'))
dev.off()


pdf(paste0(pdir, 'pca_facet_tissue_short.pdf'), width = 5.5, height = 4.5)
ggplot(data = u, aes_string(x = colnames(u)[1], y = colnames(u)[2], color = colvar)) +
  geom_point( size = 0.2) + 
  scale_color_manual(values = colorRampPalette(rev(brewer.pal(8, 'Set1')))(length(unique(u[, colvar])))) + 
  facet_wrap(~tissue_short) + 
  theme(legend.position = 'none')+
  xlab(paste0('Principal component 1 (', round(pcvar[1],4)*100, '%)')) +
  ylab(paste0('Principal component 2 (', round(pcvar[2],4)*100, '%)'))
dev.off()

