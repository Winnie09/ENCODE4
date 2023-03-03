rm(list=ls())
## check sample size 
# ddir <- '/Users/wenpinhou/Dropbox/encode4/doc/'
# tb = read.csv(paste0(ddir, 'sc_coordination_spreadsheet.csv')) ## 915
#str(tb)
#human <- tb[tb$organism=='human', ] ## 317
#str(human)
#dim(human)
#mouse <- tb[tb$organism=='mouse ', ] ## 161
#dim(mouse)

##### rockfish
#s <- readRDS('/data/hji7/zji/encode/data/scrna/filelist/summary.rds')

## read in pseudobulk
source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
pdir <- '/home/whou10/scratch4/whou10/encode4/pseudobulk_umap/snrna/pb_larger_than_50_cell/plot/'
rdir <- '/home/whou10/scratch4/whou10/encode4/pseudobulk_umap/snrna/pb_larger_than_50_cell/result/'
pd <- readRDS('/home/whou10/scratch4/whou10/encode4/data/pseudobulk/snrna/pb/pb.rds')  
str(pd)

## retain pseudobulks with > 50 cells
num_cell <- as.numeric(sub('.*:','',colnames(pd)))
summary(num_cell)
pd <- pd[, num_cell > 50]
## pd <- pd[rowMeans(pd > mean(apply(pd, 1, median))) > 0.1, ] ## [1:26561, 1:350]
pd <- pd[rowMeans(pd > median(pd[pd > 0])) > 0.1, ] ## [1:17686, 1:231]
str(pd)


## PCA and UMAP
pc <- PCA(genebycellmat = pd, save.pca = TRUE, plot.statistics=FALSE, 
          plot.dir = pdir, result.dir = rdir, PC_for_columns = TRUE, 
          findVariableGenes = TRUE, findVariableGenesMethod = 'lm', 
          maxVariableGenes = 2e3, numPC = 10, 
          smoothFittingMethod = 'loess') ## 17 pc selected
str(pc)
u <- UMAP(samplebyfeature_mat = pc, save.umap = TRUE, result.dir = rdir)

## plot UMAP on var
library(ggplot2)
library(RColorBrewer)

u <- readRDS(paste0(rdir, '/umap/umap.rds'))
colnames(u) <- c('UMAP1', 'UMAP2')
rownames(u) <- sapply(rownames(u), tools::toTitleCase)
u <-
  data.frame(
    u,
    tissue = sapply(sub(':.*', '', rownames(u)), tools::toTitleCase),
    celltype = sapply(sapply(rownames(u), function(i) strsplit(i,':')[[1]][2]) , tools::toTitleCase),
    num_cell = as.numeric(sapply(rownames(u), function(i) strsplit(i,':')[[1]][3])),
    sex = sapply(rownames(u), function(i)
      ifelse(grepl('Male', i, ignore.case = FALSE), 'Male', 'Female')
    )
  )
u[,3] <- sapply(u[,3], tools::toTitleCase)
u[,4] <- sapply(u[,4], tools::toTitleCase)
u[u[,4]=='t', 4] <- 'T'
table(u[,4])

table(u[,6])

for (s in c(' Cells$', ' Cell$', ' $', 's$', ' [0-9]$', '[0-9]$')){
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
    width = 4,
    height = 2
  )
  
  ggplot(data = u) +
    geom_point(aes_string(
      x = colnames(u)[1],
      y = colnames(u)[2],
      color = var
    ), size = 0.1) +
    scale_color_brewer(palette = 'Set1')
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

### DEGs: between liver t cells and other t cells
library(limma)
var <- 'celltype'
sel <- 'T(18)'
expr <-  pd[, u[,var] == sel]
str(expr)
var2 <- 'tissue_short'
des <- cbind(1,ifelse(grepl('Liver', u[colnames(expr),var2]),1,0))
res <- topTable(eBayes(lmFit(expr, des)),n=nrow(expr),coef=2)
res <- res[res[,'adj.P.Val']<0.05,]
gs <- rownames(res)
str(gs)
write.csv(res, paste0(rdir, '/DEG_of_t_cells_between_liver_and_others.csv'))
saveRDS(res, paste0(rdir, '/DEG_of_t_cells_between_liver_and_others.rds'))

### GO enrichment
sigres <- myGO(sub(':.*', '',gs), sub(':.*', '', rownames(expr)))
sigres <- sigres[sigres[,'FDR'] < 0.05 & sigres[,'FC'] > 2, ]
str(sigres)
saveRDS(sigres, paste0(rdir, '/DEG_of_t_cells_between_liver_and_others_GOEnrich.rds'))
write.csv(sigres, paste0(rdir, '/DEG_of_t_cells_between_liver_and_others_GOEnrich.csv'))

## hist of num_cell in pseudobulk samples
pdf(paste0(pdir, 'hist_pseudobulk_num_cell.pdf'), width = 1.6, height = 1.3)
df <- data.frame(celltype = u[,'celltype'], num_cell = u[,'num_cell'])
ggplot(df, aes(x=num_cell)) + 
  geom_histogram(colour="black", fill="white", binwidth = 50) +
  xlab('Number of cells') + ylab('Count of pseudobulk') +
  xlim(c(0, max(df$num_cell) + 1000))
dev.off()

