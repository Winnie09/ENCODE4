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

length(unique(u[,3]))
length(unique(u[,4]))
length(unique(u[,6]))
table(u[,6])



u[u[,4]=='t', 4] <- 'T'
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
    paste0(pdir, 'pca_colored_by_', var, '.pdf'),
    width = 2.8,
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


for (id in 1:length(myvar)){
  assign('var', myvar[id])
  assign('colvar', 'tissue_short')
  pdf(paste0(pdir, 'pca_facet_', var, '.pdf'), width = 11, height = 9)
  print(ggplot(data = u, aes_string(x = colnames(u)[1], y = colnames(u)[2], color = colvar)) +
          geom_point( size = 0.2) + 
          scale_color_manual(values = colorRampPalette(rev(brewer.pal(8, 'Set1')))(length(unique(u[, colvar])))) + 
          facet_wrap(~celltype) + 
          theme(legend.position = 'bottom')+
          xlab(paste0('Principal component 1 (', round(pcvar[1],4)*100, '%)')) +
          ylab(paste0('Principal component 2 (', round(pcvar[2],4)*100, '%)'))
          )
  dev.off()
}


