rm(list=ls())
library(Matrix)
d = readRDS('/home/whou10/scratch4/whou10/encode4/data/celltype/proc/final.rds')
d <- d[!is.na(d$rna_dataset), ]
d$lifestage <-
  ifelse(d$lifestage %in% c('adult', 'child'),
         'adult_child',
         'embryo_postnatal')

l = readRDS('/home/whou10/scratch4/whou10/encode4/data/pseudobulk/snrna/pb/merge/levels.rds')
d2 = cbind(d, level1 = l[,2], level2 = l[,3], level3 = l[,4])


level = 'level1' ##########
rdir <- paste0('/home/whou10/scratch4/whou10/encode4/data/pseudobulk/snrna/pb/', level, '/')
rdir2 <- '/home/whou10/scratch4/whou10/encode4/data/pseudobulk/snrna/pb/merge/'

af = list.files(rdir)
str(af)
if (level == 'level1'){
  enc = sub('-.*','',sub('.*ENC', 'ENC', af))
  spe = d2[match(enc,d2$rna_dataset), 'species']
} else {
  spe = rep('Mus_musculus', length(af))
  spe[grep('Homo_sapiens', af)] = 'Homo_sapiens'
}
table(spe)  
for (uspe in unique(spe)){
  pb <- sapply(af[spe == uspe], function(f){
    print(f)
    tmp = readRDS(paste0(rdir, f))
    
    
    
  })
  str(pb)
  saveRDS(pb, paste0(rdir2, level, '_', sub(' ', '_',uspe), '_count.rds'))
  
  pb = log2(t(t(pb)/(colSums(pb) / 1e6) + 1))
  saveRDS(pb, paste0(rdir2, level, '_', sub(' ', '_', uspe), '_log2CPM.rds'))
}
  




### ==================
##### = plot
### ==================
source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
source('/home/whou10/scratch16/whou10/resource/ggplot_theme.R')
theme_set(.new_theme)
library(ggplot2)
library(RColorBrewer)



for (level in paste0('level', 1:3)){
  print(level)
  for (uspe in c("Homo_sapiens", "Mus_musculus")){
    print(uspe)
    pb = readRDS(paste0(rdir2, level, '_', sub(' ', '_', uspe), '_log2CPM.rds'))
    print(str(pb)) 
  }
}
  


level = 'level1'
uspe = c("Homo_sapiens", "Mus_musculus")[2]
pb = readRDS(paste0(rdir2, level, '_', sub(' ', '_', uspe), '_log2CPM.rds'))
print(str(pb))

pdir = paste0('/home/whou10/scratch4/whou10/encode4/data/pseudobulk/snrna/pb/plot/', level, '/', sub(' ', '_',uspe),'/')
dir.create(pdir)
pr <-
  PCA(
    genebycellmat = pb,
    save.pca = T,
    plot.statistics = F,
    plot.dir = pdir,
    result.dir = pdir,
    PC_for_columns = TRUE,
    findVariableGenes = TRUE,
    findVariableGenesMethod = 'lm',
    maxVariableGenes = 2e3,
    numPC = 50,
    smoothFittingMethod = 'loess'
  )

str(pr)
prfull <- readRDS(paste0(pdir, '/prfull.rds'))

str(prfull)
pr <- prfull$x
pcvar = (prfull[[1]])^2
pcvar <- pcvar/sum(pcvar)
sum(pcvar)
pdf(paste0(pdir, 'PC_variance_explained.pdf'), width = 4, height = 4)
plot(x = seq(1, length(pcvar)), y = cumsum(pcvar), xlab = 'Top PCs', ylab = 'Variance explained', pch = 20)
abline(v = 50, col = 'red')
dev.off()

pr = pr[,1:2]
colnames(pr) <- c(paste0('Principal component 1(', round(pcvar[1]*100,2), '%)'), 
                  paste0('Principal component 2(', round(pcvar[2]*100,2), '%)'))




str(pr)  

if (level == 'level1'){
  rdir <- paste0('/home/whou10/scratch4/whou10/encode4/data/pseudobulk/snrna/pb/', level, '/')
  af = list.files(rdir)
  enc = sub('-.*','',sub('.*ENC', 'ENC', af))
  spe = d2[match(enc,d2$rna_dataset), 'species']
  enc = enc[sub(' ', '_', spe) == uspe]
  u = data.frame(
    pr,
    tissue = d2[match(enc,d2$rna_dataset), 'generaltissue'],
    celltype = d2[match(enc,d2$rna_dataset), 'celltype'],
    lifestage = d2[match(enc,d2$rna_dataset), 'lifestage'])
} else {
  u = data.frame(
    pr,
    tissue = sapply(rownames(pr), function(i)
      strsplit(i, '-')[[1]][3]),
    celltype = sapply(rownames(pr), function(i)
      strsplit(i, '-')[[1]][4]),
    lifestage = sapply(rownames(pr), function(i)
      sub('.rds', '', strsplit(i, '-')[[1]][5])))
}
str(u)  
  
for (var in colnames(u)[3:5]){
  tab <- table(u[,var])
  u[,var] <- paste0(u[,var], '(', tab[u[,var]],')')  
}
str(u)



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



