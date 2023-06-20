rm(list=ls())
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(vegan) # for `reorder.hclust` (may be masked by the library `seriation`)
library(dendextend) # for `color_branches`

# ad <- readRDS('/home/zji4/data-hji7/zji/encode/data/celltype/proc/final.rds')
ad <- readRDS('/Users/wenpinhou/Dropbox/encode4/data/celltype/proc/final.rds')
d <- ad[!is.na(ad$rna_barcode)&!is.na(ad$atac_barcode),]## retains cells with multiome
# d <- ad[is.na(ad$rna_barcode)&!is.na(ad$atac_barcode),] ## cells with only atac]
# d <- ad[!is.na(ad$rna_barcode)&is.na(ad$atac_barcode),] # cells with only rna
# d.human <- d[d$species == 'Homo sapiens', ]
colnames(d)[colnames(d)== 'generaltissue'] <- 'celltype.tmp'
colnames(d)[colnames(d)== 'celltype'] <- 'generaltissue'
colnames(d)[colnames(d)== 'celltype.tmp'] <- 'celltype'

at <- table(d$generaltissue, d$celltype, d$species)

t <- as.matrix(as.data.frame.matrix(at[,,1]))
# t2 <- t[rowSums(t) > 100,] ## retains celltype with >1e3 cells
# t <- log2(t/max(t)+1) ## control height
t.human <- t/rowSums(t) ## proportion data
t.human[is.na(t.human >= 0)] <- 0
numCell.human <- rowSums(t)/max(rowSums(t))

t <- as.matrix(as.data.frame.matrix(at[,,2]))
# t3 <- t[rowSums(t) > 100,] ## retains celltype with >1e3 cells
# t <- log2(t/max(t)+1) ## control height
t.mouse <- t/rowSums(t)
t.mouse[is.na(t.mouse >= 0)] <- 0
numCell.mouse <- rowSums(t)/max(rowSums(t))

hc.human=reorder(hclust(dist(t.human)),-as.matrix(t.human)%*%seq(ncol(t.human))^2)

hc.mouse=reorder(hclust(dist(t.mouse)),-as.matrix(t.mouse)%*%seq(ncol(t.mouse))^2)

labelcolor=hcl(c(260,90,120,60,0,210,180,310)+15,60,70)
barcolor=colorRampPalette(readRDS('/Users/wenpinhou/Dropbox/resource/color74.rds'))(length(unique(ad$generaltissue)))
names(barcolor) <- sort(unique(ad$generaltissue))
# barcolor=readRDS('/Users/wenpinhou/Dropbox/resource/color74.rds')[1:length(unique(c(colnames(t.human), colnames(t.mouse))))]
labels=hc.human$labels[hc.human$order]
# cut=cutree(hc.human,7) ## from 8 to 6
# dend=color_branches(as.dendrogram(hc.human),k=length(unique(cut)),
#   col=labelcolor[unique(cut[labels])])


## both human an mouse
pdf('/Users/wenpinhou/Dropbox/encode4/circlize/cir_multiome_all.pdf',width=9,height=9)
circos.clear()
circos.par(cell.padding=c(0,0,0,0))
circos.initialize('a', xlim=c(0,nrow(t.mouse)))
circos.track(ylim=c(0,1),track.height=.2,track.margin=c(0,0),bg.border=NA,
  panel.fun=function(x,y)for(i in 1:nrow(t.human))circos.text(i-.5,0,labels[i],adj=c(0,.5),
    facing="clockwise",niceFacing=T,cex=.75))
# col=labelcolor[cut[labels[i]]]

## human
circos.track(ylim=c(0,1),track.height=.1,track.margin=c(0,.01),bg.border=NA,
  panel.fun=function(x,y)circos.barplot(as.matrix(t.human)[hc.human$order,],-.5+1:nrow(t.human),
    col=barcolor,bar_width=1,lwd=1.2,border="grey20"))


## human
circos.track(ylim=c(1,0),track.height=.1,track.margin=c(0,.01),bg.border=NA,
  panel.fun=function(x,y)circos.barplot(numCell.human[hc.human$order],-.5+1:nrow(t.human),
    col='grey20',bar_width=1,lwd=1.2,border="grey20"))

## mouse
circos.track(ylim=c(0,1),track.height=.1,track.margin=c(0,.01),bg.border=NA,
  panel.fun=function(x,y)circos.barplot(as.matrix(t.mouse)[hc.human$order,],-.5+1:nrow(t.mouse),
    col=barcolor,bar_width=1,lwd=1.8,border="skyblue"))

## mouse
circos.track(ylim=c(1,0),track.height=.1,track.margin=c(0,.01),bg.border=NA,
  panel.fun=function(x,y)circos.barplot(numCell.mouse[hc.human$order],-.5+1:nrow(t.mouse),
    col='skyblue',bar_width=1,lwd=1.8,border="skyblue"))

lgd_points = Legend(at = colnames(t.mouse), type = "points", 
    legend_gp = gpar(col = barcolor), title_position = "topleft", 
    title = "Tissue")

draw(lgd_points, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))

circos.clear()

dev.off()


