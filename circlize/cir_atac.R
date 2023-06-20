library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(vegan) # for `reorder.hclust` (may be masked by the library `seriation`)
library(dendextend) # for `color_branches`

# ad <- readRDS('/home/zji4/data-hji7/zji/encode/data/celltype/proc/final.rds')
ad <- readRDS('/Users/wenpinhou/Dropbox/encode4/data/celltype/proc/final.rds')
# d <- ad[!is.na(ad$rna_barcode)&!is.na(ad$atac_barcode),]## retains cells with both rna and atac
d <- ad[is.na(ad$rna_barcode)&!is.na(ad$atac_barcode),] ## cells with only atac]
# d <- ad[!is.na(ad$rna_barcode)&is.na(ad$atac_barcode),] # cells with only rna
# d.human <- d[d$species == 'Homo sapiens', ]
colnames(d)[colnames(d)== 'generaltissue'] <- 'celltype.tmp'
colnames(d)[colnames(d)== 'celltype'] <- 'generaltissue'
colnames(d)[colnames(d)== 'celltype.tmp'] <- 'celltype'

at <- table(d$generaltissue, d$celltype, d$species)

t <- as.matrix(as.data.frame.matrix(at[,,1]))
t2 <- t[rowSums(t) > 100,] ## retains celltype with >1e3 cells
# t <- log2(t/max(t)+1) ## control height
t <- t2/rowSums(t2)
t.human <- t
numCell.human <- log2(rowSums(t2))/max(log2(rowSums(t2)))

# t <- as.matrix(as.data.frame.matrix(at[,,2]))
# t3 <- t[rowSums(t) > 100,] ## retains celltype with >1e3 cells
# # t <- log2(t/max(t)+1) ## control height
# t <- t3/rowSums(t3)
# t.mouse <- t
# numCell.mouse <- log2(rowSums(t3))/max(log2(rowSums(t3)))

hc.human=reorder(hclust(dist(t.human)),-as.matrix(t.human)%*%seq(ncol(t.human))^2)

# hc.mouse=reorder(hclust(dist(t.mouse)),-as.matrix(t.mouse)%*%seq(ncol(t.mouse))^2)

labelcolor=hcl(c(260,90,120,60,0,210,180,310)+15,60,70)
barcolor=colorRampPalette(readRDS('/Users/wenpinhou/Dropbox/resource/color74.rds'))(length(unique(c(colnames(t.human)))))
labels=hc.human$labels[hc.human$order]
cut=cutree(hc.human,7) ## from 8 to 6
dend=color_branches(as.dendrogram(hc.human),k=length(unique(cut)),
  col=labelcolor[unique(cut[labels])])

# pdf('/home/zji4/data-hji7/zji/encode/analysis/cellnum/combine/cir.pdf',width=12,height=12)
pdf('/Users/wenpinhou/Dropbox/encode4/data/celltype/proc/cir_atac_human.pdf',width=9,height=9)
circos.clear()
circos.par(cell.padding=c(0,0,0,0))
circos.initialize('a', xlim=c(0,nrow(t.human)))
circos.track(ylim=c(0,1),track.height=.2,track.margin=c(0,0),bg.border=NA,
  panel.fun=function(x,y)for(i in 1:nrow(t.human))circos.text(i-.5,0,labels[i],adj=c(0,.5),
    facing="clockwise",niceFacing=T,cex=.75))
# col=labelcolor[cut[labels[i]]]

circos.track(ylim=c(0,1),track.height=.1,track.margin=c(0,.01),bg.border=NA,
  panel.fun=function(x,y)circos.barplot(as.matrix(t.human)[hc.human$order,],-.5+1:nrow(t.human),
    col=barcolor,bar_width=1,lwd=.3,border="gray20"))

circos.track(ylim=c(1,0),track.height=.1,track.margin=c(0,.01),bg.border=NA,
  panel.fun=function(x,y)circos.barplot(numCell.human[hc.human$order],-.5+1:nrow(t2),
    col='grey86',bar_width=1,lwd=.3,border="grey20"))
lgd_points = Legend(at = colnames(t), type = "points", 
    legend_gp = gpar(col = barcolor), title_position = "topleft", 
    title = "Tissue")

draw(lgd_points, x = unit(1, "mm"), y = unit(1, "mm"), just = c("left", "bottom"))

circos.clear()

dev.off()





