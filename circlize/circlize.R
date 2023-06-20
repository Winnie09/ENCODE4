library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(vegan) # for `reorder.hclust` (may be masked by the library `seriation`)
library(dendextend) # for `color_branches`

# ad <- readRDS('/home/zji4/data-hji7/zji/encode/data/celltype/proc/final.rds')
ad <- readRDS('/Users/wenpinhou/Dropbox/encode4/data/celltype/proc/final.rds')
d <- ad[!is.na(ad$rna_barcode)&!is.na(ad$atac_barcode),] ## retains cells with both rna and atac
t <- table(d$generaltissue,d$celltype)
t <- as.matrix(as.data.frame.matrix(t))
t <- t[,colSums(t) > 1000] ## retains celltype with >1e3 cells
t <- t/max(t)/3 ## control height
tmulti <- t

d <- ad[is.na(ad$rna_barcode)&!is.na(ad$atac_barcode),] ## cells with only atac
t <- table(d$generaltissue,d$celltype)
t <- as.matrix(as.data.frame.matrix(t))
t <- t[,colSums(t) > 1000]
t <- t/max(t)/3
tatac <- tmulti

d <- ad[!is.na(ad$rna_barcode)&is.na(ad$atac_barcode),]
t <- table(d$generaltissue,d$celltype)
t <- as.matrix(as.data.frame.matrix(t))
t <- t[,colSums(t) > 1000]
t <- t/max(t)/3
trna <- tmulti

hc=reorder(hclust(dist(t)),-as.matrix(t)%*%seq(ncol(t))^2)

labelcolor=hcl(c(260,90,120,60,0,210,180,310)+15,60,70)
barcolor=rainbow(ncol(t))

labels=hc$labels[hc$order]
cut=cutree(hc,6) ## from 8 to 6
dend=color_branches(as.dendrogram(hc),k=length(unique(cut)),
  col=labelcolor[unique(cut[labels])])

# pdf('/home/zji4/data-hji7/zji/encode/analysis/cellnum/combine/cir.pdf',width=12,height=12)
pdf('/Users/wenpinhou/Dropbox/encode4/data/celltype/proc/cir.pdf',width=11,height=11)

circos.clear()
circos.par(cell.padding=c(0,0,0,0))
circos.initialize("a",xlim=c(0,nrow(t)))

circos.track(ylim=c(0,1),track.height=.2,track.margin=c(0,0),bg.border=NA,
  panel.fun=function(x,y)for(i in 1:nrow(t))circos.text(i-.5,0,labels[i],adj=c(0,.5),
    facing="clockwise",niceFacing=T,cex=.75))
#col=labelcolor[cut[labels[i]]]

circos.track(ylim=c(0,1),track.height=.3,track.margin=c(0,.01),bg.border=NA,
  panel.fun=function(x,y)circos.barplot(as.matrix(t)[hc$order,],-.5+1:nrow(t),
    col=barcolor,bar_width=1,lwd=.3,border="gray20"))

circos.track(ylim=c(0,1),track.height=.3,track.margin=c(0,.01),bg.border=NA,
  panel.fun=function(x,y)circos.barplot(as.matrix(t)[hc$order,],-.5+1:nrow(t),
    col=barcolor,bar_width=1,lwd=.3,border="gray20"))

circos.track(ylim=c(0,1),track.height=.3,track.margin=c(0,.01),bg.border=NA,
  panel.fun=function(x,y)circos.barplot(as.matrix(t)[hc$order,],-.5+1:nrow(t),
    col=barcolor,bar_width=1,lwd=.3,border="gray20"))


lgd_points = Legend(at = colnames(t), type = "points", 
    legend_gp = gpar(col = barcolor), title_position = "topleft", 
    title = "Cell types")

draw(lgd_points, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))


circos.clear()

dev.off()


