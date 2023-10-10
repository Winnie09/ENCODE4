rm(list=ls())
options(scipen = 9999999)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(vegan) # for `reorder.hclust` (may be masked by the library `seriation`)
library(dendextend) # for `color_branches`
library(ggplot2)
source('/Users/wenpinhou/Dropbox/resource/ggplot_theme.R')
# ad <- readRDS('/home/zji4/data-hji7/zji/encode/data/celltype/proc/final.rds')
ad.bk <- ad <- readRDS('/Users/wenpinhou/Dropbox/encode4/data/celltype/proc/final.rds')
length(unique(ad$celltype)) ## 121

## check celltypes
unique(ad.bk[grep('mesothelial', ad.bk$celltype),c('species','celltype','tissue','generaltissue','CL')])

matchtb = read.csv('/Users/wenpinhou/Dropbox/encode4/circlize/circlize_plot_summarize_celltype.csv', sep = ',')
matchtb[,1] = tolower(matchtb[,1])
matchtb[matchtb[,1] %in% c('hippocampus', 'layer of hippocampus', 'left cerebral cortex'),1] = 'brain'
                                  
## add cell subtype numbers to cell type names
ctv = unique(matchtb[,4])
for (i in 1:length(ctv)){
  print(i)
  len1 = length(which(matchtb[,4] == ctv[i] & matchtb[,2] == 'Homo sapiens')) 
  len2 = length(which(matchtb[,4] == ctv[i] & matchtb[,2] == 'Mus musculus'))
  if (length(which(ad[,7] == ctv[i] & ad[,1] == 'Homo sapiens'))>0) len1 <- len1 + 1 # count in original cell types
  if (length(which(ad[,7] == ctv[i] & ad[,1] == 'Mus musculus'))>0) len2 <- len2 + 1 # count in original cell types
  # if (len1 == 0) len1 = 1
  # if (len2 == 0)  {
  #   if (length(which(ad[,7] == ctv[i] & matchtb[,1] == 'Mus musculus'))>0) len2 = 1
  # }
  if (len1 == 0 | len2 == 0) print('still 0 ...')
  text = paste0('(',len1,
                         ',',len2,')')
  print(paste0(sub('\\(.*','',ctv[i]), text))
  ad[ad$celltype == ctv[i], 'celltype']  <- matchtb[matchtb[,4] == ctv[i], 4] <- paste0(sub('\\(.*','',ctv[i]), text)
  Sys.sleep(0.2)
} 


for (i in 1:nrow(matchtb)){
  if (matchtb[i,2] == 'Mus musculus'){
    ## if in mouse, the tissue can match "tissue" or "generaltissue"
    if (length(which(ad$generaltissue==matchtb[i,1] & ad$species == 'Mus musculus' & ad$celltype == (matchtb[i,3]))) > 0){
    ad[which(ad$generaltissue==matchtb[i,1] & ad$species == 'Mus musculus' & ad$celltype == (matchtb[i,3])), 'celltype'] = matchtb[i,4]    
    } else if (length(which(ad$tissue==matchtb[i,1] & ad$species == 'Mus musculus' & ad$celltype == (matchtb[i,3]))) > 0){
    ad[which(ad$tissue==matchtb[i,1] & ad$species == 'Mus musculus' & ad$celltype == (matchtb[i,3])), 'celltype'] = matchtb[i,4]    
    } else {
      print('dont find.')
      print(i)
      print(matchtb[i,])
    }
  } else {
    ## if in human, the tissue match is omitted. There is an "all" in tissue names. 
    ad[which(ad$species == 'Homo sapiens' & ad$celltype == (matchtb[i,3])), 'celltype'] = matchtb[i,4] 
  }
}
length(unique(ad$celltype)) ## 80

## check celltypes
# unique(ad[grep('monocyte', ad$celltype),c('celltype','tissue','generaltissue','CL')])

## check human and mouse cells difference
# c1 = unique(ad[ad$species == 'Mus musculus', 'celltype'])
# c2 = unique(ad[ad$species == 'Homo sapiens', 'celltype'])
# sort(setdiff(c1, c2))




##
tech <- rep('sc-multiome', nrow(ad))
tech[is.na(ad$rna_barcode)&!is.na(ad$atac_barcode)] <- 'snATAC-seq'
tech[!is.na(ad$rna_barcode)&is.na(ad$atac_barcode)] <- 'scRNA-seq'
ad <- data.frame(ad, technology = tech)
d <- ad ## use all cells


at <- table(d$celltype, d$generaltissue, d$species)
t1 <- as.matrix(as.data.frame.matrix(at[,,1]))
# t2 <- t[rowSums(t) > 100,] ## retains celltype with >1e3 cells
# t <- log2(t/max(t)+1) ## control height
t.human <- t1/rowSums(t1) ## proportion data
t.human[is.na(t.human >= 0)] <- 0
numCell.human <-log10(rowSums(t1)+1)/max(log10(rowSums(t1)+1))


t2 <- as.matrix(as.data.frame.matrix(at[,,2]))
# t3 <- t[rowSums(t) > 100,] ## retains celltype with >1e3 cells
# t <- log2(t/max(t)+1) ## control height
t.mouse <- t2/rowSums(t2)
t.mouse[is.na(t.mouse >= 0)] <- 0
numCell.mouse <- log10(rowSums(t2)+1)/max(log10(rowSums(t2)+1))



## order the celltypes:[note: the length of values is 79, but the reording numberic vector has length 80, the result order does not put brain celltypes at the end]
hc.human=reorder(hclust(dist(t.human)),-as.matrix(t.human)%*%seq(ncol(t.human))^2)
# hc.human=reorder(hclust(dist(hc.human)),-as.matrix(t.mouse[hc.human$labels,])%*%seq(ncol(t.human))^2)

#overallorder <- rev(names(sort((-as.matrix(t.human)%*%seq(ncol(t.human))^2)[,1]))) ###### !!


#hclust(dist(t.human))$order

overallorder <- hc.human$labels[hc.human$order]
mouseterm <- names(which(rowSums(t.human)==0))
overallorder <- c(overallorder[overallorder%in%mouseterm],overallorder[!overallorder%in%mouseterm])
###### !!



#########################
at.tech <- table(d$celltype, d$technology, d$species)
t.tech <- as.matrix(as.data.frame.matrix(at.tech[,,1]))
# t2 <- t[rowSums(t) > 100,] ## retains celltype with >1e3 cells
# t <- log2(t/max(t)+1) ## control height
t.human.tech <- t.tech/rowSums(t.tech) ## proportion data
t.human.tech[is.na(t.human.tech >= 0)] <- 0


t.tech <- as.matrix(as.data.frame.matrix(at.tech[,,2]))
# t3 <- t[rowSums(t) > 100,] ## retains celltype with >1e3 cells
# t <- log2(t/max(t)+1) ## control height
t.mouse.tech <- t.tech/rowSums(t.tech)
t.mouse.tech[is.na(t.mouse.tech >= 0)] <- 0
########################

## define and order bar colors according to show up order
labelcolor=hcl(c(260,90,120,60,0,210,180,310)+15,60,70)
# barcolor=colorRampPalette(readRDS('/Users/wenpinhou/Dropbox/resource/color74.rds'))(length(unique(ad$generaltissue[!is.na(ad$generaltissue)])))
# names(barcolor) <- sort(unique(ad$generaltissue))

coldf2 = readRDS('/Users/wenpinhou/Dropbox/resource/color32.rds')
colv = coldf2[,2]
barcolor=colv[1:length(unique(ad$generaltissue[!is.na(ad$generaltissue)]))]
names(barcolor) <- sort(unique(ad$generaltissue))
techcolor=RColorBrewer::brewer.pal(n = 6, 'Set3')[4:6]
names(techcolor)=unique(d$technology)
orange=RColorBrewer::brewer.pal(n = 5, 'Oranges')[2:5]

write.csv(data.frame(Tissue = names(barcolor), Color = barcolor, stringsAsFactors = F), '/Users/wenpinhou/Dropbox/encode4/circlize/color_scheme.csv', row.names = F)


tmp = read.csv('/Users/wenpinhou/Dropbox/encode4/circlize/color_scheme.csv')
str(tmp)
barcolor = tmp[,2]
names(barcolor) = tmp[,1]


colm = reshape2::melt(t.human[overallorder, ])
colm = colm[colm[,3]!=0,]
colm[,1] = factor(colm[,1],levels=overallorder)
colm <- colm[rev(order(as.numeric(colm[,1]))),]
colm[,2] = as.character(colm[,2])

barcolor = c(barcolor[unique(colm[,2])], barcolor[setdiff(names(barcolor), unique(colm[,2]))])


##
techcolor=RColorBrewer::brewer.pal(n = 12, 'Set3')[c(4,5,8)]
names(techcolor)=unique(d$technology)
orange=RColorBrewer::brewer.pal(n = 5, 'Oranges')[2:5]


# barcolor=readRDS('/Users/wenpinhou/Dropbox/resource/color74.rds')[1:length(unique(c(colnames(t.human), colnames(t.mouse))))]
# labels=hc.human$labels[hc.human$order]



#(-as.matrix(t.human)%*%seq(ncol(t.human))^2)[,1]
#(-as.matrix(t.mouse)%*%seq(ncol(t.mouse))^2)[,1]

# Convert the matrix to a data frame for ggplot2
data_df.human <- as.data.frame(as.table(t.human))
colnames(data_df.human) <- c("Celltype", "Tissue", "Proportion")
data_df.human$Celltype <-factor(as.character(data_df.human$Celltype),levels=overallorder) ###### !!
# data_df.human$Proportion <- data_df.human$Proportion * 100
p1 <- ggplot(data_df.human, aes(y = Celltype, x = Proportion, fill = Tissue)) +
  geom_bar(stat = "identity") +
  labs(title = "Human: tissue proportion", y = "", x = "Proportion") +
  #scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(values = barcolor, breaks = names(barcolor))+
  theme( legend.position = 'none', plot.title = element_text(size=8))

tech.human <- as.data.frame(as.table(t.human.tech))
colnames(tech.human) <- c("Celltype", "Technology", "Proportion")
tech.human$Celltype <-factor(as.character(tech.human$Celltype),levels=overallorder) ###### !!
# tech.human$Proportion <- tech.human$Proportion * 100
p2 <- ggplot(tech.human, aes(y = Celltype, x = Proportion, fill = Technology)) +
  geom_bar(stat = "identity") +
  labs(title = "Human: sequencing technology proportion", y = "Celltype", x = "Proportion") +
  #scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(values = techcolor)+
  theme(axis.text.y = element_blank())

num_df.human = data.frame(Celltype = names(numCell.human), 
                    Number = numCell.human, 
                    stringsAsFactors = T)
levels(num_df.human[,1]) = names(numCell.human)


num_df.human$Celltype <-factor(as.character(num_df.human$Celltype),levels=overallorder) ###### !!



p3 <- ggplot(num_df.human, 
       aes(y = Celltype, x = Number, fill = 'grey20')) +
  geom_bar(stat = "identity") +
  labs(title = "Human: cell number", y = "", x = "Cell Number (scaled)") +
  scale_fill_manual(values = techcolor)+
  theme(axis.text.y = element_blank()) +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_vline(xintercept = log10(1e2+1)/max(log10(rowSums(t1)+1)), slope = 0, color = orange[1], linetype = 'dashed') +
  geom_vline(xintercept = log10(1e3+1)/max(log10(rowSums(t1)+1)), slope = 0, color = orange[2], linetype = 'dashed') +
  geom_vline(xintercept = log10(1e4+1)/max(log10(rowSums(t1)+1)), slope = 0, color = orange[3], linetype = 'dashed') +
  geom_vline(xintercept = log10(1e5+1)/max(log10(rowSums(t1)+1)), slope = 0, color = orange[4], linetype = 'dashed')

#####
data_df.mouse <- as.data.frame(as.table(t.mouse))
colnames(data_df.mouse) <- c("Celltype", "Tissue", "Proportion")
data_df.mouse$Celltype <-factor(as.character(data_df.mouse$Celltype),levels=overallorder) ###### !!


# data_df.mouse$Proportion <- data_df.mouse$Proportion * 100
p4 <- ggplot(data_df.mouse, aes(y = Celltype, x = Proportion, fill = Tissue)) +
  geom_bar(stat = "identity") +
  labs(title = "Mouse: tissue proportion", y = "", x = "Proportion") +
  #scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(values = barcolor, breaks = names(barcolor))+
  theme(axis.text.y = element_blank())

tech.mouse <- as.data.frame(as.table(t.mouse.tech))
colnames(tech.mouse) <- c("Celltype", "Technology", "Proportion")
# tech.mouse$Proportion <- tech.mouse$Proportion * 100
tech.mouse$Celltype <-factor(as.character(tech.mouse$Celltype),levels=overallorder) ###### !!

p5 <- ggplot(tech.mouse, aes(y = Celltype, x = Proportion, fill = Technology)) +
  geom_bar(stat = "identity") +
  labs(title = "Mouse: sequencing techonology proportion", y = "", x = "Proportion") +
  #scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(values = techcolor)+
  theme(axis.text.y = element_blank(), legend.position = 'none')

####
numCell.mouse2 = numCell.mouse
str(numCell.mouse2)
num_df = data.frame(Celltype = names(numCell.mouse2), 
                    Number = numCell.mouse2, 
                    stringsAsFactors = F)
num_df$Celltype <-factor(as.character(num_df$Celltype),levels=overallorder) ###### !!

# num_df[,1] = factor(num_df[,1],levels=names(numCell.mouse2))
p6 <- ggplot(num_df, 
       aes(y = Celltype, x = Number, fill = 'grey20')) +
  geom_bar(stat = "identity") +
  labs(title = "Mouse: cell number", y = "", x = "Cell Number (scaled)") +
  scale_fill_manual(values = techcolor)+
  theme(axis.text.y = element_blank())+
  geom_vline(xintercept = log10(1e2+1)/max(log10(rowSums(t2)+1)), slope = 0, color = orange[1], linetype = 'dashed') +
  geom_vline(xintercept = log10(1e3+1)/max(log10(rowSums(t2)+1)), slope = 0, color = orange[2], linetype = 'dashed') +
  geom_vline(xintercept = log10(1e4+1)/max(log10(rowSums(t2)+1)), slope = 0, color = orange[3], linetype = 'dashed') +
  geom_vline(xintercept = log10(1e5+1)/max(log10(rowSums(t2)+1)), slope = 0, color = orange[4], linetype = 'dashed')


library(gridExtra)
library(egg)


#grid.arrange(p1,p2,p3,p4,p5,p6, ncol = 1)
p <- egg::ggarrange(p1,p4, p2,p5,p3,p6, nrow = 1)
dev.off()
pdf('/Users/wenpinhou/Dropbox/encode4/circlize/barplot_all.pdf',width=11,height=10)
grid::grid.draw(p)
dev.off()

