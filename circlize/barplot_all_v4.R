rm(list=ls())
options(scipen = 9999999)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(vegan) # for `reorder.hclust` (may be masked by the library `seriation`)
library(dendextend) # for `color_branches`
library(ggplot2)
# ad <- readRDS('/home/zji4/data-hji7/zji/encode/data/celltype/proc/final.rds')
ad <- readRDS('/Users/wenpinhou/Dropbox/encode4/data/celltype/proc/final.rds')
length(unique(ad$celltype)) ## 121


unique(ad[grep('neuron',ad$celltype), c('celltype', 'species', 'CL','tissue')]) #"CL:0002132" = stromal cell of ovary, Under stromal cell

# unique(ad[ad$CL=='CL:0000188', c('celltype', 'tissue')]) #"CL:0002132" = stromal cell of ovary, Under stromal cell

matchtb = read.csv('/Users/wenpinhou/Dropbox/encode4/circlize/circlize_plot_summarize_celltype.csv', sep = ',')
ctv = matchtb[,4]
for (i in 1:length(ctv)){
  print(i)
  text = paste0('(',length(which(matchtb[,4] == ctv[i] & matchtb[,2] == 'Mus musculus')),
                         ',',length(which(matchtb[,4] == ctv[i] & matchtb[,2] == 'Homo sapiens')),')')
  ad[ad$celltype == ctv[i], 'celltype']  <- matchtb[matchtb[,4] == ctv[i], 4] <- paste0(sub('\\(.*','',ctv[i]), text)
}


'alpha cell, gamma cell, epsilon cell'
'chromaffin cell'

matchtb[,1] = tolower(matchtb[,1])
matchtb[matchtb[,1] %in% c('hippocampus', 'layer of hippocampus', 'left cerebral cortex'),1] = 'brain'
for (i in 1:nrow(matchtb)){
  if (matchtb[i,2] == 'Mus musculus'){
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
    ad[which(ad$species == 'Homo sapiens' & ad$celltype == (matchtb[i,3])), 'celltype'] = matchtb[i,4]  
  }
}
length(unique(ad$celltype)) ## 80

sort(unique(ad$celltype)) 

c1 = unique(ad[ad$species == 'Mus musculus', 'celltype'])
c2 = unique(ad[ad$species == 'Homo sapiens', 'celltype'])
sort(setdiff(c1, c2))
  

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

hc.human=reorder(hclust(dist(t.human)),-as.matrix(t.human)%*%seq(ncol(t.human))^2)


# hc.mouse=reorder(hclust(dist(t.mouse)),-as.matrix(t.mouse)%*%seq(ncol(t.mouse))^2)

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

labelcolor=hcl(c(260,90,120,60,0,210,180,310)+15,60,70)
barcolor=colorRampPalette(readRDS('/Users/wenpinhou/Dropbox/resource/color74.rds'))(length(unique(ad$generaltissue[!is.na(ad$generaltissue)])))
names(barcolor) <- sort(unique(ad$generaltissue))
techcolor=RColorBrewer::brewer.pal(n = 6, 'Set3')[4:6]
names(techcolor)=unique(d$technology)
orange=RColorBrewer::brewer.pal(n = 5, 'Oranges')[2:5]


coldf = read.csv('/Users/wenpinhou/Dropbox/encode4/circlize/color_scheme.csv', row.names = 1)
colv = coldf[,1]
names(colv) = rownames(coldf)

# barcolor=readRDS('/Users/wenpinhou/Dropbox/resource/color74.rds')[1:length(unique(c(colnames(t.human), colnames(t.mouse))))]
labels=hc.human$labels[hc.human$order]

# Convert the matrix to a data frame for ggplot2
data_df.human <- as.data.frame(as.table(t.human[hc.human$order, ]))
colnames(data_df.human) <- c("Celltype", "Tissue", "Proportion")

data_df.human$Proportion <- data_df.human$Proportion * 100
p1 <- ggplot(data_df.human, aes(y = Celltype, x = Proportion, fill = Tissue)) +
  geom_bar(stat = "identity") +
  labs(title = "Human: tissue proportion", y = "", x = "Proportion") +
  #scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(values = barcolor)+
  theme_classic() +
  theme( legend.position = 'none')

tech.human <- as.data.frame(as.table(t.human.tech[hc.human$order, ]))
colnames(tech.human) <- c("Celltype", "Technology", "Proportion")
tech.human$Proportion <- tech.human$Proportion * 100
p2 <- ggplot(tech.human, aes(y = Celltype, x = Proportion, fill = Technology)) +
  geom_bar(stat = "identity") +
  labs(title = "Human: sequencing technology proportion", y = "Celltype", x = "Proportion") +
  #scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(values = techcolor)+
  theme_classic() +
  theme(axis.text.y = element_blank())

num_df.human = data.frame(Celltype = names(numCell.human), 
                    Number = numCell.human, 
                    stringsAsFactors = T)
levels(num_df.human[,1]) = names(numCell.human[hc.human$order])
p3 <- ggplot(num_df.human, 
       aes(y = Celltype, x = Number, fill = 'grey20')) +
  geom_bar(stat = "identity") +
  labs(title = "Human: cell number", y = "", x = "Number") +
  scale_fill_manual(values = techcolor)+
  theme_classic() +
  theme(axis.text.y = element_blank()) +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_vline(xintercept = log10(1e2+1)/max(log10(rowSums(t1)+1)), slope = 0, color = orange[1], linetype = 'dashed') +
  geom_vline(xintercept = log10(1e3+1)/max(log10(rowSums(t1)+1)), slope = 0, color = orange[2], linetype = 'dashed') +
  geom_vline(xintercept = log10(1e4+1)/max(log10(rowSums(t1)+1)), slope = 0, color = orange[3], linetype = 'dashed') +
  geom_vline(xintercept = log10(1e5+1)/max(log10(rowSums(t1)+1)), slope = 0, color = orange[4], linetype = 'dashed')

#####
data_df.mouse <- as.data.frame(as.table(t.mouse[hc.human$order, ]))
colnames(data_df.mouse) <- c("Celltype", "Tissue", "Proportion")
data_df.mouse$Proportion <- data_df.mouse$Proportion * 100
p4 <- ggplot(data_df.mouse, aes(y = Celltype, x = Proportion, fill = Tissue)) +
  geom_bar(stat = "identity") +
  labs(title = "Mouse: tissue proportion", y = "", x = "Proportion") +
  #scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(values = barcolor)+
  theme_classic() +
  theme(axis.text.y = element_blank())

tech.mouse <- as.data.frame(as.table(t.mouse.tech[hc.human$order, ]))
colnames(tech.mouse) <- c("Celltype", "Technology", "Proportion")
tech.mouse$Proportion <- tech.mouse$Proportion * 100
p5 <- ggplot(tech.mouse, aes(y = Celltype, x = Proportion, fill = Technology)) +
  geom_bar(stat = "identity") +
  labs(title = "Mouse: sequencing techonology proportion", y = "", x = "Proportion") +
  #scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(values = techcolor)+
  theme_classic() +
  theme(axis.text.y = element_blank(), legend.position = 'none')

####
numCell.mouse2 = numCell.mouse[hc.human$order]
str(numCell.mouse2)
num_df = data.frame(Celltype = names(numCell.mouse2), 
                    Number = numCell.mouse2, 
                    stringsAsFactors = F)
num_df[,1] = factor(num_df[,1],levels=names(numCell.mouse2))
p6 <- ggplot(num_df, 
       aes(y = Celltype, x = Number, fill = 'grey20')) +
  geom_bar(stat = "identity") +
  labs(title = "Mouse: cell number", y = "", x = "Number") +
  scale_fill_manual(values = techcolor)+
  theme_classic() +
  theme(axis.text.y = element_blank())+
  geom_vline(xintercept = log10(1e2+1)/max(log10(rowSums(t2)+1)), slope = 0, color = orange[1], linetype = 'dashed') +
  geom_vline(xintercept = log10(1e3+1)/max(log10(rowSums(t2)+1)), slope = 0, color = orange[2], linetype = 'dashed') +
  geom_vline(xintercept = log10(1e4+1)/max(log10(rowSums(t2)+1)), slope = 0, color = orange[3], linetype = 'dashed') +
  geom_vline(xintercept = log10(1e5+1)/max(log10(rowSums(t2)+1)), slope = 0, color = orange[4], linetype = 'dashed')


library(gridExtra)
library(egg)

pdf('/Users/wenpinhou/Dropbox/encode4/circlize/barplot_all.pdf',width=14,height=14)
#grid.arrange(p1,p2,p3,p4,p5,p6, ncol = 1)
p <- egg::ggarrange(p1,p4, p2,p5,p3,p6, nrow = 1)
grid::grid.draw(p)
dev.off()
dev.off()


