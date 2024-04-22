# Peaks are in /home/zji4/data-hji7/zji/encode/data/snatac/caper/peak/peak
# Cell type info is in /home/zji4/data-hji7/zji/encode/data/snatac/splitbam/ct/ct.rds
# calculate and visualize length of peaks for each pseudobulk (sample&celltype)
# these pseudobulk has passed > 5 million reads criteria 
# 20230313

## calculate numPeak
setwd('/home/whou10/data/zji/encode/data/snatac/caper/peak/peak')
af <- list.files()
d <- data.table::fread(af[1])
str(d)
lenPeak <- sapply(af, function(f){
  tmp <- fread(f)
  tmp[,3] - tmp[,2]
}, USE.NAMES = F)
names(lenPeak) <- af
names(lenPeak) <- sapply(names(lenPeak), function(i) sub('.narrowPeak.gz', '', i), USE.NAMES = F)
saveRDS(lenPeak, '/home/whou10/scratch4/whou10/encode4/cCRE/samplectpb/lenPeak/length_of_peaks.rds')

## find pseudobulk's cell types
ct <- readRDS('/home/whou10/data/zji/encode/data/snatac/splitbam/ct/ct.rds')
rownames(ct) <- ct[,1]
ov <- intersect(ct[,1], names(lenPeak))

## plot numPeak vs. celltype for pseudobulks
pd <- data.frame(file = names(lenPeak), 
                 lenPeak = sapply(lenPeak, median), 
                 celltype = ct[names(lenPeak), 2])
tmp <- tapply(pd[,2], pd[,3], mean, na.rm = T)
pd[,3] <- factor(pd[,3], levels = names(sort(tmp)))
library(ggplot2)
library(RColorBrewer)
source('/home/whou10/scratch16/whou10/resource/ggplot_theme.R')
theme_set(.new_theme)

setwd('/home/whou10/scratch4/whou10/encode4/cCRE/samplectpb/lenPeak/')
pdf('length_of_peak_median_vs_celltype_boxplot.pdf', width = 10, height = 4)
ggplot(data = pd, aes(x = celltype, y = lenPeak, fill = celltype)) + 
  geom_boxplot(alpha = 0.3) +
  theme(legend.position = 'none') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()


## max
pd <- data.frame(file = names(lenPeak), 
                 lenPeak = sapply(lenPeak, max), 
                 celltype = ct[names(lenPeak), 2])
tmp <- tapply(pd[,2], pd[,3], mean, na.rm = T)
pd[,3] <- factor(pd[,3], levels = names(sort(tmp)))
pdf('length_of_peak_max_vs_celltype_boxplot.pdf', width = 10, height = 4)
ggplot(data = pd, aes(x = celltype, y = lenPeak, fill = celltype)) + 
  geom_boxplot(alpha = 0.3) +
  theme(legend.position = 'none') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

## min
pd <- data.frame(file = names(lenPeak), 
                 lenPeak = sapply(lenPeak, min), 
                 celltype = ct[names(lenPeak), 2])
tmp <- tapply(pd[,2], pd[,3], mean, na.rm = T)
pd[,3] <- factor(pd[,3], levels = names(sort(tmp)))
pdf('length_of_peak_min_vs_celltype_boxplot.pdf', width = 10, height = 4)
ggplot(data = pd, aes(x = celltype, y = lenPeak, fill = celltype)) + 
  geom_boxplot(alpha = 0.3) +
  theme(legend.position = 'none') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

