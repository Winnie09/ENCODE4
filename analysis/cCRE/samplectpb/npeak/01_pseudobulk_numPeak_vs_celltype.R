# Peaks are in /home/zji4/data-hji7/zji/encode/data/snatac/caper/peak/peak
# Cell type info is in /home/zji4/data-hji7/zji/encode/data/snatac/splitbam/ct/ct.rds
# count and visualize number of peaks for each pseudobulk (sample&celltype)
# these pseudobulk has passed > 5 million reads criteria 
# 20230313

## calculate numPeak
setwd('/home/whou10/data/zji/encode/data/snatac/caper/peak/peak')
af <- list.files()
d <- data.table::fread(af[1])
npeak <- sapply(af, function(f){
  tmp <- fread(f)
  nrow(tmp)
})
names(npeak) <- sapply(names(npeak), function(i) sub('.narrowPeak.gz', '', i), USE.NAMES = F)
saveRDS(npeak, '/home/whou10/scratch4/whou10/encode4/cCRE/samplectpb/npeak/number_of_peaks.rds')

## load numRead
v = readLines('/home/whou10/scratch4/whou10/encode4/cCRE/samplectpb/numRead/numread.txt')
nRead <- tmp <- NULL
for (i in seq(2, length(v), 2)){
  nRead = c(nRead, v[i])
}
nRead <- as.numeric(nRead)
for (i in seq(1, length(v), 2)){
  tmp = c(tmp, sub('.bam', '', v[i]))
}
names(nRead) <- tmp
str(nRead)


## find pseudobulk's cell types
ct <- readRDS('/home/whou10/data/zji/encode/data/snatac/splitbam/ct/ct.rds')
rownames(ct) <- ct[,1]
ov <- intersect(ct[,1], names(npeak))
str(ct)
## plot numPeak vs. celltype for pseudobulks
pd <- data.frame(file = names(npeak), 
                 numPeak = npeak, 
                 celltype = ct[names(npeak), 2],
                 numRead = nRead[ov])
str(pd)
tmp <- tapply(pd[,2], pd[,3], mean, na.rm = T)
pd[,3] <- factor(pd[,3], levels = names(sort(tmp)))
library(ggplot2)
library(RColorBrewer)
source('/home/whou10/scratch16/whou10/resource/ggplot_theme.R')
theme_set(.new_theme)

setwd('/home/whou10/scratch4/whou10/encode4/cCRE/samplectpb/numRead')
pdf('number_of_peak_vs_celltype_boxplot.pdf', width = 10, height = 4)
ggplot(data = pd, aes(x = celltype, y = numPeak, fill = celltype)) + 
  geom_boxplot(alpha = 0.3) +
  theme(legend.position = 'none') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

## numRead by ct
pdf('number_of_read_vs_celltype_boxplot.pdf', width = 10, height = 4)
ggplot(data = pd, aes(x = celltype, y = numRead, fill = celltype)) + 
  geom_boxplot(alpha = 0.3) +
  theme(legend.position = 'none') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()


## numRead vs. numPeak
pdf('number_of_read_vs_number_of_peak.pdf', width = 2, height = 2)
ggplot(data = pd, aes(x = numPeak, y = numRead, color = celltype)) + 
  geom_point(alpha = 0.5, size = 1, stroke = 0) +
  theme(legend.position = 'none') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  xlab('Number of peaks') + 
  ylab('Number of reads') +
  ggtitle('SampleCt_Pseudobulk')
dev.off() 

