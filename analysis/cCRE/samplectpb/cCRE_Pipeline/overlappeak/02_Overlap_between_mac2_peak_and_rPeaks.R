## ==============================================
## overlapped peaks of macs2 peaks and rPeaks
## ==============================================
species <- 'human'
peakDir <- paste0('/home/whou10/scratch4/whou10/encode4/cCRE/cCRE_Pipeline/rPeak/', species,'/')
peakDir.s <- paste0('/home/whou10/scratch4/whou10/encode4/cCRE/cCRE_Pipeline/data/narrowPeak/processed/', species,'/')
rdir <- paste0('/home/whou10/scratch4/whou10/encode4/cCRE/cCRE_Pipeline/findOverlapPeak/res/', species, '/')
pdir <- paste0('/home/whou10/scratch4/whou10/encode4/cCRE/cCRE_Pipeline/findOverlapPeak/plot/', species, '/')

## read in peaks
suppressMessages(library(GenomicRanges))
peak <- read.table(paste0(peakDir, 'Final-rPeaks.bed'), sep = '\t')
n <- peak[,4]
high <- GRanges(seqnames=peak[,1],IRanges(start=round((peak[,2]+peak[,3])/2),end=round((peak[,2]+peak[,3])/2)))
names(high) <- n

## findOverlaps
af <- list.files(peakDir.s)
f <- af[1]

numPeak <- parallel::mclapply(af, function(f){
  peak.s <- read.table(paste0(peakDir.s, f), sep = '\t')
  #n <- peak.s[,4]
  peak.s <- GRanges(seqnames=peak.s[,1],IRanges(start=peak.s[,2],end=peak.s[,3]))
  #names(peak.s) <- n
  o <- as.matrix(findOverlaps(peak.s, high))
  length(unique(o[,1]))
}, mc.cores = 48)

numPeak <- unlist(numPeak)
names(numPeak) <- sub('.bed', '', af)
str(numPeak)

df <- data.frame(encsr = names(numPeak), numPeakOverlap = numPeak, pd.tissue[names(numPeak), 2:6])
saveRDS(df, paste0(rdir, 'numPeakOverlap.rds'))

df <- readRDS('/scratch4/hji7/whou10/encode4/cCRE/samplectpb/cCRE_Pipeline/overlappeak/human.rds')
str(df)

library(ggplot2)
library(RColorBrewer)
source('/home/whou10/scratch16/whou10/resource/ggplot_theme.R')
theme_set(.new_theme)

# numPeakOverlap vs. tissue
pdir <- '/scratch4/hji7/whou10/encode4/cCRE/samplectpb/cCRE_Pipeline/overlappeak/'
pdf(paste0(pdir, 'numPeakOverlap_human.pdf'), width = 10, height = 4)
ggplot(df, aes(x = tissue, y = peaknumber, color = tissue, alpha = 0.2, fill = tissue)) + 
  geom_violin(aes(alpha = 0.1), width=1, scale = 'width') +
  geom_boxplot(aes(alpha = 0.3), width=0.4, color="black", alpha=0.6) +
  geom_jitter(color = 'black', alpha = 0.5, size = 0.8, width = 0.3)+
  theme_classic()+
  theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylab('number of peaks overlapped with rPeak')+
  xlab('tissue')
dev.off()

df <- readRDS('/scratch4/hji7/whou10/encode4/cCRE/samplectpb/cCRE_Pipeline/overlappeak/mouse.rds')

# numPeakOverlap vs. tissue
pdf(paste0(pdir, 'numPeakOverlap_mouse.pdf'), width = 5, height = 5)
ggplot(df, aes(x = tissue, y = peaknumber, color = tissue, alpha = 0.2, fill = tissue)) + 
  geom_violin(aes(alpha = 0.1), width=1, scale = 'width') +
  geom_boxplot(aes(alpha = 0.3), width=0.4, color="black", alpha=0.6) +
  geom_jitter(color = 'black', alpha = 0.5, size = 0.8, width = 0.3)+
  theme_classic()+
  theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylab('number of peaks overlapped with rPeak')+
  xlab('tissue')
dev.off()



