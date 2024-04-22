library(data.table)
suppressMessages(library(GenomicRanges))
ct <- readRDS('/home/whou10/scratch4/whou10/encode4/cCRE/samplectpb/cCRE_Pipeline/data/meta/ct.rds')
af <- sub('.narrowPeak.gz','',list.files('/home/whou10/data/zji/encode/data/snatac/caper/peak/peak'))

peak <- sapply(af,function(f) {
  k <- fread(paste0('/home/whou10/data/zji/encode/data/snatac/caper/peak/peak/',f,'.narrowPeak.gz'),data.table=F)
  GRanges(seqnames=k[,1],IRanges(k[,2],k[,3]))
},simplify = F)

for (spe in c('human','mouse')) {
  if (spe=='human') {
    taf <- intersect(af,ct[ct$species=='Homo sapiens',1])
  } else {
    taf <- intersect(af,ct[ct$species=='Mus musculus',1])
  }
  p <- fread(paste0('/home/whou10/scratch4/whou10/encode4/cCRE/samplectpb/cCRE_Pipeline/res/',spe,'/Final-rPeaks.bed'),data.table=F)
  p <- GRanges(seqnames=p[,1],IRanges(p[,2],p[,3]))
  o <- sapply(taf,function(i) {
    sum(countOverlaps(peak[[i]],p))
  })
  res <- data.frame(peaknumber=o,ct[match(names(o),ct[,1]),])
  res <- res[order(-res[,1]),]
  saveRDS(res,file=paste0('/scratch4/hji7/whou10/encode4/cCRE/samplectpb/cCRE_Pipeline/overlappeak/',spe,'.rds'))
}


