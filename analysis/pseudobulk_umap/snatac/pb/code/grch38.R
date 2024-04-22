suppressMessages(library(GenomicRanges))
gr <- readRDS("/home/zji4/data-hji7/zji/encode/data/snatac/pb/genomebin/grch38.rds")
af <- list.files('/home/zji4/data-hji7/zji/encode/data/snatac/pb/rds/')
s <- sub('-.*','',af)
ct <- readRDS('/home/zji4/data-hji7/zji/encode/data/celltype/proc/final.rds')
af <- af[s%in%unique(ct[ct[,1]=='Homo sapiens','atac_dataset'])]

m <- sapply(af,function(f) {
  print(f)
  d <- readRDS(paste0('/home/zji4/data-hji7/zji/encode/data/snatac/pb/rds/',f))
  log2(countOverlaps(gr,d)/length(d)*1e6 + 1)
})

saveRDS(m,file='/home/zji4/data-hji7/zji/encode/data/snatac/pb/mat/grch38.rds')

