suppressMessages(library(GenomicAlignments))
ct <- readRDS('/home/zji4/data-hji7/zji/encode/data/celltype/proc/final.rds')

af <- list.files('/home/zji4/data-hji7/zji/encode/data/snatac/sgr/rds')
for (f in af) {
  d <- readRDS(paste0('/home/zji4/data-hji7/zji/encode/data/snatac/sgr/rds/',f))
  sct <- ct[which(ct$atac_dataset==sub('.rds','',f)),]
  names(d) <- sub('.*:','',names(d))
  if (sum(duplicated(names(d))) > 0) {
    stop(f)
  }
  for (k in unique(sct$celltype)) {
    if (mean(sct[sct$celltype==k,'atac_barcode'] %in% names(d)) < 1) {
      stop(f)
    }
    sd <- unlist(d[sct[sct$celltype==k,'atac_barcode']])
    saveRDS(sd,file=paste0('/home/zji4/data-hji7/zji/encode/data/snatac/pbgr/rds/',sub('.rds','',f),'-',gsub(' ','_',k),'.rds'))
  }
}
