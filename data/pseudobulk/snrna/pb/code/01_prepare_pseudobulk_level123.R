library(Matrix)
d = readRDS('/home/whou10/scratch4/whou10/encode4/data/celltype/proc/final.rds')
d <- d[!is.na(d$rna_dataset), ]
d$lifestage <-
  ifelse(d$lifestage %in% c('adult', 'child'),
         'adult_child',
         'embryo_postnatal')

l <-
  list(gsub(' ', '_', paste0('level1-', d$rna_dataset, '-', d$celltype)),
       gsub(
         ' ',
         '_',
         paste0(
           'level2-',
           d$species,
           '-',
           d$tissue,
           '-',
           d$celltype,
           '-',
           d$sex,
           '-',
           d$lifestage
         )
       ),
       gsub(
         ' ',
         '_',
         paste0(
           'level3-',
           d$species,
           '-',
           d$generaltissue,
           '-',
           d$celltype,
           '-',
           d$lifestage
         )
       ))
d2 = cbind(d, level1 = l[[1]], level2 = l[[2]], level3 = l[[3]])
for (level in paste0('level', 1:3)){
  rdir <- paste0('/home/whou10/scratch4/whou10/encode4/data/pseudobulk/snrna/pb/', level, '/')
  for (i in unique(d2[, level])){
    print(i)
    id = which(d2[, level]==i)
    bc = paste0(d2[id,'rna_dataset'], ':', d2[id, 'rna_library'], ':', d2[id,'rna_barcode'])
    pb <- sapply(unique(d[id,'rna_dataset']), function(f){
      print(f)
      m <- readRDS(paste0('/home/whou10/data/zji/encode/data/scrna/mat/mat/',f,'.rds'))
      tmp = rowSums(m[,colnames(m) %in% bc])
    }, simplify = T)
    
    if (ncol(pb) > 1) {
      pb <- rowSums(pb)
    } else {
      pb2 <- as.vector(pb)
      names(pb2) <- rownames(pb)
      pb <- pb2
    }
    saveRDS(pb, paste0(rdir, i, '.rds'))
  }
}
  
