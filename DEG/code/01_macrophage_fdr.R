## for a celltpye, identify DEG across tissue types
## one versus all
rm(list=ls())
library(Matrix)
ddir <- '/home/whou10/data/zji/encode/data/scrna/mat/mat/'

allsp = readRDS('/home/whou10/scratch4/whou10/encode4/data/samplemeta/scrna/summary.rds')
allct = readRDS('/home/whou10/scratch4/whou10/encode4/data/celltype/proc/final.rds')
allct = allct[allct$species == 'Homo sapiens', ]

## prepare a matrix, columns are cells from all macrophages of all tissues
tmp = allct[allct$celltype=='macrophage', ]
str(tmp)
tmp = tmp[!is.na(tmp$rna_dataset),]
cell = paste0(tmp$rna_dataset, ':', tmp$rna_library, ':', tmp$rna_barcode)
af = list.files(ddir)
m <- lapply(unique(tmp$rna_dataset), function(sp){ ##
  print(sp)
  d = as.matrix(readRDS(paste0(ddir, sp,'.rds')))
  print(sum(colnames(d) %in% cell))
  d[, colnames(d) %in% cell]
})
tab = table(unlist(sapply(m, function(i) rownames(i))))
rn = names(tab[tab == length(m)])
m <- lapply(m, function(i) i[rn, ])
m = do.call(cbind, m)

v = sub(':.*','',colnames(m))
str(v)
ts = tmp[tmp$rna_dataset %in% v, 'tissue']
table(ts)
gts = tmp[tmp$rna_dataset %in% v, 'generaltissue']
str(ts)
str(m)


## identify DEG
fdr <- sapply(unique(ts), function(ts1){ ##
  pval <- sapply(1:nrow(m), function(i) { ##
    v1 = m[i, ts == ts1]
    v2 = m[i, ts != ts1]
    if (length(unique(c(v1,v2)))==1){
      1
    } else {
      wilcox.test(v1, v2)$p.value  
    }
  }, USE.NAMES = T)
  fdrv <- p.adjust(pval,method='fdr')
  names(fdrv) = rownames(m)
  summary(fdrv)
  fdrv
})
rdir = '/home/whou10/scratch4/whou10/encode4/DEG/res/'
saveRDS(fdr, paste0(rdir, 'fdr_macrophage.rds'))
