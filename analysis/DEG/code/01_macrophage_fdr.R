## for a celltpye, identify DEG across tissue types
## one versus all
rm(list=ls())
library(Matrix)
ddir <- '/home/whou10/data/zji/encode/data/scrna/mat/mat/'
rdir = '/home/whou10/scratch4/whou10/encode4/DEG/res/'

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
  d[, colnames(d) %in% cell, drop = F]
})
pdf(paste0(rdir, 'macrophage_number_of_cells.pdf'), width = 5, height = 4)
hist(sapply(m, length), breaks = 100, col = 'grey', xlab = 'Number of cells in a sample', ylab = 'Number of samples',
     main = '')
dev.off()
tab = table(unlist(sapply(m, function(i) rownames(i))))
rn = names(tab[tab == length(m)])
m <- lapply(m, function(i) i[rn, ])
m = do.call(cbind, m)
m = log2(t(t(m)/(colSums(m) / 1e6) + 1))
saveRDS(m, paste0(rdir, 'macrophage_expr.rds'))

v = sub(':.*','',colnames(m))
ts = tmp[tmp$rna_dataset %in% v, 'tissue']
gts = tmp[tmp$rna_dataset %in% v, 'generaltissue']


## identify DEG
fdr <- parallel::mclapply(unique(gts), function(ts1){ ##
  print(ts1)
  pval <- sapply(1:nrow(m), function(i) { ##
    v1 = m[i, gts == ts1]
    v2 = m[i, gts != ts1]
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
}, mc.cores = 12)
saveRDS(fdr, paste0(rdir, 'fdr_macrophage_list.rds'))

fdr2 = do.call(cbind, fdr)
colnames(fdr2) = unique(gts)
saveRDS(fdr2, paste0(rdir, 'fdr_macrophage.rds'))