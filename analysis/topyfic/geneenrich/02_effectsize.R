rm(list=ls())
ddir = '/home/whou10/scratch4/whou10/encode4/topyfic/tissue_hvg/res/'
at = list.files(ddir)

d = readRDS('/home/whou10/scratch4/whou10/encode4/data/celltype/proc/final.rds')
ac = paste0(d$rna_dataset, ':', d$rna_library, ':', d$rna_barcode)

ctall <- c('macrophage', 'endothelial cell','fibroblast', 'T cell')
for (ct in ctall){
  print(ct)
  ct %in% d$celltype
  ts =  paste0(sub(' ','_', d[d$celltype == ct, 'species']), ':', sub(' ', '_', d[d$celltype == ct, 'generaltissue']))
  c = ac[d$celltype == ct]
  names(c) = ts
  
  failts = c("Homo_sapiens:pancreas", "Mus_musculus:adrenal_gland", "Homo_sapiens:lung", "Mus_musculus:heart")
  runts = intersect(sort(unique(ts)), setdiff(at, failts))
  res <- lapply(runts, function(ts.sel){
    print(ts.sel)
    # ts.sel = sort(unique(ts))[1]
    cp = read.csv(paste0(ddir, ts.sel, '/figures/','cell_participation.csv'), row.names = 1)
    cm = colMeans(cp[rownames(cp) %in% c, ])
    tb = read.csv(paste0(ddir, ts.sel, '/figures/','GSEA_topic_', sub('Topic_','',names(sort(cm, decreasing = T)[1])), '.csv'))
    sum(tb$FWER.p.val < 0.05)
    es = tb$ES
    names(es) = tb$Term
    es
  })
  names(res) = runts
  
  
  allterm = unique(unlist(sapply(res, function(i){
    names(i)
  })))
  m = matrix(NA, nrow = length(allterm), ncol = length(res))
  dimnames(m) = list(allterm, runts)
  
  
  for (i in runts){
    tmp = res[[i]]
    m[names(tmp), i] = tmp
  }
  
  library(pheatmap)
  anncol = data.frame(row.names = runts, PCC = apply(m, 2, cor, m[,1], use = 'pairwise.complete.obs'))
  paletteLength <- 100
  # length(breaks) == length(paletteLength) + 1
  # use floor and ceiling to deal with even/odd length pallettelengths
  myBreaks <- c(seq(min(as.vector(m), na.rm = T), 0, length.out=ceiling(paletteLength/2) + 2), 
                seq(0, max(as.vector(m), na.rm = T), length.out=floor(paletteLength/2)))
  myBreaks = unique(myBreaks)
  
  pdir = '/home/whou10/scratch4/whou10/encode4/topyfic/geneenrich/'
  pdf(paste0(pdir, ct, '_effectsize.pdf'), width = 10, height = 15)
  pheatmap(m, cluster_rows = F, cluster_cols = F, scale = 'none', annotation_col = anncol, na_col = 'grey',
           show_rownames = F, breaks = myBreaks)
  dev.off()
  
}

  
