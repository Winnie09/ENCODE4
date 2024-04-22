rm(list=ls())
ddir = '/home/whou10/scratch4/whou10/encode4/topyfic/tissue_hvg/res/'
at = list.files(ddir)

d = readRDS('/home/whou10/scratch4/whou10/encode4/data/celltype/proc/final.rds')
ac = paste0(d$rna_dataset, ':', d$rna_library, ':', d$rna_barcode)

str(d)
# ct = 'T cell'
# ct = 'macrophage'
# ct = 'basal epithelial cell'
ct = 'macrophage'

ctall <- c('macrophage', 'endothelial cell','fibroblast', 'T cell')
for (ct in ctall){
  print(ct)
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
    FWER.p.val = tb$FWER.p.val
    # FWER.p.val = FWER.p.val[tb$FWER.p.val < 0.05]
    names(FWER.p.val) = tb$Term
    FWER.p.val
  })
  names(res) = runts
  
  allterm = unique(unlist(sapply(res, function(i){
    names(i)
  })))
  m = matrix(NA, nrow = length(allterm), ncol = length(res))
  dimnames(m) = list(allterm, runts)
  
  for (i in runts){
    print(i)
    tmp = res[[i]]
    str(tmp)
    m[names(tmp), i] = tmp
  }
  
  m = -log10(m+10e-100)
  m = m[order(rowMeans(m, na.rm = T), decreasing = T), ]
  
  library(pheatmap)
  anncol = data.frame(row.names = runts, PCC = apply(m, 2, cor, m[,1], use = 'pairwise.complete.obs'))
  
  paletteLength <- 100
  # length(breaks) == length(paletteLength) + 1
  # use floor and ceiling to deal with even/odd length pallettelengths
  myBreaks <- c(seq(min(as.vector(m), na.rm = T), -log10(0.05), length.out=ceiling(paletteLength/2) + 2), 
                seq(-log10(0.05), max(as.vector(m), na.rm = T), length.out=floor(paletteLength/2)))
  myBreaks = unique(myBreaks)
  
  pdir = '/home/whou10/scratch4/whou10/encode4/topyfic/geneenrich/'
  
  pdf(paste0(pdir, ct, '_-log10FWER_pval.pdf'), width = 8 + 1*ncol(m), height = 4+0.1*nrow(m))
  pheatmap(m, cluster_rows = F, cluster_cols = F, scale = 'none', na_col = 'grey', annotation_col = anncol,
           breaks = myBreaks)
  dev.off()
  
  
}

  
