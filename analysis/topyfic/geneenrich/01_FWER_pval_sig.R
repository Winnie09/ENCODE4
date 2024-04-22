rm(list=ls())
ddir = '/home/whou10/scratch4/whou10/encode4/topyfic/tissue_hvg/res/'
at = list.files(ddir)

d = readRDS('/home/whou10/scratch4/whou10/encode4/data/celltype/proc/final.rds')
ac = paste0(d$rna_dataset, ':', d$rna_library, ':', d$rna_barcode)

str(d)
ctall <- c('macrophage', 'endothelial cell','fibroblast', 'T cell')
ct = ctall[4]
for (ct in ctall){
  print(ct)
  ct %in% d$celltype
  sort(unique(d$celltype))
  ts =  paste0(sub(' ','_', d[d$celltype == ct, 'species']), ':', sub(' ', '_', d[d$celltype == ct, 'generaltissue']))
  table(ts)
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
    FWER.p.val = FWER.p.val[tb$FWER.p.val < 0.05]
    names(FWER.p.val) = tb$Term[tb$FWER.p.val < 0.05]
    FWER.p.val
  })
  names(res) = runts
  
  allterm = unique(unlist(sapply(res, function(i){
    names(i)
  })))
  str(allterm)
  
  m <- sapply(runts, function(ts.sel){
    print(ts.sel)
    # ts.sel = sort(unique(ts))[1]
    cp = read.csv(paste0(ddir, ts.sel, '/figures/','cell_participation.csv'), row.names = 1)
    cm = colMeans(cp[rownames(cp) %in% c, ])
    tb = read.csv(paste0(ddir, ts.sel, '/figures/','GSEA_topic_', sub('Topic_','',names(sort(cm, decreasing = T)[1])), '.csv'))
    FWER.p.val = tb$FWER.p.val
    names(FWER.p.val) = tb$Term
    FWER.p.val[allterm]
  })
  rownames(m) = allterm
  str(m)
  m = -log10(m+10e-100)
  str(m)
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
  
  pdf(paste0(pdir, ct, '_-log10FWER_pval_sig.pdf'), width = 8 + 1*ncol(m), height = 4+0.1*nrow(m))
  pheatmap(m, cluster_rows = F, cluster_cols = F, scale = 'none', na_col = 'grey', annotation_col = anncol,
           breaks = myBreaks, 
           legend_breaks = round(c(min(as.vector(m), na.rm = T), seq(-log10(0.05), max(as.vector(m), na.rm = T), length.out=5)),2))
  dev.off()
}

