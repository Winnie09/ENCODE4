## =========================================
## pool pseudobulks using SEACells metacells
## =========================================
rm(list=ls())
## result directory
rdir.atac <- '/home/whou10/scratch4/whou10/encode4/SEACells/metacell/atac/processed/'  
rdir.rna <- '/home/whou10/scratch4/whou10/encode4/SEACells/metacell/rna/processed/'  
dir.create(rdir.atac)
dir.create(rdir.rna)
source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')

## Seacells metacell
af1 <- list.files('/home/whou10/data/zji/encode/data/multiome/seacells/csv/')

## ATAC count (peak by cell)
af2 <- list.files('/home/whou10/data/zji/encode/data/multiome/mat/rds/atac/')

## RNA count (gene by cell)
af3 <- list.files('/home/whou10/data/zji/encode/data/multiome/mat/rds/rna/')

## ====================================
## save metacell information as a list
## ====================================
f = af1[1]
f
metalist <- lapply(af1, function(f){
  print(f)
  meta = read.csv(paste0('/home/whou10/data/zji/encode/data/multiome/seacells/csv/', f))
  v <- meta[,2]
  names(v) <- meta[,1]
  v
})  
names(metalist) <- af1
saveRDS(metalist, '/home/whou10/scratch4/whou10/encode4/SEACells/metacell/metacell_list/metacell_list.rds')


## =======
## rna
## =======
f = af1[1]
f
for (f in af1){
  print(f)
  meta = read.csv(paste0('/home/whou10/data/zji/encode/data/multiome/seacells/csv/', f))
  rna = readRDS(paste0('/home/whou10/data/zji/encode/data/multiome/mat/rds/rna/', sub('.csv','.rds', f)))
  rs <- t(rowsum(t(rna), meta[match(colnames(rna), meta[,1]),2]))
  mat <- log2CPM_from_10x_count(rs)
  mat2 <- mat[rowMeans(mat > 0) > 0.1, ]
  saveRDS(mat2, paste0(rdir.rna, sub('.csv','_metacell.rds',f)))
}

## =========
## atac
## =========
f = af1[1]
f
## Because we  only have peaks there it’s only a small subset of everything,
## we need to use original libsize
## load atac encsr and encff table
libsize <- readRDS('/home/whou10/data/zji/encode/data/snatac/libsize/raw.rds')
str(libsize)
length(libsize)
str(libsize[[1]])

# Use this to map across ENCFF and ENCSR
tb <- readRDS('/home/whou10/data/zji/encode/data/samplemeta/snatac/summary.rds')
str(tb)
tb <- tb[tb$`File accession` %in% names(libsize), ]
encff <- tb$`File accession`
names(encff) <- tb$`Experiment accession`
encff <- encff[which(encff %in% names(libsize))]
str(encff)

## load encsr matching table: mltiome, rna, atac
tb.encsr <- readRDS('/home/whou10/data/zji/encode/data/samplemeta/multiomeacc/list.rds')
str(tb.encsr)
f = af1[1]

## mark mouse and human samples
species <- sapply(af1, function(f){
  encsr.f = tb.encsr[which(tb.encsr[,1] == sub(':.*', '', f)), 3] ## multiome encsr to atac encsr
  tb[tb[,3]==encsr.f,2]
})

## 
for (f in af1){
  print(f)
  meta = read.csv(paste0('/home/whou10/data/zji/encode/data/multiome/seacells/csv/', f))
  atac = readRDS(paste0('/home/whou10/data/zji/encode/data/multiome/mat/rds/atac/', sub('.csv','.rds', f)))
  rs <- t(rowsum(t(atac), meta[match(colnames(atac), meta[,1]),2]))
  rs2 <- (rs > 0) + 0
  saveRDS(rs2, paste0(rdir.atac, 'binary/',sub('.csv','_metacell.rds',f)))
  # ## normalize
  # encsr.f = tb.encsr[which(tb.encsr[,1] == sub(':.*', '', f)), 3] ## multiome encsr to atac encsr
  # lib <- libsize[encff[encsr.f]][[1]]
  # s <- tapply(lib[colnames(atac)], meta[match(colnames(atac), meta[,1]),2], sum)
  # log2(t(t(rs)/(lib[rownames(rs)] / factor) + 1))
  # 
  # ## filtering
  # mat2 <- mat[rowMeans(mat > 0) > 0.1, ]
  # saveRDS(mat2, paste0(rdir.atac, 'binary/',sub('.csv','_metacell.rds',f)))
}


## gene peak pairs
f = af1[1]
atac <- readRDS(paste0(rdir.atac, 'binary/',sub('.csv','_metacell.rds',f)))
rna <- readRDS(paste0(rdir.rna, sub('.csv','_metacell.rds',f)))
suppressMessages(library(bsseq))
### mouse
gtf <- data.table::fread('/home/whou10/scratch16/whou10/resource/gtf/grcm38.gtf',data.table=F)
gtf <- gtf[gtf[,3]=='gene',]
gene <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand=gtf[,7])
names(gene) <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
pro <- promoters(gene,upstream=500,downstream=500)
peaks <- GRanges(seqnames=sapply(rownames(atac), function(i) sub(':.*','',i)),
                 IRanges(start=as.numeric(sapply(rownames(atac), function(i) sub('_.*','',sub('.*:','',i)))),
                         end=as.numeric(sapply(rownames(atac), function(i) sub('.*_','',sub('.*:','',i))))),
                 strand='*')

o <- as.matrix(findOverlaps(pro, peaks))
str(o)
peakgenepair <- list()
for (i in o[,1]){
  peakgenepair[[names(gene)[i]]] <- rownames(atac)[o[o[,1]==i,2,drop=FALSE]]
}
length(peakgenepair)
saveRDS(peakgenepair, '/home/whou10/scratch4/whou10/encode4/SEACells/metacell/peakgenepair/mouse_grcm38.rds')

### human
gtf <- data.table::fread('/home/whou10/scratch16/whou10/resource/gtf/grch38.gtf',data.table=F)
gtf <- gtf[gtf[,3]=='gene',]
gene <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand=gtf[,7])
names(gene) <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
pro <- promoters(gene,upstream=500,downstream=500)
peaks <- GRanges(seqnames=sapply(rownames(atac), function(i) sub(':.*','',i)),
                 IRanges(start=as.numeric(sapply(rownames(atac), function(i) sub('_.*','',sub('.*:','',i)))),
                         end=as.numeric(sapply(rownames(atac), function(i) sub('.*_','',sub('.*:','',i))))),
                 strand='*')

o <- as.matrix(findOverlaps(pro, peaks))
peakgenepair <- list()
for (i in o[,1]){
  peakgenepair[[names(gene)[i]]] <- rownames(atac)[o[o[,1]==i,2,drop=FALSE]]
}
length(peakgenepair)
saveRDS(peakgenepair, '/home/whou10/scratch4/whou10/encode4/SEACells/metacell/peakgenepair/human_grch38.rds')



## calculate peak-gene pair correlation: metacell
rdir <- '/home/whou10/scratch4/whou10/encode4/SEACells/metacell/metacell_correlation/'
f = af1[2]

for (f in af1){
  print(f)
  if (species[f] == 'mm10'){
    peakgenepair <- readRDS('/home/whou10/scratch4/whou10/encode4/SEACells/metacell/peakgenepair/mouse_grcm38.rds')
    str(peakgenepair[1:3])
  } else {
    peakgenepair <- readRDS('/home/whou10/scratch4/whou10/encode4/SEACells/metacell/peakgenepair/human_grch38.rds')
    str(peakgenepair[1:3])
  }
  atac <- readRDS(paste0(rdir.atac, 'binary/',sub('.csv','_metacell.rds',f)))
  rna <- readRDS(paste0(rdir.rna, sub('.csv','_metacell.rds',f)))
  str(atac)
  str(rna)
  corlist <- list()
  for (g in rownames(rna)){
    if (sub(":.*", '', g) %in% names(peakgenepair)) {
      pk = peakgenepair[[sub(":.*", '', g)]]
      if (sum(pk %in% rownames(atac)) > 0) {
        atac.tmp =   t(atac[pk[pk %in% rownames(atac)], , drop = F])
        corlist[[g]] <-
          as.vector(corfunc(atac.tmp, t(rna[g, , drop = F])))
      }
    }
  }
  saveRDS(corlist, paste0(rdir, sub('.csv', '_metacell.rds', f)))    
}
str(corlist;1:3)  


## calculate single-cell level correlation
rdir <- '/home/whou10/scratch4/whou10/encode4/SEACells/metacell/singlecell_correlation/'
f = af1[1]
for (f in af1){
  print(f)
  if (species[f] == 'mm10'){
    peakgenepair <- readRDS('/home/whou10/scratch4/whou10/encode4/SEACells/metacell/peakgenepair/mouse_grcm38.rds')
  } else {
    peakgenepair <- readRDS('/home/whou10/scratch4/whou10/encode4/SEACells/metacell/peakgenepair/human_grch38.rds')
  }
  
  ## process sc data
  ## scRNA-seq log2CPM and gene QC
  rna = readRDS(paste0('/home/whou10/data/zji/encode/data/multiome/mat/rds/rna/', sub('.csv','.rds', f)))
  mat <- log2CPM_from_10x_count(rna)
  rna <- mat[rowMeans(mat > 0) > 0.1, ]
  
  ## scATAC-seq binarized
  atac = readRDS(paste0('/home/whou10/data/zji/encode/data/multiome/mat/rds/atac/', sub('.csv','.rds', f)))
  atac <- (atac > 0) + 0
  
  corlist <- list()
  for (g in rownames(rna)){
    if (sub(":.*",'', g) %in% names(peakgenepair)){
      pk = peakgenepair[[sub(":.*",'', g)]]
      if (sum(pk %in% rownames(atac)) > 0){
        atac.tmp =   t(atac[pk[pk %in% rownames(atac)], ,drop=F])
        corlist[[g]] <- as.vector(corfunc(atac.tmp, t(rna[g,,drop=F])))
      }
    }  
  }
  saveRDS(corlist, paste0(rdir, sub('.csv', '.rds', f)))
}



## =========
##### plot
## =========
library(ggplot2)
source('/home/whou10/scratch16/whou10/resource/ggplot_theme.R')
theme_set(.new_theme)
pdir <- '/home/whou10/scratch4/whou10/encode4/SEACells/metacell/metacell_list/'

# number of cells in each metacell
# all histogram
# stragify by tissue types (paste #of total cells)
tissue <- gsub('tissue.*', '', gsub('.*musculus ', '', gsub('.*sapiens ','', tb.encsr[,4])))
tb.encsr$tissue = tissue
ct <- gsub('.csv', '', sub('.*:','', names(metalist)))
table(tissue)
table(ct)
t = ct[1]
id <- which(ct == t)
i = id[1]
pd <- lapply(1:length(metalist), function(i){
  tmp <- table(metalist[[i]])
  
  tmpdf <- data.frame(species = species[i], 
             tissue = tb.encsr[which(tb.encsr[,1] %in% sub(':.*','',names(metalist)[i]) ), 5],
             celltype = ct[i],
             numCell = as.vector(tmp),
             stringsAsFactors = FALSE)
})
pd <- do.call(rbind, pd)
str(pd)
saveRDS(pd, paste0(pdir, 'plotdata.rds'))


library(RColorBrewer)

my_palette <- colorRampPalette(brewer.pal(n = 8, name =  "Dark2"))(20)

pdf(
  paste0(pdir, 'distribution_of_numCell_in_metacell.pdf'),
  width = 9,
  height = 3.5
)

ggplot(data = pd) +
  geom_boxplot(aes(x = celltype, y = numCell), 
               outlier.shape = NA) +
  geom_jitter(
    aes(
      x = celltype,
      y = numCell,
      color = tissue,
      shape = species
    ),
    size = 0.5,
    stroke = 0,
    alpha = 0.5,
    width =0.2
  ) +
  scale_y_log10() +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  )) +
  scale_alpha(guide = "none") +
  scale_size(guide = "none") +
  guides(
    size = guide_legend(override.aes = list(size = 6)),
    alpha = guide_legend(override.aes = list(alpha = 1))
  ) + 
  scale_color_manual(values = my_palette) +
  ylab('Number of cells in metacells (log10)') +
  xlab('Celltype')
  dev.off()

