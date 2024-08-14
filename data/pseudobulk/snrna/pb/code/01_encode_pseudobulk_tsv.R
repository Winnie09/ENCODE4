source('/home/whou10/scratch4/whou10/encode4/data/pseudobulk/snrna/pb/code/00_gene_quantification_pseudobulk.R')
## ========================================================
## load ct annotation file and prepare pseudobulk level123
## ========================================================
library(Matrix)
# d = readRDS('/home/whou10/scratch4/whou10/encode4/data/celltype/proc/old/final.rds')

d = readRDS('/home/whou10/scratch4/whou10/encode4/data/celltype/proc/step3.rds') ## 	1520779 obs. of  18 variables:

d <- d[!is.na(d$rna_dataset), ] # 	870041 obs. of  18 variables
d$lifestage <-
  ifelse(d$lifestage %in% c('adult', 'child'),
         'adult_child',
         'embryo_postnatal')
colnames(d) <- tolower(colnames(d))
l <-
  data.frame(
    barcode = paste0(d$rna_dataset, ':', d$rna_library, ':', d$rna_barcode),
    pseudobulk_level1 = gsub(' ', '_', paste0('level1-', d$rna_dataset, '-', d$celltype)),
    pseudobulk_level2 = gsub(
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
    pseudobulk_level3 = gsub(
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
    )
  )
str(l) ## 	870041 obs. of  4 variables:
saveRDS(l, '/home/whou10/scratch4/whou10/encode4/data/pseudobulk/snrna/pb/merge/levels.rds')


## =============================================
## generate pseudobulks for uploading to ENCODE
## =============================================
for (level in paste0('pseudobulk_level', 1:3)) {
  rdir <-
    paste0(
      '/home/whou10/scratch4/whou10/encode4/data/pseudobulk/snrna/pb/uploadfile/',
      sub('.*_','',level),
      '/'
    )
  alli = unique(l[,level])
  for (i in alli) {
    gene_quantification_pseudobulk(
      dataset_path = '/home/whou10/data/zji/encode/data/scrna/mat/mat/',
      dataset = paste0(unique(sapply(l[l[,level]==i, 1], function(i){sub(':.*','',i)})),'.rds'),
      barcode = l[l[,level]==i, 1],
      savefilename = paste0(rdir, i, '.tsv')
    )
  }
}

