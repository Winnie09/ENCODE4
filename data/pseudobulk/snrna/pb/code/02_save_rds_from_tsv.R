setwd('/home/whou10/scratch4/whou10/encode4/data/pseudobulk/snrna/pb/')
for (level in paste0('level', 1:3)){
  ddir <- paste0('uploadfile/', level, '/')
  rdir <- paste0(level, '/')
  i = 'level1-ENCSR480RWU-endothelial_cell.tsv'
  for (i in list.files(rdir)){
    tb = read.table(paste0(ddir, i), sep = '\t', header = T)
    v = tb[,3]
    names(v) = paste0(tb[,1], ':', tb[,2])
    saveRDS(v, paste0(rdir, sub('.tsv','.rds',i)))
  }
}
  