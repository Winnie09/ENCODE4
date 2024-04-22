## check results
m = readRDS('/home/whou10/scratch4/whou10/encode4/data/pseudobulk/snrna/pb/merge/level1_Homo_sapiens_count.rds')
i = 'level1-ENCSR398PIW-macrophage'
v = m[,paste0(i, '.rds')]
d = read.table(paste0('/home/whou10/scratch4/whou10/encode4/data/pseudobulk/snrna/pb/uploadfile/level1/', i, '.tsv'), 
                header = T)
pb = d[,3]
names(pb) = paste0(d[,2],':',d[,1])
str(pb)
str(v)
identical(v, pb)
setdiff(v,pb)
plot(v, pb[names(v)])
identical(names(v), names(pb))
summary(v-pb[names(v)])
id = which(v-pb[names(v)] != 0)
g= names(id[1])
v[g] - pb[g]


ddir = '/home/whou10/scratch4/whou10/encode4/data/pseudobulk/snrna/pb/merge/'
sink('/home/whou10/scratch4/whou10/encode4/data/pseudobulk/snrna/pb/summary/non_integer_info.txt')
for (f in list.files(ddir, pattern = '_count')){
  print(f)
  m = readRDS(paste0(ddir, f))
  a = (m%%1==0)
  print('%integers in each gene')
  print(summary(rowMeans(a)))
  print('%integers in each pseudobulk')
  print(summary(colMeans(a)))
}
sink()



m = readRDS('/home/whou10/scratch4/whou10/encode4/data/pseudobulk/snrna/pb/merge/level1_Homo_sapiens_log2CPM.rds')
identical(v, log2(cpm + 1))
str(cpm)
str(v)
summary(v)
summary(cpm)
summary(log2(cpm + 1))



str(m)
head(m[,1])





