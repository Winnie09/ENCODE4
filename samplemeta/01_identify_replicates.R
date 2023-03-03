ddir <- '/home/whou10/data/zji/encode/data/snatac/sgr/rds/'
rdir <- '/home/whou10/scratch4/whou10/encode4/cCRE/cCRE_Pipeline/signal/res/'
af <- list.files(ddir)

## snatac-seq replicates
m <- readRDS('/home/whou10/data/zji/encode/data/snatac/filelist/summary.rds')
tab <- table(m[,'Experiment accession'])
encsr <- names(tab)[tab>1]
m2 <- m[m[, 'Experiment accession'] %in% encsr, ]
dim(m2)
write.csv(m2, '/home/whou10/scratch4/whou10/encode4/samplemeta/snatac.rep.csv')

## scrna-seq replicates
m <- readRDS('/home/whou10/data/zji/encode/data/scrna/filelist/summary.rds')
tab <- table(m[,'Experiment accession'])
encsr <- names(tab)[tab>1]
m3 <- m[m[, 'Experiment accession'] %in% encsr, ]
dim(m3)
write.csv(m3, '/home/whou10/scratch4/whou10/encode4/samplemeta/scrna.rep.csv')

