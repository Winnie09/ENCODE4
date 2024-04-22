ct <- readRDS('/home/zji4/data-hji7/zji/encode/data/snatac/splitbam/ct/ct.rds')
n <- readLines('/home/zji4/data-hji7/zji/encode/data/snatac/splitbam/numread/numread')
nc <- as.numeric(n[c(1:(length(n)/2))*2])
names(nc) <- sub('.bam','',n[c(1:(length(n)/2))*2-1])
nc <- nc[order(names(nc))]
ct <- ct[order(ct[,1]),]
identical(ct[,1],names(nc))
ct <- data.frame(ct,nc)
ct[,1] <- sub('_.*','',ct[,1])
colnames(ct) <- c('sample','celltype','tissue','readnumber')
ct <- ct[order(-ct$readnumber),]
write.csv(ct,'/home/zji4/data-hji7/zji/encode/data/snatac/splitbam/numread/numwithct.csv',row.names=T)

