library(reshape2)
library(data.table)
anno <- readRDS('/home/zji4/data-hji7/zji/encode/data/celltype/proc/final.rds')
anno$ctn <- paste0(anno$rna_dataset,':',anno$rna_library,':',anno$rna_barcode)
af <- list.files('/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/tissue/text')
res <- do.call(rbind,sapply(af,function(sf) {
	print(sf)
	d <- fread(paste0('/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/tissue/text/',sf,'/figures/cell_participation.csv'),data.table=F)
	rownames(d) <- d[,1]
	d <- as.matrix(d[,-1])
	colnames(d) <- tolower(colnames(d))
	ct <- anno[match(rownames(d),anno$ctn),'celltype']
	score <- rowsum(d,ct)
	score <- score/as.vector(table(ct)[rownames(score)])
	score <- sapply(colnames(score),function(i) {
		paste0(names(head(sort(score[,i],decreasing = T),5)),collapse = ';')
	})
	
	d2 <- fread(paste0('/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/tissue/text/',sf,'/figures/gene_rank_weights.to_csv'),data.table=F,header=T)[,-1]
	d2 <- sapply(colnames(d2),function(i) paste0(sub(';.*','',d2[1:5,i]),collapse = ';'))
	names(d2) <- tolower(sub('model_','',names(d2)))
	ag <- list.files(paste0('/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/tissue/text/',sf,'/figures'),pattern = 'GSEA')
	d3 <- sapply(ag,function(sg) {
		dt <- fread(paste0('/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/tissue/text/',sf,'/figures/',sg),data.table=F)
		dt <- dt[order(dt[,'FDR q-val'],-dt$NES),]
		dt <- dt[dt[,'FDR q-val'] < 0.05,]
		paste0(dt[1:5,'Term'],collapse = ';')
	})
	names(d3) <- gsub('GSEA_|.csv','',names(d3))
	tp <- paste0('topic_',1:length(score))
	data.frame(tissue=sf,topic=tp,celltype=score[tp],gene=d2[tp],GSEA=d3[tp])
},simplify = F))

write.csv(res,file='/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/tissue/compile/res.csv',row.names=F)
