suppressMessages(library(SingleCellExperiment))
suppressMessages(library(zellkonverter))

library(Matrix)
act <- readRDS('/home/whou10/scratch4/whou10/encode4/topyfic/final.rds')
act <- act[!is.na(act$rna_dataset),]
act$k <- paste0(act$species,':',act$tissue)
for (sct in unique(act$k)) {
	tem <- readLines('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_hvg/code/model/template')
	print(sct)
	tem <- gsub('samplenamereserve',gsub(' ','_',sct),tem)
	af <- unique(act[act$k==sct,'rna_dataset'])
	tem <- sub('pathreadreserve1',paste0(paste0(af,'=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/',af,'.p")'),collapse ='\n'),tem)
	tem <- sub('pathreadreserve2',paste0(af,collapse = ','),tem)
	writeLines(tem,paste0('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_hvg/code/model/',gsub(' ','_',sct),'.py'))
}

