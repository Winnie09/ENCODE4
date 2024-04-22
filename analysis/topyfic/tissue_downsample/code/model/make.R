suppressMessages(library(SingleCellExperiment))
suppressMessages(library(zellkonverter))

library(Matrix)
act <- readRDS('/home/whou10/scratch4/whou10/encode4/data/celltype/proc/final.rds')
act <- act[!is.na(act$rna_dataset),]
act$k <- paste0(act$species,':',act$tissue)
selct <- c('Homo_sapiens:heart_right_ventricle.py',
'Homo_sapiens:adrenal_gland.py',
'Mus_musculus:heart.py',
'Mus_musculus:adrenal_gland.py',
'Homo_sapiens:left_colon.py')

for (sct in unique(act$k)) {
#for (sct in selct) {
	tem <- readLines('/home/whou10/scratch4/whou10/encode4/topyfic/tissue/code/model/template')
	print(sct)
	tem <- gsub('samplenamereserve',gsub(' ','_',sct),tem)
	af <- unique(act[act$k==sct,'rna_dataset'])
	tem <- sub('pathreadreserve1',paste0(paste0(af,'=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/',af,'.p")'),collapse ='\n'),tem)
	tem <- sub('pathreadreserve2',paste0(af,collapse = ','),tem)
	writeLines(tem,paste0('/home/whou10/scratch4/whou10/encode4/topyfic/tissue/code/model/',gsub(' ','_',sct),'.py'))
}






