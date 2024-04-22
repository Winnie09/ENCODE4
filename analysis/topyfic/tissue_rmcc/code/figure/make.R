library(Matrix)
for (sct in list.files('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_hvg/res/')) {
	tem <- readLines('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_hvg/code/figure/template')
	print(sct)
	tem <- gsub('samplenamereserve',sct,tem)
	writeLines(tem,paste0('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_hvg/code/figure/',sct,'.py'))
}




