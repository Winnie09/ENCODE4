library(Matrix)
for (sct in list.files('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_no_downsample/res/')) {
	tem <- readLines('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_no_downsample/code/figure/template')
	print(sct)
	tem <- gsub('samplenamereserve',sct,tem)
	writeLines(tem,paste0('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_no_downsample/code/figure/',sct,'.py'))
}



