library(data.table)
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(zellkonverter))

library(Matrix)
act <- readRDS('/home/whou10/scratch4/whou10/encode4/data/celltype/proc/final.rds')
act <- act[!is.na(act$rna_dataset),]
exp <- fread('/home/whou10/scratch4/whou10/encode4/data/celltype/proc/exp',data.table=F)
exp$platform <- exp[,'Library construction platform']
exp$platform[exp$platform==''] <- exp[exp$platform=='','Library construction method']
exp$platform <- sub("10X Genomics Chromium Controller",'10x',sub("Parse Single Cell Whole Transcriptome Kit",'Parse',exp$platform))
act$platform <- exp[match(act$rna_dataset,exp$Accession),'platform']
act$k <- paste0(act$species,':',act$tissue,':',act$platform)

for (sct in unique(act$k)) {
  tem <- readLines('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_no_downsample/code/model/template')
  print(sct)
  tem <- gsub('samplenamereserve',gsub(' ','_',sct),tem)
  af <- unique(act[act$k==sct,'rna_dataset'])
  tem <- sub('pathreadreserve1',paste0(paste0(af,'=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/',af,'.p")'),collapse ='\n'),tem)
  tem <- sub('pathreadreserve2',paste0(af,collapse = ','),tem)
  writeLines(tem,paste0('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_no_downsample/code/model/',gsub(' ','_',sct),'.py'))
}

