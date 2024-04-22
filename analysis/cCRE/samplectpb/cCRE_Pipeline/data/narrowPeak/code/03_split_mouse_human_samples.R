meta <- readRDS('/home/whou10/scratch4/whou10/encode4/cCRE/samplectpb/cCRE_Pipeline/data/meta/ct.rds')
str(meta)
ddir <- '/home/whou10/scratch4/whou10/encode4/cCRE/samplectpb/cCRE_Pipeline/data/narrowPeak/raw/'
rdir <- '/home/whou10/scratch4/whou10/encode4/cCRE/samplectpb/cCRE_Pipeline/data/narrowPeak/processed/'
dir.create(paste0(ddir, 'human/'))
dir.create(paste0(ddir, 'mouse/'))
dir.create(paste0(rdir, 'human/'))
dir.create(paste0(rdir, 'mouse/'))
af <- list.files(ddir, pattern = '.gz$')

## subset the meta of all pseudobulks available
meta <- meta[meta[, 1] %in% intersect(sub('\\..*','',af), meta[,1]), ]

## split human and mouse files
mo <- meta[meta[,'species'] == 'Mus musculus', ]
hm <- meta[meta[,'species'] == 'Homo sapiens', ]

for (i in hm[,1]){
  system(paste0('mv ', ddir, i, '.narrowPeak.gz ', ddir, 'human/'))
  system(paste0('mv ', rdir, i, '.bed ', rdir, 'human/'))
}

for (i in mo[,1]){
  system(paste0('mv ', ddir, i, '.narrowPeak.gz ', ddir, 'mouse/'))
  system(paste0('mv ', rdir, i, '.bed ', rdir, 'mouse/'))
}
