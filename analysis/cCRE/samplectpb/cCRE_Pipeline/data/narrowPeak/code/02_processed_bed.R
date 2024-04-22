rdir <-
  '/home/whou10/scratch4/whou10/encode4/cCRE/samplectpb/cCRE_Pipeline/data/narrowPeak/processed/'
ddir <-
  '/home/whou10/scratch4/whou10/encode4/cCRE/samplectpb/cCRE_Pipeline/data/narrowPeak/raw/'
af <- list.files(ddir, pattern = 'narrowPeak.gz')
f = af[1]
parallel::mclapply(af, function(f) {
  print(f)
  d <- data.table::fread(paste0(ddir, f))
  str(d)
  d[, 4] <-
    sapply(d[, 4], function(i)
      paste0(sub('_.*', '', f), '_', i))
  write.table(
    d[, 1:5],
    paste0(rdir, sub('.narrowPeak.gz', '.bed', f)),
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE,
    sep = '\t'
  )
  return(0)
}, mc.cores = 24)

