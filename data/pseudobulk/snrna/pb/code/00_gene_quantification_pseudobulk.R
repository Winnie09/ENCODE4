gene_quantification_pseudobulk <-
  function(dataset_path,
           dataset,
           barcode,
           savefilename) {
    pb <- sapply(dataset, function(f) {
      print(f)
      m <- readRDS(paste0(dataset_path, f))
      tmp = rowSums(m[, colnames(m) %in% bc, drop = F])
    }, simplify = T)
    
    if (ncol(pb) > 1) {
      pb <- rowSums(pb)
    } else {
      pb2 <- as.vector(pb)
      names(pb2) <- rownames(pb)
      pb <- pb2
    }
    cpm = pb / (sum(pb) / 1e6)
    
    ## save file
    gene_id = sub('.*:', '', names(pb))
    gene_name = sub(':.*', '', names(pb))
    df = data.frame(
      gene_id = gene_id,
      gene_name = gene_name,
      count = pb,
      CPM = cpm,
      stringsAsFactors = FALSE
    )
    colnames(df) = paste0('#', colnames(df))
    write.table(
      df,
      file = savefilename,
      sep = '\t',
      row.names = F,
      col.names = T
    )
    return('File saved.')
  }

