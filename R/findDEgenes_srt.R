edgeRDE_srt <- function (input, groups, batch, sizefactor, lib_size, minCells, 
                         pVal = 1, contrast = list(c(-1, 1))) {  
  
  rvar <- apply(input@assays$RNA@counts, 1, var)
  idx.keep <- which(rvar > 0)
  input <- input[idx.keep, ]
  
  gfrac <- apply(input@assays$RNA@counts, 1, function(a) length(a[a > 
                                                                    0])/length(a))
  idx.keep <- which(gfrac > minCells)
  input <- input[idx.keep, ]
  
  y <- edgeR::DGEList(counts = input@assays$RNA@counts, group = groups)
  
  y$samples$lib.size <- lib_size
  
  if (is.null(sizefactor)) {
    message(sprintf("Calculating norm factors..."))
    flush.console()
    y <- edgeR::calcNormFactors(y)
  }
  else {
    y$samples$norm.factors = input@meta.data$sizefactor/lib_size
  }
  if (length(unique(batch)) > 1) {
    design <- model.matrix(~0 + batch + group, data = y$samples)
  }
  else {
    design <- model.matrix(~0 + group, data = y$samples)
  }
  message("Estimating dispersion...")
  flush.console()
  y <- edgeR::estimateDisp(y, design)
  message("Doing likelihood ratio fit...")
  flush.console()
  
  if (is.na(y$common.dispersion) == TRUE) {
    tab <- NA
  }
  if (is.na(y$common.dispersion) == FALSE){
    fit <- edgeR::glmFit(y, design)
    tab <- list()
    if (length(unique(batch)) > 1) {
      lrt <- edgeR::glmLRT(fit)
      tab[["contrast_1"]] <- as.data.frame(edgeR::topTags(lrt, 
                                                          p.value = pVal, n = Inf, sort.by = "logFC"))
    }
    else {
      for (j in 1:length(contrast)) {
        cStr <- sprintf("contrast_%d", j)
        lrt <- edgeR::glmLRT(fit, contrast = contrast[[j]])
        tab[[cStr]] <- as.data.frame(edgeR::topTags(lrt, 
                                                    p.value = pVal, n = Inf, sort.by = "logFC"))
      }
    }
    return(tab)
  }
}

findDEgenes_srt <- function (input, DEgroup, contrastID, sizefactor = NULL, 
                             lib_size = NULL, facet_by, minCells = 0.1, batchID = NULL, 
                             outdir = "DEresults/", outsuffix = "DEresults.tsv", pVal = 1, 
                             contrast = list(c(-1, 1))) {
  
  if (is.null(lib_size)) {
    libsize = colSums(input@assays$RNA@counts)
    input@meta.data$lib_size = libsize
    lib_size = "lib_size"
  }
  for (i in 1:length(unique(input@meta.data[, facet_by]))) {
    name = sort(unique(as.character(input@meta.data[, facet_by])))[i]
    print(paste0("Performing DE for ", name))
    
    idx = rownames(input@meta.data)[which(input@meta.data[, facet_by] == name)]
    cntr = rownames(input@meta.data)[which(input@meta.data[, facet_by] == name & input@meta.data[, DEgroup] == contrastID)]
    z = as.matrix(input@assays$RNA@counts)[, idx]
    if (!is.null(dim(z))) {
      if (ncol(z) > 2) {
        if (is.null(batchID)) {
          batch = as.factor(rep(1, ncol(z)))
        }
        else {
          batch_table = as.matrix(table(input@meta.data[idx, batchID], 
                                        input@meta.data[idx, DEgroup]))
          if (length(batch_table[batch_table == 0]) < 2) {
            batch = as.factor(input@meta.data[idx, batchID])
          }
          else {
            message(sprintf("Warning: not enough cells to include batch in model, setting batches to 1"))
            flush.console()
            batch = as.factor(rep(1, ncol(z)))
          }
        }
        groupList = rep(0, times = ncol(z))
        groupList[which(colnames(z) %in% cntr)] = 1
        if (length(unique(groupList)) > 1) {
          group = factor(groupList)
          
          Idents(input) <- input@meta.data[,facet_by]
          z <- subset(input, idents = name)
          z_libsize <- z@meta.data$lib_size
          
          tab = edgeRDE_srt(input = z, groups = group, batch = batch, 
                            sizefactor = sizefactor, lib_size = z_libsize, 
                            minCells = minCells, pVal = pVal, contrast = contrast)
          
          DEtablecntr = tab[["contrast_1"]]
          # outfile = paste(name, contrastID, outsuffix, sep = "_")
          outfile = paste(name, outsuffix, sep = "_")
          write.table(DEtablecntr, paste(outdir, outfile, 
                                         sep = ""), sep = "\t")
        }
      }
    }
  }
}
