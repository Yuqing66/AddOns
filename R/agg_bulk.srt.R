# seurat object as input. Return bulk matrix without putting in fData
agg_bulk.srt <- function (input, assay = "RNA", aggregate_by, group_by = FALSE, cutoff_frac = FALSE, 
                      cutoff_num = FALSE, cutoff_cpm = FALSE, rm.0CellsCondition = T) 
{
  input@active.assay <- assay
  if (group_by != FALSE) {
    ind <- match(group_by, aggregate_by)
    if (ind != 1) {
      stop("Please provide the group_by value first in the aggreggate_by argument")
    }
    groups <- sort(unique(input@meta.data[,group_by]))
  }
  to_expand <- vector("list", length(aggregate_by))
  for (i in 1:length(aggregate_by)) {
    var <- aggregate_by[i]
    vars <- unique(input@meta.data[,var])
    vars <- sort(vars)
    to_expand[[i]] <- vars
  }
  names(to_expand) <- aggregate_by
  bulks <- expand.grid(to_expand, stringsAsFactors = FALSE)
  colnames(bulks) <- c(aggregate_by)
  
  upm_vals <- c()
  num_cells_vals <- c()
  num_genes_vals <- c()
  rem_genes <- vector(mode = "list", length = nrow(bulks))
  
  for (j in 1:nrow(bulks)) {
    # get indices of cells meeting conditions in row j of bulks
    int <- bulks[j, ]
    full_match <- c()
    for (k in 1:length(int)) {
      ind <- which(input@meta.data[, colnames(bulks)[k]] == int[[k]])
      if (k == 1) {
        full_match <- c(full_match, ind)
      }
      else {
        full_match <- intersect(full_match, ind)
      }
    }
    num_cells <- length(full_match)
    num_cells_vals <- c(num_cells_vals, num_cells)
    if (length(full_match) > 1) {
      # Within all cells in a bulk column, sum each gene and divide by sum of all * 1million
      tmp <- input@assays$RNA@counts[, full_match]
      upm <- (apply(tmp, 1, sum)/sum(tmp)) * 1e+06
      if (cutoff_frac != FALSE || cutoff_num != FALSE) {
        zero_out_frac <- c()
        zero_out_num <- c()
        tmp2 <- tmp
        tmp2[which(tmp2 > 0)] <- 1
        gSums <- apply(tmp2, 1, sum)
        if (cutoff_frac != FALSE) {
          frac <- gSums/num_cells
          zero_out_frac <- which(frac < cutoff_frac)
        }
        if (cutoff_num != FALSE) {
          zero_out_num <- which(gSums < cutoff_num)
        }
        if (group_by != FALSE) {
          rem_genes[[j]] <- unique(names(c(zero_out_num, zero_out_frac)))
        }else {
          upm[unique(c(zero_out_num, zero_out_frac))] <- 0
        }
      }
    } else {
      # If not > 1 cells in this condition, assign 0 to every gene
      upm <- rep(0, nrow(input))
    }
    expressed <- length(which(upm > 0))
    upm_vals <- c(upm_vals, upm)
    num_genes_vals <- c(num_genes_vals, expressed)
  }
  names(rem_genes) <- apply(bulks, 1, FUN = paste, collapse = "_")
  bulks$numcells <- num_cells_vals
  bulks$numgenes <- num_genes_vals
  bulks$proportion <- 0
  if (group_by != FALSE) {
    for (i in 1:length(groups)) {
      ind <- grep(groups[i], bulks[, group_by])
      total <- sum(bulks$numcells[ind])
      bulks$proportion[ind] <- round((bulks$numcells[ind]/total) * 100, 2)
    }
  }else {
    # proportion of cells in each condition comparing to all
    total <- sum(bulks$numcells)
    bulks$proportion <- round((bulks$numcells/total) * 100, 2)
  }
  bulk <- matrix(upm_vals, nrow = nrow(input))
  rownames(bulk) <- rownames(input)
  colnames(bulk) <- seq(1:ncol(bulk))
  for (l in 1:nrow(bulks)) {
    cname <- bulks[l, -c(match(c("numcells", "proportion", 
                                 "numgenes"), colnames(bulks)))]
    cname2 <- c()
    for (i in 1:length(cname)) {
      cint <- as.character(cname[[i]])
      cname2 <- c(cname2, cint)
    }
    cname <- cname2
    cnum <- bulks[l, "numcells"]
    cpro <- bulks[l, "proportion"]
    cgen <- bulks[l, "numgenes"]
    cname <- paste0(c(cname, "num_genes", cgen, "num_cells", 
                      cnum, "percent", cpro, "bulk"), collapse = "_")
    colnames(bulk)[l] <- cname
  }
  if (group_by != FALSE) {
    if (!is.null(unlist(rem_genes))) {
      for (i in 1:length(vars)) {
        int_cell <- vars[i]
        ind <- grep(int_cell, names(rem_genes))
        vals <- table(unlist(rem_genes[ind]))
        zero_out <- names(which(vals == max(vals)))
        bulk[zero_out, ind] <- 0
      }
    }
  }
  if (cutoff_cpm) {
    for (i in 1:length(vars)) {
      int_cell <- vars[i]
      ind <- grep(int_cell, colnames(bulk))
      if (length(ind) > 1) {
        gCount <- apply(bulk[, ind], 1, function(x) length(which(x >= cutoff_cpm)))
        zero_out <- which(gCount == 0)
      }
      else {
        zero_out <- which(bulk[, ind] < cutoff_cpm)
      }
      bulk[zero_out, ind] <- 0
    }
  }
  # remove conditions that do not exist
  if (rm.0CellsCondition){
    rm <- which(apply(bulk, 2, sum) == 0)
    if (length(rm) != 0) bulk <- bulk[,-rm]
  }
  return(bulk)
}

# bulk_val <- agg_bulk.srt(srt, aggregate_by = c("Disease", "Skin", "CellType"))

