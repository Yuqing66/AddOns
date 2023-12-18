# In this new version, default environment colors are used

# In seurat DimPlot function, split.by argument allows to facet data points into different plots,
# but unlike SingallingSingleCell package, the cells belong to other plots are not shown instead of gray dots,
# making it hard to relate cells to full data map.
# This function is plotting split.by with grey dots for seurat objects

# had trouble with DefaultDimReduc and SingleDimPlot functions

DimPlot.grey <- function(
  object,
  dims = c(1, 2),
  pt.size = 0.2,
  reduction = "umap",
  group.by = NULL,
  split.by = NULL,
  shape.by = NULL,
  shuffle = F,
  ncol = NULL) {
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  data <- Embeddings(object = object[[reduction]])[, dims]
  data <- as.data.frame(x = data)
  dims <- paste0(Key(object = object[[reduction]]), dims)
  object[['ident']] <- Idents(object = object)
  orig.groups <- group.by
  group.by <- group.by %||% 'ident'
  data <- cbind(data, object[[group.by]][, , drop = FALSE])
  group.by <- colnames(x = data)[3:ncol(x = data)] # seurat supports multiple group.by, but not in this function
  for (group in group.by) {
    if (!is.factor(x = data[, group])) {
      data[, group] <- factor(x = data[, group])
    }
  }
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  if (!is.null(x = split.by)) {
    data[, split.by] <- object[[split.by, drop = TRUE]]
  }
  if (isTRUE(x = shuffle)) {
    set.seed(seed = seed)
    data <- data[sample(x = 1:nrow(x = data)), ]
  }
  g <- ggplot() +
    geom_point(data = data[, dims], mapping = aes(x=get(dims[1]), y=get(dims[2])), col = "lightgrey", size = 0.1)

  g <- g + geom_point(data = data, mapping = aes_string(x=dims[1], y=dims[2], fill=group.by),colour="black", shape = 21, stroke=0.2) +
    xlab(dims[1]) + ylab(dims[2]) + theme_classic() + labs(col=group.by)
  if (!is.null(split.by)){
    g <- g + facet_wrap(reformulate(split.by), ncol = ncol)
  }
  g <- g + scale_fill_manual(values=colors.celltype$V2[match(levels(data$subCellType.ificc.main),
                                                             colors.celltype$V1)])
  return(g)
}
