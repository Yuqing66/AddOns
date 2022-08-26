#' DimPlot with grey background in split panels
#'
#' @description In Seurat DimPlot function, split.by parameter allows cells to be splited into
#' separate panels according to a meta column. However, it's sometimes hard to compare across panels because
#' some clusters are not existed in some conditions, and therefore the look of embeddings in each panel can
#' look drastically different.
#'
#' To locate cells on the full embedding easier, this function plot cells in other panels in grey at the back,
#' providing a shape of the full embedding.
#'
#' @param ... same as Seurat DimPlot function.
#' @examples
#' DimPlot.grey(srt, group.by = "Disease", split.by = "Disease")
#' @import ggplot2
#' @export



DimPlot.grey <- function(
    object,
    dims = c(1, 2),
    pt.size = 0.2,
    reduction = "umap",
    group.by = NULL,
    split.by = NULL,
    shape.by = NULL,
    shuffle = F,
    ncol = NULL,
    back.col = "lightgrey") {
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
    geom_point(data = data[, dims], mapping = aes(x=get(dims[1]), y=get(dims[2])), col = back.col, size = 0.1)

  g <- g + geom_point(data = data, mapping = aes(x=get(dims[1]), y=get(dims[2]), col=get(group.by)), size = pt.size) +
    xlab(dims[1]) + ylab(dims[2]) + theme_classic() + labs(col=group.by)
  if (!is.null(split.by)){
    g <- g + facet_wrap(reformulate(split.by), ncol = ncol)
  }

  return(g)
}
