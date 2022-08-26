#' DimPlot.patient
#'
#' @description DimPlot with grey background in split panels, and color each panel independently.
#'  When there are too many catagories in group.by, colors can be too similar and therefore hard
#'  to see their location on embedding plot. When we can put them into some meaningful groups, e.g.
#'  Group patients by disease, we can split.by group and use a whole color scale for each panel.
#'
#' @param ... same as Seurat DimPlot function.
#' @param grey Add cells belong to other panels as grey at the back, if TRUE.
#' @examples
#' DimPlot.patient(srt, group.by = "Patient", split.by = "Disease")
#' @import patchwork wrap_plots
#' @import ggplot2
#' @export

# This function subset data by facet argument, then treat each as an individual plot for coloring
# Default using DimPlot.grey

DimPlot.patient <- function(object,
                            dims = c(1, 2),
                            pt.size = 0.2,
                            reduction = "umap",
                            group.by = NULL,
                            split.by = NULL,
                            shape.by = NULL,
                            shuffle = F,
                            ncol = NULL,
                            grey = T){
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


  plots <- lapply(
    X = unique(data[, split.by]), # each split as a separate plot
    FUN = function(x) {
      plot <- ggplot()
      if (grey){
        plot <- plot + geom_point(mapping = aes(x=get(dims[1]), y=get(dims[2])), data = data, col = "grey", size = pt.size)
      }
      plot <- plot + geom_point(mapping = aes(x=get(dims[1]), y=get(dims[2]), col=get(group.by)), data = data[data[,split.by] == x,], size = pt.size) +
        xlab(dims[1]) + ylab(dims[2]) + theme_classic() + labs(col=group.by)
    }
  )
  plots <- wrap_plots(plots, ncol = ncol)

  return(plots)
}

