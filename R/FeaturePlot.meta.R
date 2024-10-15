#' FeaturePlot.meta, FeaturePlot for genelist with metagenes in.
#'
#' @description Although ToMeta can easily transform any gene to its meta name, the saved FeaturePlot can be confusing not
#' knowing what was the gene name wanted to plot at first. This function embed ToMeta inside, and display the original gene name.
#'
#' @param ... same as Seurat FeaturePlot function.
#' @param data.ToMeta meta gene table passed to ToMeta function.
#' @param to.ToMeta to parameter in ToMeta function. The gene name format used in this Seurat object.
#' @examples
#' FeaturePlot.meta(srt, "FOXP3")
#' @importFrom cowplot theme_cowplot
#' @importFrom patchwork wrap_plots
#' @importFrom Seurat AutoPointSize SetQuantile SingleDimPlot
#' @import ggplot2
#' @export



FeaturePlot.meta <- function(
    object,
    features,
    data.ToMeta = NULL,
    to.ToMeta = "mergedname",
    dims = c(1, 2),
    cells = NULL,
    cols = if (blend) {
      c('lightgrey', '#ff0000', '#00ff00')
    } else {
      c('lightgrey', 'blue')
    },
    pt.size = NULL,
    order = FALSE,
    min.cutoff = NA,
    max.cutoff = NA,
    reduction = NULL,
    split.by = NULL,
    keep.scale = "feature",
    shape.by = NULL,
    slot = 'data',
    blend = FALSE,
    blend.threshold = 0.5,
    label = FALSE,
    label.size = 4,
    label.color = "black",
    repel = FALSE,
    ncol = NULL,
    coord.fixed = FALSE,
    by.col = TRUE,
    sort.cell = NULL,
    interactive = FALSE,
    combine = TRUE,
    raster = NULL,
    raster.dpi = c(512, 512)
) {
  # TODO: deprecate fully on 3.2.0
  if (!is.null(x = sort.cell)) {
    warning(
      "The sort.cell parameter is being deprecated. Please use the order ",
      "parameter instead for equivalent functionality.",
      call. = FALSE,
      immediate. = TRUE
    )
    if (isTRUE(x = sort.cell)) {
      order <- sort.cell
    }
  }
  if (interactive) {
    return(IFeaturePlot(
      object = object,
      feature = features[1],
      dims = dims,
      reduction = reduction,
      slot = slot
    ))
  }
  # Check keep.scale param for valid entries
  if (!(is.null(x = keep.scale)) && !(keep.scale %in% c("feature", "all"))) {
    stop("`keep.scale` must be set to either `feature`, `all`, or NULL")
  }
  # Set a theme to remove right-hand Y axis lines
  # Also sets right-hand Y axis text label formatting
  no.right <- theme(
    axis.line.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    axis.title.y.right = element_text(
      face = "bold",
      size = 14,
      margin = margin(r = 7)
    )
  )
  # Get the DimReduc to use
  reduction <- reduction %||% DefaultDimReduc(object = object)
  if (length(x = dims) != 2 || !is.numeric(x = dims)) {
    stop("'dims' must be a two-length integer vector")
  }
  # Figure out blending stuff
  if (blend && length(x = features) != 2) {
    stop("Blending feature plots only works with two features")
  }
  # Set color scheme for blended FeaturePlots
  if (blend) {
    default.colors <- eval(expr = formals(fun = FeaturePlot)$cols)
    cols <- switch(
      EXPR = as.character(x = length(x = cols)),
      '0' = {
        warning(
          "No colors provided, using default colors",
          call. = FALSE,
          immediate. = TRUE
        )
        default.colors
      },
      '1' = {
        warning(
          "Only one color provided, assuming specified is double-negative and augmenting with default colors",
          call. = FALSE,
          immediate. = TRUE
        )
        c(cols, default.colors[2:3])
      },
      '2' = {
        warning(
          "Only two colors provided, assuming specified are for features and agumenting with '",
          default.colors[1],
          "' for double-negatives",
          call. = FALSE,
          immediate. = TRUE
        )
        c(default.colors[1], cols)
      },
      '3' = cols,
      {
        warning(
          "More than three colors provided, using only first three",
          call. = FALSE,
          immediate. = TRUE
        )
        cols[1:3]
      }
    )
  }
  if (blend && length(x = cols) != 3) {
    stop("Blending feature plots only works with three colors; first one for negative cells")
  }
  # Name the reductions
  dims <- paste0(Key(object = object[[reduction]]), dims)
  cells <- cells %||% colnames(x = object)

  ###### ------ CHANGE ------ ######
  # Get plotting data
  features.meta <- ToMeta(features, data = data.ToMeta, to = to.ToMeta)
  data <- FetchData(
    object = object,
    vars = c(dims, 'ident', features.meta), # transfer genes to meta names
    cells = cells,
    slot = slot
  )
  # transfer meta names back to gene names
  features.selected <- colnames(data)[4:ncol(data)] # some input genes may not exist in dataset. Used found ones instead of query vector
  ind <- which(features.selected %in% data.ToMeta[,to.ToMeta]) # index of returned genes in meta data
  if (length(ind) > 0){
    features.selected[ind] <- features[match(features.selected[ind], features.meta)] # change gene name back after grepping values
  }
  colnames(data)[4:ncol(data)] <- features.selected
  ###### ------ END ------ ######

  # Check presence of features/dimensions
  if (ncol(x = data) < 4) {
    stop(
      "None of the requested features were found: ",
      paste(features, collapse = ', '),
      " in slot ",
      slot,
      call. = FALSE
    )
  } else if (!all(dims %in% colnames(x = data))) {
    stop("The dimensions requested were not found", call. = FALSE)
  }
  features <- colnames(x = data)[4:ncol(x = data)]
  # Determine cutoffs
  min.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = min(data[, feature]),
        no = cutoff
      ))
    },
    cutoff = min.cutoff,
    feature = features
  )
  max.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = max(data[, feature]),
        no = cutoff
      ))
    },
    cutoff = max.cutoff,
    feature = features
  )
  check.lengths <- unique(x = vapply(
    X = list(features, min.cutoff, max.cutoff),
    FUN = length,
    FUN.VALUE = numeric(length = 1)
  ))
  if (length(x = check.lengths) != 1) {
    stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  }
  brewer.gran <- ifelse(
    test = length(x = cols) == 1,
    yes = brewer.pal.info[cols, ]$maxcolors,
    no = length(x = cols)
  )
  # Apply cutoffs
  data[, 4:ncol(x = data)] <- sapply(
    X = 4:ncol(x = data),
    FUN = function(index) {
      data.feature <- as.vector(x = data[, index])
      min.use <- SetQuantile(cutoff = min.cutoff[index - 3], data.feature)
      max.use <- SetQuantile(cutoff = max.cutoff[index - 3], data.feature)
      data.feature[data.feature < min.use] <- min.use
      data.feature[data.feature > max.use] <- max.use
      if (brewer.gran == 2) {
        return(data.feature)
      }
      data.cut <- if (all(data.feature == 0)) {
        0
      }
      else {
        as.numeric(x = as.factor(x = cut(
          x = as.numeric(x = data.feature),
          breaks = brewer.gran
        )))
      }
      return(data.cut)
    }
  )
  colnames(x = data)[4:ncol(x = data)] <- features
  rownames(x = data) <- cells
  # Figure out splits (FeatureHeatmap)
  data$split <- if (is.null(x = split.by)) {
    RandomName()
  } else {
    switch(
      EXPR = split.by,
      ident = Idents(object = object)[cells, drop = TRUE],
      object[[split.by, drop = TRUE]][cells, drop = TRUE]
    )
  }
  if (!is.factor(x = data$split)) {
    data$split <- factor(x = data$split)
  }
  # Set shaping variable
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  # Make list of plots
  plots <- vector(
    mode = "list",
    length = ifelse(
      test = blend,
      yes = 4,
      no = length(x = features) * length(x = levels(x = data$split))
    )
  )
  # Apply common limits
  xlims <- c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[, dims[1]])))
  ylims <- c(floor(min(data[, dims[2]])), ceiling(x = max(data[, dims[2]])))
  # Set blended colors
  if (blend) {
    ncol <- 4
    color.matrix <- BlendMatrix(
      two.colors = cols[2:3],
      col.threshold = blend.threshold,
      negative.color = cols[1]
    )
    cols <- cols[2:3]
    colors <- list(
      color.matrix[, 1],
      color.matrix[1, ],
      as.vector(x = color.matrix)
    )
  }
  # Make the plots
  for (i in 1:length(x = levels(x = data$split))) {
    # Figure out which split we're working with
    ident <- levels(x = data$split)[i]
    data.plot <- data[as.character(x = data$split) == ident, , drop = FALSE]
    # Blend expression values
    if (blend) {
      features <- features[1:2]
      no.expression <- features[colMeans(x = data.plot[, features]) == 0]
      if (length(x = no.expression) != 0) {
        stop(
          "The following features have no value: ",
          paste(no.expression, collapse = ', '),
          call. = FALSE
        )
      }
      data.plot <- cbind(data.plot[, c(dims, 'ident')], BlendExpression(data = data.plot[, features[1:2]]))
      features <- colnames(x = data.plot)[4:ncol(x = data.plot)]
    }
    # Make per-feature plots
    for (j in 1:length(x = features)) {
      feature <- features[j]
      # Get blended colors
      if (blend) {
        cols.use <- as.numeric(x = as.character(x = data.plot[, feature])) + 1
        cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
      } else {
        cols.use <- NULL
      }
      data.single <- data.plot[, c(dims, 'ident', feature, shape.by)]
      # Make the plot
      plot <- SingleDimPlot(
        data = data.single,
        dims = dims,
        col.by = feature,
        order = order,
        pt.size = pt.size,
        cols = cols.use,
        shape.by = shape.by,
        label = FALSE,
        raster = raster,
        raster.dpi = raster.dpi
      ) +
        scale_x_continuous(limits = xlims) +
        scale_y_continuous(limits = ylims) +
        theme_cowplot() +
        CenterTitle()
      # theme(plot.title = element_text(hjust = 0.5))
      # Add labels
      if (label) {
        plot <- LabelClusters(
          plot = plot,
          id = 'ident',
          repel = repel,
          size = label.size,
          color = label.color
        )
      }
      # Make FeatureHeatmaps look nice(ish)
      if (length(x = levels(x = data$split)) > 1) {
        plot <- plot + theme(panel.border = element_rect(fill = NA, colour = 'black'))
        # Add title
        plot <- plot + if (i == 1) {
          labs(title = feature)
        } else {
          labs(title = NULL)
        }
        # Add second axis
        if (j == length(x = features) && !blend) {
          suppressMessages(
            expr = plot <- plot +
              scale_y_continuous(
                sec.axis = dup_axis(name = ident),
                limits = ylims
              ) +
              no.right
          )
        }
        # Remove left Y axis
        if (j != 1) {
          plot <- plot + theme(
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y.left = element_blank()
          )
        }
        # Remove bottom X axis
        if (i != length(x = levels(x = data$split))) {
          plot <- plot + theme(
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank()
          )
        }
      } else {
        plot <- plot + labs(title = feature)
      }
      # Add colors scale for normal FeaturePlots
      if (!blend) {
        plot <- plot + guides(color = NULL)
        cols.grad <- cols
        if (length(x = cols) == 1) {
          plot <- plot + scale_color_brewer(palette = cols)
        } else if (length(x = cols) > 1) {
          unique.feature.exp <- unique(data.plot[, feature])
          if (length(unique.feature.exp) == 1) {
            warning("All cells have the same value (", unique.feature.exp, ") of ", feature, ".")
            if (unique.feature.exp == 0) {
              cols.grad <- cols[1]
            } else{
              cols.grad <- cols
            }
          }
          plot <- suppressMessages(
            expr = plot + scale_color_gradientn(
              colors = cols.grad,
              guide = "colorbar"
            )
          )
        }
      }
      if (!(is.null(x = keep.scale)) && keep.scale == "feature" && !blend) {
        max.feature.value <- max(data[, feature])
        min.feature.value <- min(data[, feature])
        plot <- suppressMessages(plot & scale_color_gradientn(colors = cols, limits = c(min.feature.value, max.feature.value)))
      }
      # Add coord_fixed
      if (coord.fixed) {
        plot <- plot + coord_fixed()
      }
      # I'm not sure why, but sometimes the damn thing fails without this
      # Thanks ggplot2
      plot <- plot
      # Place the plot
      plots[[(length(x = features) * (i - 1)) + j]] <- plot
    }
  }
  # Add blended color key
  if (blend) {
    blend.legend <- BlendMap(color.matrix = color.matrix)
    for (ii in 1:length(x = levels(x = data$split))) {
      suppressMessages(expr = plots <- append(
        x = plots,
        values = list(
          blend.legend +
            scale_y_continuous(
              sec.axis = dup_axis(name = ifelse(
                test = length(x = levels(x = data$split)) > 1,
                yes = levels(x = data$split)[ii],
                no = ''
              )),
              expand = c(0, 0)
            ) +
            labs(
              x = features[1],
              y = features[2],
              title = if (ii == 1) {
                paste('Color threshold:', blend.threshold)
              } else {
                NULL
              }
            ) +
            no.right
        ),
        after = 4 * ii - 1
      ))
    }
  }
  # Remove NULL plots
  plots <- Filter(f = Negate(f = is.null), x = plots)
  # Combine the plots
  if (is.null(x = ncol)) {
    ncol <- 2
    if (length(x = features) == 1) {
      ncol <- 1
    }
    if (length(x = features) > 6) {
      ncol <- 3
    }
    if (length(x = features) > 9) {
      ncol <- 4
    }
  }
  ncol <- ifelse(
    test = is.null(x = split.by) || blend,
    yes = ncol,
    no = length(x = features)
  )
  legend <- if (blend) {
    'none'
  } else {
    split.by %iff% 'none'
  }
  # Transpose the FeatureHeatmap matrix (not applicable for blended FeaturePlots)
  if (combine) {
    if (by.col && !is.null(x = split.by) && !blend) {
      plots <- lapply(
        X = plots,
        FUN = function(x) {
          return(suppressMessages(
            expr = x +
              theme_cowplot() +
              ggtitle("") +
              scale_y_continuous(sec.axis = dup_axis(name = ""), limits = ylims) +
              no.right
          ))
        }
      )
      nsplits <- length(x = levels(x = data$split))
      idx <- 1
      for (i in (length(x = features) * (nsplits - 1) + 1):(length(x = features) * nsplits)) {
        plots[[i]] <- suppressMessages(
          expr = plots[[i]] +
            scale_y_continuous(
              sec.axis = dup_axis(name = features[[idx]]),
              limits = ylims
            ) +
            no.right
        )
        idx <- idx + 1
      }
      idx <- 1
      for (i in which(x = 1:length(x = plots) %% length(x = features) == 1)) {
        plots[[i]] <- plots[[i]] +
          ggtitle(levels(x = data$split)[[idx]]) +
          theme(plot.title = element_text(hjust = 0.5))
        idx <- idx + 1
      }
      idx <- 1
      if (length(x = features) == 1) {
        for (i in 1:length(x = plots)) {
          plots[[i]] <- plots[[i]] +
            ggtitle(levels(x = data$split)[[idx]]) +
            theme(plot.title = element_text(hjust = 0.5))
          idx <- idx + 1
        }
        ncol <- 1
        nrow <- nsplits
      } else {
        nrow <- split.by %iff% length(x = levels(x = data$split))
      }
      plots <- plots[c(do.call(
        what = rbind,
        args = split(x = 1:length(x = plots), f = ceiling(x = seq_along(along.with = 1:length(x = plots)) / length(x = features)))
      ))]
      # Set ncol to number of splits (nrow) and nrow to number of features (ncol)
      plots <- wrap_plots(plots, ncol = nrow, nrow = ncol)
      if (!is.null(x = legend) && legend == 'none') {
        plots <- plots & NoLegend()
      }
    } else {
      plots <- wrap_plots(plots, ncol = ncol, nrow = split.by %iff% length(x = levels(x = data$split)))
    }
    if (!is.null(x = legend) && legend == 'none') {
      plots <- plots & NoLegend()
    }
    if (!(is.null(x = keep.scale)) && keep.scale == "all" && !blend) {
      max.feature.value <- max(data[, features])
      min.feature.value <- min(data[, features])
      plots <- suppressMessages(plots & scale_color_gradientn(colors = cols, limits = c(min.feature.value, max.feature.value)))
    }
  }
  return(plots)
}




#
#
# SingleDimPlot <- function(
#     data,
#     dims,
#     col.by = NULL,
#     cols = NULL,
#     pt.size = NULL,
#     shape.by = NULL,
#     alpha.by = NULL,
#     order = NULL,
#     label = FALSE,
#     repel = FALSE,
#     label.size = 4,
#     cells.highlight = NULL,
#     cols.highlight = '#DE2D26',
#     sizes.highlight = 1,
#     na.value = 'grey50',
#     raster = NULL,
#     raster.dpi = NULL
# ) {
#   pt.size <- pt.size %||% AutoPointSize(data = data, raster = raster)
#   if ((nrow(x = data) > 1e5) & !isFALSE(raster)){
#     message("Rasterizing points since number of points exceeds 100,000.",
#             "\nTo disable this behavior set `raster=FALSE`")
#   }
#   raster <- raster %||% (nrow(x = data) > 1e5)
#   if (!is.null(x = raster.dpi)) {
#     if (!is.numeric(x = raster.dpi) || length(x = raster.dpi) != 2)
#       stop("'raster.dpi' must be a two-length numeric vector")
#   }
#   if (length(x = dims) != 2) {
#     stop("'dims' must be a two-length vector")
#   }
#   if (!is.data.frame(x = data)) {
#     data <- as.data.frame(x = data)
#   }
#   if (is.character(x = dims) && !all(dims %in% colnames(x = data))) {
#     stop("Cannot find dimensions to plot in data")
#   } else if (is.numeric(x = dims)) {
#     dims <- colnames(x = data)[dims]
#   }
#   if (!is.null(x = cells.highlight)) {
#     highlight.info <- SetHighlight(
#       cells.highlight = cells.highlight,
#       cells.all = rownames(x = data),
#       sizes.highlight = sizes.highlight %||% pt.size,
#       cols.highlight = cols.highlight,
#       col.base = cols[1] %||% '#C3C3C3',
#       pt.size = pt.size
#     )
#     order <- highlight.info$plot.order
#     data$highlight <- highlight.info$highlight
#     col.by <- 'highlight'
#     pt.size <- highlight.info$size
#     cols <- highlight.info$color
#   }
#   if (!is.null(x = order) && !is.null(x = col.by)) {
#     if (typeof(x = order) == "logical") {
#       if (order) {
#         data <- data[order(!is.na(x = data[, col.by]), data[, col.by]), ]
#       }
#     } else {
#       order <- rev(x = c(
#         order,
#         setdiff(x = unique(x = data[, col.by]), y = order)
#       ))
#       data[, col.by] <- factor(x = data[, col.by], levels = order)
#       new.order <- order(x = data[, col.by])
#       data <- data[new.order, ]
#       if (length(x = pt.size) == length(x = new.order)) {
#         pt.size <- pt.size[new.order]
#       }
#     }
#   }
#   if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
#     warning("Cannot find ", col.by, " in plotting data, not coloring plot")
#     col.by <- NULL
#   } else {
#     # col.index <- grep(pattern = col.by, x = colnames(x = data), fixed = TRUE)
#     col.index <- match(x = col.by, table = colnames(x = data))
#     if (grepl(pattern = '^\\d', x = col.by)) {
#       # Do something for numbers
#       col.by <- paste0('x', col.by)
#     } else if (grepl(pattern = '-', x = col.by)) {
#       # Do something for dashes
#       col.by <- gsub(pattern = '-', replacement = '.', x = col.by)
#     }
#     colnames(x = data)[col.index] <- col.by
#   }
#   if (!is.null(x = shape.by) && !shape.by %in% colnames(x = data)) {
#     warning("Cannot find ", shape.by, " in plotting data, not shaping plot")
#   }
#   if (!is.null(x = alpha.by) && !alpha.by %in% colnames(x = data)) {
#     warning(
#       "Cannot find alpha variable ",
#       alpha.by,
#       " in data, setting to NULL",
#       call. = FALSE,
#       immediate. = TRUE
#     )
#     alpha.by <- NULL
#   }
#
#   plot <- ggplot(data = data)
#   plot <- if (isTRUE(x = raster)) {
#     plot + geom_scattermore(
#       mapping = aes_string(
#         x = dims[1],
#         y = dims[2],
#         color = paste0("`", col.by, "`"),
#         shape = shape.by,
#         alpha = alpha.by
#       ),
#       pointsize = pt.size,
#       pixels = raster.dpi
#     )
#   } else {
#     plot + geom_point(
#       mapping = aes_string(
#         x = dims[1],
#         y = dims[2],
#         color = paste0("`", col.by, "`"),
#         shape = shape.by,
#         alpha = alpha.by
#       ),
#       size = pt.size
#     )
#   }
#   plot <- plot +
#     guides(color = guide_legend(override.aes = list(size = 3))) +
#     labs(color = NULL, title = col.by) +
#     CenterTitle()
#   if (label && !is.null(x = col.by)) {
#     plot <- LabelClusters(
#       plot = plot,
#       id = col.by,
#       repel = repel,
#       size = label.size
#     )
#   }
#   if (!is.null(x = cols)) {
#     if (length(x = cols) == 1 && (is.numeric(x = cols) || cols %in% rownames(x = brewer.pal.info))) {
#       scale <- scale_color_brewer(palette = cols, na.value = na.value)
#     } else if (length(x = cols) == 1 && (cols %in% c('alphabet', 'alphabet2', 'glasbey', 'polychrome', 'stepped'))) {
#       colors <- DiscretePalette(length(unique(data[[col.by]])), palette = cols)
#       scale <- scale_color_manual(values = colors, na.value = na.value)
#     } else {
#       scale <- scale_color_manual(values = cols, na.value = na.value)
#     }
#     plot <- plot + scale
#   }
#   plot <- plot + theme_cowplot()
#   return(plot)
# }
#
#
#
#
#
#
#
DefaultDimReduc <- function(object, assay = NULL) {
  object <- UpdateSlots(object = object)
  assay <- assay %||% DefaultAssay(object = object)
  drs.use <- c('umap', 'tsne', 'pca')
  dim.reducs <- FilterObjects(object = object, classes.keep = 'DimReduc')
  drs.assay <- Filter(
    f = function(x) {
      return(DefaultAssay(object = object[[x]]) == assay)
    },
    x = dim.reducs
  )
  if (length(x = drs.assay) > 0) {
    index <- lapply(
      X = drs.use,
      FUN = grep,
      x = drs.assay,
      ignore.case = TRUE
    )
    index <- Filter(f = length, x = index)
    if (length(x = index) > 0) {
      return(drs.assay[min(index[[1]])])
    }
  }
  index <- lapply(
    X = drs.use,
    FUN = grep,
    x = dim.reducs,
    ignore.case = TRUE
  )
  index <- Filter(f = length, x = index)
  if (length(x = index) < 1) {
    stop(
      "Unable to find a DimReduc matching one of '",
      paste(drs.use[1:(length(x = drs.use) - 1)], collapse = "', '"),
      "', or '",
      drs.use[length(x = drs.use)],
      "', please specify a dimensional reduction to use",
      call. = FALSE
    )
  }
  return(dim.reducs[min(index[[1]])])
}




UpdateSlots <- function(object) {
  object.list <- sapply(
    X = slotNames(x = object),
    FUN = function(x) {
      return(tryCatch(
        expr = slot(object = object, name = x),
        error = function(...) {
          return(NULL)
        }
      ))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  object.list <- Filter(f = Negate(f = is.null), x = object.list)
  object.list <- c('Class' = class(x = object)[1], object.list)
  object <- do.call(what = 'new', args = object.list)
  for (x in setdiff(x = slotNames(x = object), y = names(x = object.list))) {
    xobj <- slot(object = object, name = x)
    if (is.vector(x = xobj) && !is.list(x = xobj) && length(x = xobj) == 0) {
      slot(object = object, name = x) <- vector(
        mode = class(x = xobj),
        length = 1L
      )
    }
  }
  return(object)
}
#
#
# SetQuantile <- function(cutoff, data) {
#   if (grepl(pattern = '^q[0-9]{1,2}$', x = as.character(x = cutoff), perl = TRUE)) {
#     this.quantile <- as.numeric(x = sub(
#       pattern = 'q',
#       replacement = '',
#       x = as.character(x = cutoff)
#     )) / 100
#     data <- unlist(x = data)
#     data <- data[data > 0]
#     cutoff <- quantile(x = data, probs = this.quantile)
#   }
#   return(as.numeric(x = cutoff))
# }
#
#
#
# AutoPointSize <- function(data, raster = NULL) {
#   return(ifelse(
#     test = isTRUE(x = raster),
#     yes = 1,
#     no = min(1583 / nrow(x = data), 1)
#   ))
# }
#
#
