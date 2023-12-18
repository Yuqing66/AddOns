#' VlnPlot.scc with group label below each condition.
#'
#' @description With more possible values in group.by column, colors of groups can become hard to tell apart.
#' This function is a variant of VlnPlot, using VlnPlot.scc format, with group labels at the bottom.
#'
#' @param ... same as VlnPlot.ssc
#' @examples
#' VlnPlot.xlabel(srt, "IFNG", group.by = "Patient", split.by = "Disease")
#' @import ggplot2
#' @export




VlnPlot.xlabel <- function (object, gene, group.by, group.order = NULL, split.order = NULL, title = "", assay = "RNA", slot = "data", log_scale = F,
                            colors = NULL, split.by = NULL, spread = NULL, jitter_pts = T,
                            plot_mean = T, size = 1, sig = 3, number_labels = T, text_sizes = c(15, 10, 7, 10, 7, 7, 2.5), alpha = 0.5, theme = "classic")
{

  df <- object@meta.data[, colnames(object@meta.data) %in% c(gene, group.by, split.by), drop = F]
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  geneExp <- FetchData(object = object, vars = gene, slot = slot)
  if (sum(geneExp) == 0){
    warning("No expression in data")
    return()
  }
  colnames(geneExp) <- "value"
  df <- cbind(df, geneExp)
  if (!is.null(group.order)) df[,group.by] <- factor(df[,group.by], levels = group.order)
  if (!is.null(split.order)){
    for (i in 1:length(split.order)){
      term <- names(split.order)[i]
      if (!is.null(split.order[[term]])) df[,term] <- factor(df[,term], levels = split.order[[term]])
    }
  }

  colnames(df) <- gsub("-", "", colnames(df))
  gene <- gsub("-", "", gene)
  if (any(!is.null(spread))) {
    others <- setdiff(unique(df[, spread[1]]), spread[2])
    ind <- which(df[, spread[1]] == spread[2])
    rmdf <- df[ind, ]
    df <- df[-ind, ]
    for (i in 1:length(others)) {
      rmdf[, spread[1]] <- others[i]
      df <- rbind(df, rmdf)
    }
  }
  if (log_scale == T) {
    df$plot <- log2(df$value + 1)
  }
  else {
    df$plot <- df$value
  }
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  cols <- gg_color_hue(length(unique(df[, group.by])))
  g <- ggplot(df)
  if (all(!is.null(colors))) {
    g <- g + scale_color_manual(values = c(colors))
    g <- g + scale_fill_manual(values = c(colors))
  }
  if (theme == "bw") {
    g <- g + theme_bw()
  }
  else {
    g <- g + theme_classic()
  }
  if (title == "")
    title <- gene
  g <- g + labs(title = title, y = gene)
  g <- g + theme(plot.title = element_text(size = text_sizes[1]),
                 axis.title = element_text(size = text_sizes[2]), axis.text = element_text(size = text_sizes[3]),
                 legend.title = element_text(size = text_sizes[4]), legend.text = element_text(size = text_sizes[5]))
  g <- g + theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
                 axis.title.x = element_blank(),
                 axis.ticks.x = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 1, angle = 90))
  if (jitter_pts == T)
    g <- g + geom_jitter(aes_string(x = group.by, y = "plot",
                                    col = group.by), width = 0.2, size = size)
  g <- g + geom_violin(aes_string(x = group.by, y = "plot",
                                  fill = group.by), col = "black", trim = T, scale = "width", alpha = alpha)
  if (number_labels == T) {
    g <- g + stat_summary(aes_string(x = group.by, y = "value"),
                          fun.data = function(x) {
                            return(c(y = -max(df$plot)/25, label = length(x)))
                          }, colour = "black", geom = "text", size = text_sizes[7])
    g <- g + stat_summary(aes_string(x = group.by, y = "value"),
                          fun.data = function(x) {
                            return(c(y = -max(df$plot)/10, label = round(mean(as.numeric(x > 0)), sig)))
                          }, colour = "black", geom = "text", size = text_sizes[7])
  }
  if (plot_mean == TRUE) {
    scale <- max(df$plot)/max(tapply(df$value, INDEX = as.list(df[, colnames(df) %in% c(group.by, split.by), drop = F]),
                                     FUN = mean), na.rm = T)
    g <- g + suppressWarnings(stat_summary(aes_string(x = group.by, y = "value"), fun.y = function(x) mean(x) * (scale * 0.5), colour = "black", geom = "point", size = 2))
    g <- g + scale_y_continuous(sec.axis = sec_axis(~./(scale * 0.5), name = "Mean Expression"))
  }
  if (length(split.by) == 1) {
    g <- g + facet_grid(facets = reformulate(split.by), scales = "free_x", space = "free_x")
  }
  else if (length(split.by) == 2) {
    g <- g + facet_grid(facets = reformulate(split.by[1], split.by[2]), scales = "free_x", space = "free_x")
  }
  else if (length(split.by) > 2) {
    stop("Parameter split.by needs to be a string with equal or less than two variables.")
  }
  if (!is.null(split.by))
    g <- g + theme(strip.text.x = element_text(size = text_sizes[6]))
  return(g)
}

