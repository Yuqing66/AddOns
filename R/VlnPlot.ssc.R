#' VlnPlot.ssc Alternate violin plot function for Seurat object
#'
#' @description Left y-axis is for UMI level for individual cells, right y-axis is for the mean level showed in large black dots.
#'
#' @param object Seurat object to plot.
#' @param gene Gene name for plot. Wrap in ToMeta function if needed.
#' @param group.by Group and color cells.
#' @param split.by Split into panels. A vector with one or two elements. Take in order as c(row, column).
#' @param title Plot title.
#' @param assay Assay to extract gene UMI counts.
#' @param slot Slot to extract gene UMI counts.
#' @param log_scale Put y-axis into log scale if TRUE.
#' @param colors A vector of colors for scale_fill_manual / scale_color_manual to generate the gradient scale.
#' @param spread Control skin is only in HealthyControl in Disease meta column, but we want to compare them to other skins in all diseases
#' even when we split by Disease. So we spread HealthyControl condition in Disease column, by spread = c("Disease", "HealthyControl").
#' @param jitter_pts Overlay a jitter plot over violin plot if TRUE.
#' @param plot_mean Plot the mean UMI value in black dot if TRUE.
#' @param size Dot size of jitter plot.
#' @param sig Number of decimal places in fraction stats.
#' @param number_labels Show stats at the bottom of the plot if TRUE. First row is total cell number in each condition,
#' second row is fraction of cells with at least one UMI of this gene detected.
#' @param text_sizes Seven elements in order: plot title, axis title, axis text, legene title, legend text, split, stat size.
#' @param alpha Transparency of violin plot.
#' @param theme ggplot theme. classic by default.
#'
#' @examples
#' VlnPlot.ssc(srt, gene = "CXCL10", group.by = "Skin", split.by = c("CellType", "Disease"), spread = c("Disease", "HealthyControl"), assay = "RNA")
#' @import ggplot2
#' @export







VlnPlot.ssc <- function (object, gene, group.by, group.order = NULL, split.by = NULL, split.order = NULL, title = "", assay = "RNA", slot = "data", log_scale = F,
                         colors = NULL, spread = NULL, jitter_pts = T,
                         plot_mean = T, size = 1, sig = 3, number_labels = T, text_sizes = c(15, 10, 7, 10, 7, 7, 2.5), alpha = 0.5, theme = "classic")
{

  df <- object@meta.data[, colnames(object@meta.data) %in% c(gene, group.by, split.by), drop = F]
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  geneExp <- FetchData(object = object, vars = gene, slot = slot)
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
  g <- g + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5),
                 axis.title.x = element_blank(), axis.text.x = element_blank(),
                 axis.ticks.x = element_blank())
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
