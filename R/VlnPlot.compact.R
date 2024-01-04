#' compact violin plots to fit more information in less space
#'
#' @description Comparing to the regular VlnPlot.ssc function, a few elements are ignored by default, including cell count in condition, proportion with positive values, dots for individual cell value, and mean expression level.
#' Supports a few more options including horizontal violin, multiple genes,
#'
#' @param object Seurat object.
#' @param assay Assay of the seurat object for the expression values.
#' @param slot Slot of the seurat object for the expression values.
#' @param log_scale log2 scale the values or not. By default FALSE.
#' @param genes a vector of genes to plot.
#' @param group.by Group and color by this meta column.
#' @param split.by Split into different panels. A vector with one or two elements, in order of c(row, col)
#' @param spread Control skin is only in HealthyControl in Disease meta column, but we want to compare them to other skins in all diseases
#' even when we split by Disease. So we spread HealthyControl condition in Disease column, by spread = c("Disease", "HealthyControl").
#' @param split.dir The direction arranging split panels by one variable. "h" horizontally (default), "v" vertically.
#' @param split.scale Split panels share the same axis scale "fixed" (default), or use different scales "free".
#' @param split.label.pos Put the split panel label at the top and/or right as default, or choose "bottom" "left".
#' @param split.label.textonly Remove the rectangular frame from the panel labels. FALSE by default.
#' @param split.label.rotate Rotate the label for 90 degrees. "none" by default. "x" "y" "both" as options.
#' @param axis.hide Hide axis ticks and text for a simpler look. "none" by default. "x" "y" "both" as options.
#' @param legend.hide Hide legend. FALSE by default.
#' @param flip Change the violin plots to horizontal direction.
#' @param jitter_pts Dot for value of each individual cell. FALSE by default.
#' @param plot_mean Plot the mean expression level as a black dot. FALSE by default.
#' @param size Dot size of jitter plot.
#' @param sig Significant number of the cell proportion with positive values.
#' @param number_labels Show stats at the bottom of the plot if TRUE. Total cell number in that condition.
#' @param text_sizes Seven elements in order: plot title, axis title, axis text, legend title, legend text, split, stat size
#' @param plot.theme ggplot theme. classic by default.
#' @param title Plot title.
#' @param violin.linewidth Width of the violin borders.
#' @param gene.group.color Color by group of genes instead of meta data. A named list with gene names in each element
#' @param custom.colors A vector with RGB colors as values. If with names, colors will be assigned to specific groups.
#' @param violin.hide.cutoff When a condition has too few cells, the violin plot can be misleading. Hide violin plots when condition has cells fewer than this.
#' @examples
#' genes <- c("IRF8","CLEC9A","XCR1","CLEC10A","CD14","CD209","CXCL9","IL1B","CCL17","JCHAIN","CD207","CCR7")
#' VlnPlot.compact(srt, assay = "RNA", genes = genes, group.by = "subCellType.newmel", split.by = c("gene","subCellType.newmel"),
#' split.scale = "free", split.label.pos = c("left", "bottom"), split.label.textonly = T,
#' split.label.rotate = "both", axis.hide = "both", flip = T, legend.hide = T)
#'
#' VlnPlot.compact(srt, genes = "CXCL9", group.by = "subCellType.newmel", legend.hide = T) +
#'   theme(axis.text.x = element_text(angle = 90))
#' @import ggplot2 dplyr
#' @export




VlnPlot.compact <- function (object, genes, group.by, title = "", assay = NULL, slot = "data", log_scale = F,
                         split.by = NULL, spread = NULL, split.dir = "h", split.scale = "fixed", split.label.pos = c("top","right"),
                         split.label.textonly = F, split.label.rotate = "none",
                         axis.hide = "none", legend.hide = F, flip = F,
                         jitter_pts = F, plot_mean = F, size = 1, sig = 3, number_labels = F,
                         text_sizes = c(15, 0, 7, 10, 7, 7, 2.5), alpha = 0.5, theme = "classic", violin.linewidth = 1,
                         gene.group.color = NULL, custom.colors = NULL, violin.hide.cutoff = NULL)
{
  ### get data
  df <- object@meta.data[, colnames(object@meta.data) %in% c(genes, group.by, split.by), drop = F]
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  geneExp <- FetchData(object = object, vars = genes, slot = slot)
  df <- cbind(df, geneExp)
  if (length(genes) > 1){
    df <- df %>% pivot_longer(cols = which(colnames(df) %in% genes), names_to = "gene", values_to = "value")
    df$gene <- factor(df$gene, levels = genes)
  }else{
    colnames(df) <- sub(genes,"value",colnames(df))
  }

  colnames(df) <- gsub("-", "", colnames(df))
  genes <- gsub("-", "", genes)

  ### spread
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

  ### violin.hide.cutoff
  if (!is.null(violin.hide.cutoff)){
    df <- df %>% group_by_at(c(group.by, split.by)) %>%
      mutate(cellnum = length(value)) %>%
      filter(cellnum > violin.hide.cutoff)
  }

  ### scale
  if (log_scale == T) {
    df$plot <- log2(df$value + 1)
  }else {
    df$plot <- df$value
  }

  ### color by
  if (is.null(gene.group.color)){
    color.by <- group.by
  }else{
    df$genegroup <- NA
    for (i in 1:length(gene.group.color)){
      df$genegroup[df$gene %in% gene.group.color[[i]]] <- names(gene.group.color)[i]
      df$genegroup <- factor(df$genegroup, levels = names(gene.group.color))
    }
    color.by <- "genegroup"
  }


  ### plot
  g <- ggplot(df)
  ### customized colors
  if (!is.null(custom.colors)) {
    if (class(custom.colors) == "function"){
      custom.colors <- custom.colors(length(unique(df[,color.by])))
    }else if(!is.null(names(custom.colors)) & !is.null(levels(df[,color.by]))){
      custom.colors <- custom.colors[match(names(custom.colors), levels(df[,color.by]))]
    }
    g <- g + scale_color_manual(values = custom.colors)
    g <- g + scale_fill_manual(values = custom.colors)
  }
  ### theme
  if (theme == "bw") {
    g <- g + theme_bw()
  }else {
    g <- g + theme_classic()
  }
  g <- g + labs(title = title, y = genes)
  g <- g + theme(plot.title = element_text(size = text_sizes[1]),
                 axis.title = element_text(size = text_sizes[2]), axis.text = element_text(size = text_sizes[3]),
                 legend.title = element_text(size = text_sizes[4]), legend.text = element_text(size = text_sizes[5]))
  if (!is.null(split.by)){g <- g + theme(strip.text.x = element_text(size = text_sizes[6]))}
  ### jitter (points)
  if (jitter_pts == T) {
    g <- g + geom_jitter(aes_string(x = group.by, y = "plot",
                                    col = color.by), width = 0.2, size = size)
  }
  ### violin
  g <- g + geom_violin(aes_string(x = group.by, y = "plot", fill = color.by),
                       col = "black", trim = T, scale = "width", alpha = alpha, linewidth = violin.linewidth)
  ### stats
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
  ### split
  if (length(split.by) == 1) {
    if (length(split.label.pos) > 1){
      message("Please choose split.label.pos on top, bottom, left, or right")
    }
    if (split.dir == "h"){
      g <- g + facet_wrap(facets = reformulate(split.by), dir = split.dir, scales = split.scale, strip.position = split.label.pos, nrow = 1)
    }else if(split.dir == "v"){
      g <- g + facet_wrap(facets = reformulate(split.by), dir = split.dir, scales = split.scale, strip.position = split.label.pos, ncol = 1)
    }
  }else if (length(split.by) == 2) {

    if(all(c("bottom","left") %in% split.label.pos)){
      lab.switch <- "both"
    }else if ("left" %in% split.label.pos) {
      lab.switch <- "y"
    }else if("bottom" %in%split.label.pos){
      lab.switch <- "x"
    }else{
      lab.switch <- NULL
    }
    g <- g + facet_grid(facets = reformulate(split.by[1], split.by[2]), scales = split.scale, space = "free",
                        switch = lab.switch) +
      theme(strip.placement = "outside")
  }else if (length(split.by) > 2) {
    stop("Parameter split.by needs to be a string with equal or less than two variables.")
  }
  ### split appearance
  if (split.label.textonly){
    g <- g + theme(strip.background = element_blank())
  }

  if (axis.hide == "y"){
    g <- g + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }else if (axis.hide == "x"){
    g <- g + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }else if (axis.hide == "both"){
    g <- g + theme(axis.text = element_blank(), axis.ticks = element_blank())
  }

  if (split.label.rotate == "x"){
    g <- g + theme(strip.text.x.bottom = element_text(angle = 90))
  }else if (split.label.rotate == "y"){
    g <- g + theme(strip.text.y.left = element_text(angle = 0))
  }else if (split.label.rotate == "both"){
    g <- g + theme(strip.text.y.left = element_text(angle = 0), strip.text.x.bottom = element_text(angle = 90))
  }

  if (legend.hide){
    g <- g + theme(legend.position = "none")
  }

  if (flip){
    g <- g + coord_flip() +
      scale_x_discrete(limits=rev)
  }
  return(g)
}






