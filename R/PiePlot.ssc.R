#' percent plot with denominator and divider specified.
#'
#' @description
#' This function uses a pie plot to show the percentage of each group in the total. The denominator and divider can be specified to normalize the percentage.
#'
#' @param object Seurat object. Full dataset for normalization denominators.
#' @param group.by Color by this meta column.
#' @param split.by Split into different panels.
#' @param layer.by Multiple layers of rings in each pie chart.
#' @param normalize.by A vector with metadata column names for normalizing groups. Use proportion as pie chart value. For example, some diseases have more cells in total, normalized by disease and pie chart for prevalence in each disease. NULL to plot the raw value.
#' @param ind.cell Indices for cells to show in plot. NULL by default, plotting all conditions.
#' @param width Width of each ring. 1 for no gap between rings.
#' @param title Plot title.
#' @param number_labels Show raw value or percentage in each condition. FALSE, "value", "percentage"
#' @param label.digit Number of decimal places to keep.
#' @param plot.theme ggplot theme. void by default. NULL to skip.
#' @examples
#' diseases <- c("HealthyControl", "Dermatomyositis", "Lupus", "Psoriasis", "Vitiligo-active")
#' skins <- c("Control", "Lesional")
#' celltypes <- c("MEL","KC_act")
#' celltypeColName <- "subCellType.ificc"
#'
#' object <- srt[,srt$Skin %in% skins &
#'                 srt$DiseaseFiner %in% diseases]
#'
#' ind <- object@meta.data[,celltypeColName] %in% celltypes
#' PiePlot.ssc(srt, group.by="DiseaseFiner", split.by="subCellType.ificc", normalize.by=c("DiseaseFiner"), ind.cell=ind,
#' width=0.9, title=NULL, label=T, color=NULL)
#' @import ggplot2 magrittr tibble dplyr tidyr circlize
#' @export


PiePlot.ssc <- function(object, group.by, split.by=NULL, layer.by=NULL, normalize.by=NULL, ind.cell=NULL,
                        width=0.9, title=NULL, number_labels="percentage", label.digit=4, plot.theme="void", color=NULL){
  tb <- object@meta.data[, c(group.by, split.by, layer.by), drop = F]

  if (!is.null(ind.cell)){
    tb$plot <- "n"
    tb$plot[ind.cell] <- "y"
  }else{
    tb$plot <- "y"
  }

  if (is.null(normalize.by)){
    mt <- tb %>%
      mutate(splitby = get(split.by)) %>%
      unite(numerator, c(group.by, split.by, layer.by), sep=":", remove = F) %>%
      group_by(numerator) %>%
      summarise(count = n(), display = unique(plot), splitby=unique(splitby))
  }else{
    mt <- tb %>%
      mutate(splitby = get(split.by)) %>%
      unite(numerator, c(group.by, split.by, layer.by), sep=":", remove = F) %>%
      unite(denominator, normalize.by, sep=":", remove = F) %>%
      group_by(numerator) %>%
      summarise(count = n(), denominator = unique(denominator), display = unique(plot), splitby=unique(splitby)) %>%
      group_by(denominator)
  }

  mt <- mt %>%
    mutate(prop = round(count/sum(count), digits = label.digit)) %>%
    replace(is.na(.),0) %>%
    group_by(splitby) %>%
    mutate(prop = prop/sum(prop)) %>%
    separate(numerator, into = c("groupby","splitby","layerby"), sep = ":", fill = "right") %>%
    filter(display == "y")
  if (is.factor(object@meta.data[,group.by])) mt$groupby <- factor(mt$groupby, levels = levels(object@meta.data[,group.by]))
  if (is.factor(object@meta.data[,split.by])) mt$splitby <- factor(mt$splitby, levels = levels(object@meta.data[,split.by]))
  if (is.factor(object@meta.data[,layer.by])) mt$layerby <- factor(mt$layerby, levels = levels(object@meta.data[,layer.by]))

  bp<- ggplot(mt, aes(x=layerby,y=prop,fill=groupby)) +
    geom_bar(width = width, stat = "identity", position = "fill") +
    labs(fill=group.by)
  pie <- bp + coord_polar("y", start=0) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank())
  if (!is.null(split.by)) pie <- pie + facet_wrap(~splitby)
  if (number_labels == "value"){
    pie <- pie + geom_text(aes(label=count), position = position_fill(vjust = 0.5))
  }else if(number_labels == "percentage"){
    pie <- pie + geom_text(aes(label=scales::percent(prop)), position = position_fill(vjust = 0.5))
  }
  if (!is.null(plot.theme)) pie <- pie + eval(call(paste0("theme_",plot.theme)))
  if (!is.null(color)) pie <- pie + scale_fill_manual(values=color)
  return(pie)
}







