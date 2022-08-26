#' DotPlot.ssc with sample mean as each dot.
#'
#' @description To evaluate a gene's consistency across patients, it's useful to know the expression
#' level of individual samples. If set point.by to Patient, and use group.by and split.by to separate conditions, each
#' dot in this plot would represent the average mean of each patient. Color by the number of cells used for mean calculation,
#' indicating how trustable they are. Blue at the lower end, then red to orange. Set bar to TRUE to plot the condition mean
#' with all cells from all samples included.
#'
#' @details Use mean of normalized UMI of cells in a condition instead of aggregated bulk to simplify calculation.
#'
#' @param input Seurat object.
#' @param gene
#' @param group.by Group and color by this meta column
#' @param split.by Split into different panels. A vector with one or two elements, in order of c(row, col)
#' @param point.by "Patient" by default.
#' @param num.cutoff Color dots with lower than num.cutoff cells to grey.
#' @param title Plot title.
#' @param number_labels Show stats at the bottom of the plot if TRUE. First row is total cell number in each condition,
#' second row is fraction of cells with at least one UMI of this gene detected.
#' @param text_sizes Seven elements in order: plot title, axis title, axis text, legene title, legend text, split, stat size.
#' @param spread Control skin is only in HealthyControl in Disease meta column, but we want to compare them to other skins in all diseases
#' even when we split by Disease. So we spread HealthyControl condition in Disease column, by spread = c("Disease", "HealthyControl").
#' @param theme ggplot theme. classic by default.
#' @param jitter.alpha
#' @param jitter.width
#' @param jitter.size
#' @param col.trans Set to "log" to make color scale zoom in to lower values
#' @param bar Plot weighted mean as a bar if TRUE.
#' @param link Link the paired conditions in each panel. However, the starting and ending points cannot adjust according to
#' the horizontal location of jitter plot. Can set a lower value for jitter.width to shrink the blank distance.
#' @examples
#' srt$Skin <- factor(srt$Skin, levels = c("Control", "NonLesional", "Lesional"))
#' DotPlot.ssc(srt, "FOSL1", group.by = "Skin", split.by = c("subCellType","DiseaseFiner"), spread = c("DiseaseFiner","HealthyControl"))
#' @import ggplot2
#' @export




DotPlot.ssc <- function(input, gene, group.by, split.by = NULL, point.by = "Patient", num.cutoff = 10, title = "",
                        number_labels = T, text_sizes = c(20, 10, 5, 10, 5, 5, 3), spread = NULL,
                        theme = "classic", jitter.alpha = 0.5, jitter.width = 0.2, jitter.size = 2, col.trans = "", bar = T,
                        link = FALSE)
{
  # keep the order
  group.by.levels <- levels(input@meta.data[,group.by])

  # mean of RNA assay, normalized UMI, in each condition
  meta <- input@meta.data[, c(group.by, split.by, point.by)]

  # data frame for plotting
  tmp <- tapply(input@assays$RNA@data[gene,],
                apply(meta, 1, function(x) paste0(x, collapse = ":")),
                function(x) c(mean(x), length(x), sum(x>0)))

  df <- data.frame(gene = sapply(tmp, function(x) x[1]),
                   num = sapply(tmp, function(x) x[2]),
                   num_above0 = sapply(tmp, function(x) x[3]))
  df$num_cutoff <- df$num
  df$num_cutoff[df$num < num.cutoff] <- NA

  df.meta <- as.data.frame(t(as.data.frame(sapply(names(tmp), function(x) strsplit(x,":")))))
  colnames(df.meta) <- c(group.by, split.by, point.by)
  df <- cbind(df, df.meta)

  # data frame for stat text of fraction cells with UMI detected.
  tmp <- tapply(input@assays$RNA@data[gene,],
                apply(meta[,c(group.by, split.by)], 1, function(x) paste0(x, collapse = ":")),
                function(x) c(mean(x), round(sum(x>0)/length(x), digits = 3)))

  df.stat <- data.frame(bar_mean = sapply(tmp, function(x) x[1]),
                        frac_sum = sapply(tmp, function(x) x[2]))
  df.meta <- as.data.frame(t(as.data.frame(sapply(names(tmp), function(x) strsplit(x,":")))))
  colnames(df.meta) <- c(group.by, split.by)
  df.stat <- cbind(df.stat, df.meta)

  # spread Control skin to disease conditions
  if (any(!is.null(spread))) {
    others <- setdiff(unique(df[, spread[1]]), spread[2])
    ind <- which(df[, spread[1]] == spread[2])
    rmdf <- df[ind, ]
    df <- df[-ind, ]

    ind.stat <- which(df.stat[, spread[1]] == spread[2])
    rmdf.stat <- df.stat[ind.stat, ]
    df.stat <- df.stat[-ind.stat, ]
    for (i in 1:length(others)) {
      rmdf[, spread[1]] <- others[i]
      df <- rbind(df, rmdf)
      rmdf.stat[, spread[1]] <- others[i]
      df.stat <- rbind(df.stat, rmdf.stat)
    }
  }

  df[,group.by] <- factor(df[,group.by], levels = group.by.levels)
  df.stat[,group.by] <- factor(df.stat[,group.by], levels = group.by.levels)

  g <- ggplot(df)

  if (number_labels == T) {
    g <- g + stat_summary(aes_string(x = group.by, y = "num"),
                          fun.data = function(x) {
                            return(c(y = -max(df$gene)/25, label = sum(x)))
                          }, colour = "black", geom = "text", size = text_sizes[7])
    g <- g + geom_text(data = df.stat, aes_string(x = group.by, y = -max(df$gene)/10, label = "frac_sum"), size = text_sizes[7])
  }

  if (theme == "bw") {
    g <- g + theme_bw()
  }else {
    g <- g + theme_classic()
  }

  ###### ###### What is changed from the original package.
  if (title == "") {
    title <- gene
    g <- g + labs(title = title, y = gene)
  }else {
    g <- g + labs(title = title, y = gene)
  }
  ###### ######

  g <- g + theme(plot.title = element_text(size = text_sizes[1]),
                 axis.title = element_text(size = text_sizes[2]), axis.text = element_text(size = text_sizes[3]),
                 legend.title = element_text(size = text_sizes[4]), legend.text = element_text(size = text_sizes[5]))
  g <- g + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))

  # g <- g + geom_col(aes_string(x = group.by, y = "weighted_mean", fill = group.by), col = "black", alpha = alpha)
  # g <- g + geom_dotplot(aes_string(x = group.by, y = "gene", fill = group.by), binaxis='y', stackdir = 'center', stackratio = stackratio, dotsize = dotsize, alpha = 0.8, binwidth = (max(df[,"gene"])-min(df[,"gene"]))*binwidth)

  g <- g + geom_jitter(aes_string(x = group.by, y = "gene", col = "num_cutoff"), width = jitter.width, size = jitter.size, alpha = jitter.alpha)
  g <- g + scale_color_gradientn(colors = c("blue","red","orange"))

  if (bar) {g <- g + geom_col(data = df.stat, aes_string(x = group.by, y = "bar_mean", fill = group.by), alpha = 0.3)}
  if (link) {g <- g + geom_line(aes_string(x = group.by, y = "gene", group = point.by), alpha = 0.8)}


  if (length(split.by) == 1){
    g <- g + facet_grid(facets = reformulate(split.by),
                        scales = "free_x", space = "free_x")
  }else if(length(split.by) == 2){
    g <- g + facet_grid(facets = reformulate(split.by[1],split.by[2]),
                        scales = "free_x", space = "free_x")
  }else if (length(split.by) > 2) {
    stop("Parameter split.by needs to be a string with equal or less than two variables.")
  }
  g <- g + theme(strip.text.x = element_text(size = text_sizes[6]))

  # g <- g + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  return(g)
}



