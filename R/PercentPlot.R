#' percent plot with denominator and divider specified.
#'
#' @description
#' This function will plot a dotplot with each dot representing the percentage of one type of cells in another type of cells.
#'
#' @param object Seurat object.
#' @param numerator One type to calculate percentage for, e.g. CD4
#' @param numerator.column Name of the metadata column where the numerator string is stored.
#' @param denominator A vector with all types to be used as total cells for percentage calculation.
#' @param denominator.column Name of the metadata column where the denominator strings are stored. Same as numerator.column if not defined.
#' @param group.by Group and color by this meta column
#' @param split.by Split into different panels. A vector with one or two elements, in order of c(row, col)
#' @param point.by "Patient" by default.
#' @param num.cutoff Color dots with fewer cells in denominator than num.cutoff in grey.
#' @param num.cutoff.rm Keep the grey dots not passing the num.cutoff threshold or remove them.
#' @param title Plot title.
#' @param number_labels Show stats at the bottom of the plot if TRUE. Total cell number in that condition.
#' @param spread Control skin is only in HealthyControl in Disease meta column, but we want to compare them to other skins in all diseases
#' even when we split by Disease. So we spread HealthyControl condition in Disease column, by spread = c("Disease", "HealthyControl").
#' @param theme ggplot theme. classic by default.
#' @param jitter.alpha Alpha of dots.
#' @param jitter.width Spread out width of dots.
#' @param jitter.size Size of dots.
#' @param jitter.col Whether to color dots based on number of cells within, or leave dots black.
#' @param col.trans Set to "log" to make color scale zoom in to lower values
#' @param bar Plot weighted mean as a bar if TRUE.
#' @param link Link the paired conditions in each panel. However, the starting and ending points cannot adjust according to
#' the horizontal location of jitter plot. Can set a lower value for jitter.width to shrink the blank distance.

#' @examples
#' PercentPlot.single(srt, numerator="CD8_RDH10", numerator.column="subCellType.ificc", denominator="CD8", denominator.column="subCellType.ificc.main", group.by = "Skin", split.by = "DiseaseFiner", ttest=F, link = T, num.cutoff.rm=T)
#' @import ggplot2
#' @export

PercentPlot.single <- function(
object, numerator, numerator.column, denominator, denominator.column = NULL,
group.by, split.by = NULL, point.by = "Patient", spread = NULL, num.cutoff = 20, num.cutoff.rm = F,
number.labels = T, xaxis.labels = F, title = "", theme = "classic",
text.size.plot.title = 20, text.size.axis.title = 10, text.size.axis.text = 10,
text.size.legend.title = 10, text.size.legend.text = 5, text.size.panel = 5,
jitter.alpha = 0.5, jitter.width = 0.1, jitter.size = 2, jitter.col = T, col.trans = "identity",
bar = T, link = FALSE, ttest = F){
# get meta data
meta <- data.frame(numer=object@meta.data[,numerator.column],
denom=object@meta.data[,denominator.column],
group=object@meta.data[,group.by],
split=object@meta.data[,split.by],
point=object@meta.data[,point.by])
colnames(meta) <- c("numer","denom",group.by, split.by, point.by)

# keep the order
group.by.levels <- levels(object@meta.data[, group.by])
if (!is.null(group.by.levels)) meta[, group.by] <- factor(meta[, group.by], levels = group.by.levels)

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

# calculate proportion and number of cells in condition
meta.dot <- meta %>% filter(denom %in% denominator) %>%
group_by_at(vars(one_of(c(group.by, split.by, point.by)))) %>%
# group_by(group, split, point) %>%
summarise(prop=sum(numer %in% numerator)/length(denom),
denomNum=length(denom)) %>%
mutate(numColor=replace(denomNum, which(denomNum < num.cutoff), NA))

# remove dots not passing cell number threshold from downstream process.
# mean calculation, plotting and t-test
if (num.cutoff.rm){
meta.dot <- meta.dot %>% filter(!is.na(numColor))
}

# spread Control skin to disease conditions
# spread = c("DiseaseFiner", "HealthyControl")
if (any(!is.null(spread))) {
others <- setdiff(unique(meta[, spread[1]]), spread[2])
ind <- which(meta.dot[, spread[1]] == spread[2])
rmmeta <- meta.dot[ind, ]
meta.dot <- meta.dot[-ind, ]
for (i in 1:length(others)) {
rmmeta[, spread[1]] <- others[i]
meta.dot <- rbind(meta.dot, rmmeta)
}
}

# calculate proportion and number of cells per bar
meta.bar <- meta.dot %>% group_by_at(vars(one_of(c(group.by, split.by)))) %>%
summarise(prop=round(sum(prop*denomNum)/sum(denomNum), 3),
denomNum=sum(denomNum))

g <- ggplot(meta.dot)
  if (number.labels == T) {
g <- g + geom_text(data = meta.bar, aes_string(x = group.by, y = -max(meta.bar$prop)/25, label = "prop"), size = 2.5)
g <- g + geom_text(data = meta.bar, aes_string(x = group.by, y = -max(meta.bar$prop)/10, label = "denomNum"), size = 2.5)
}
  if (!is.null(theme)) g <- g + eval(call(paste0("theme_",theme)))
  if (title == "") title <- paste0("Percentage of ", numerator, " in ", denominator)
g <- g + labs(title = title, y = "Percentage")
g <- g + theme(plot.title = element_text(size = text.size.plot.title),
axis.title = element_text(size = text.size.axis.title), axis.text = element_text(size = text.size.axis.text),
legend.title = element_text(size = text.size.legend.title), legend.text = element_text(size = text.size.legend.text))
g <- g + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  if (jitter.col){
g <- g + geom_jitter(aes_string(x = group.by, y = "prop", col = "numColor"),
width = jitter.width, size = jitter.size, alpha = jitter.alpha)
g <- g + scale_color_gradientn(colors = c("blue", "red", "orange"), trans=col.trans)
}else{
g <- g + geom_jitter(aes_string(x = group.by, y = "prop"),
width = jitter.width, size = jitter.size, alpha = jitter.alpha)
}
g <- g + scale_y_continuous(labels = scales::percent)
  if (!(xaxis.labels)) g <- g + theme(axis.text.x = element_blank())
  if (bar) {
g <- g + geom_col(data = meta.bar, aes_string(x = group.by, y = "prop", fill = group.by), alpha = 0.3)
}
  if (link) {
g <- g + geom_line(aes_string(x = group.by, y = "prop", group = point.by), alpha = 0.8)
}
  if (length(split.by) == 1) {
g <- g + facet_grid(facets = reformulate(split.by),
scales = "free_x", space = "free_x") +
theme(strip.text.x = element_text(size = text.size.panel))
}else if (length(split.by) == 2) {
g <- g + facet_grid(facets = reformulate(split.by[1], split.by[2]),
scales = "free_x", space = "free_x") +
theme(strip.text.x = element_text(size = text.size.panel),
strip.text.y = element_text(size = text.size.panel))
}else if (length(split.by) > 2) {
stop("Parameter split.by needs to be a string with equal or less than two variables.")
}
  if (ttest == F){
return(g)
}else{
    # do t-test on dots pass filter (with more than a certain number of cells)
df.ttest <- as.data.frame(meta.dot)
list.ttest <- list()
    if (is.null(split.by)){
      # t-test all
list.ttest[[1]] <- df.ttest
}else{
      # plot split in multiple panels. Do t-test across groups with in each panel
list.ttest <- split(df.ttest, df.ttest[,split.by])
}
ttestresult <- vector("list", length(list.ttest))
names(ttestresult) <- names(list.ttest)
    for (i.tb in 1:length(list.ttest)){
      tb <- list.ttest[[i.tb]]
      groups <- unique(as.character(tb[,group.by]))
      if (length(groups) == 1) {
        ttestresult <- ttestresult[[i.tb]] <- NULL
        next
      }
      comparisons <- as.data.frame(t(combn(groups,2)))
      colnames(comparisons) <- c("group1","group2")
      t.res <- data.frame()
      for (i in 1:nrow(comparisons)){
        tb2 <- split(tb, as.character(tb[,group.by]))
        if (nrow(tb2[[1]]) == 1 | nrow(tb2[[2]]) == 1) next
        if (link){
          # if link is true, then do paired t-test with connected points.
          pointpaired <- intersect(tb2[[1]][,point.by], tb2[[2]][,point.by])
          if (length(pointpaired) < 2) message("Less than 3 points paired")
          tb2[[1]] <- tb2[[1]][match(pointpaired, tb2[[1]][,point.by]),]
          tb2[[2]] <- tb2[[2]][match(pointpaired, tb2[[2]][,point.by]),]
          tmp <- t.test(x=tb2[[comparisons[i,1]]][,"prop"], y=tb2[[comparisons[i,2]]][,"prop"], paired = T)

        }else{
          # do with all samples with regular t-test
          tmp <- t.test(x=tb2[[comparisons[i,1]]][,"prop"], y=tb2[[comparisons[i,2]]][,"prop"], paired = F)
        }
        t.res <- rbind(t.res,
                       data.frame(t = tmp$statistic,
                                  df = tmp$parameter,
                                  p = tmp$p.value,
                                  confint.low = tmp$conf.int[1],
                                  confint.high = tmp$conf.int[2],
                                  mean.1 = tmp$estimate[1],
                                  mean.2 = tmp$estimate[2],
                                  stderr = tmp$stderr))
      }
      if (nrow(t.res) > 0){
        res <- cbind(comparisons, t.res)
        ttestresult[[i.tb]] <- res
      }else{
        ttestresult[[i.tb]] <- NULL
      }
    }
    names(ttestresult) <- names(list.ttest)
    output <- list(g=g, ttest=ttestresult)
    return(output)
  }

}






