#' ColorPlot display multiple RGB colors
#'
#' @description To make picking RGB color set easier, this function display multiple colors in a matrix. Custom labels
#'
#' @param colors.ved a vector with RGB colors in format "#354C31"
#' @param n.col number of columns for the color matrix
#' @examples
#' colors.custom <- c("#354C31", "#354C45", "#313E4C",
#'                    "#342F4B", "#414749", "#3E6C72")
#' ColorPlot(colors.custom, n.col = 3, n.row = 2)
#' colors.celltype <- as.data.frame(matrix(data = c("KC","#3F837A",
#'                                                  "TC","#A55F6F",
#'                                                  "IC","#8659AC",
#'                                                  "MEL","#3E6D95"), ncol = 2, byrow = T))
#' ColorPlot(colors.celltype$V2, n.col = 2, n.row = 2, label = colors.celltype$V1)
#' @import ggplot2
#' @export


ColorPlot <- function(colors.vec, n.col, n.row, label=NULL){
  coor <- expand.grid(1:n.col, 1:n.row)
  coor <- coor[1:length(colors.vec),]
  df <- data.frame(col=colors.vec, x=coor$Var1, y=coor$Var2, lab=colors.vec)
  ind <- !duplicated(df$col)
  df$col <- factor(df$col, levels = df$col[ind])


  if (!is.null(label)){
    df$lab <- paste(df$lab, label, sep = "\n")
  }

  g <- ggplot(df, aes(x=x,y=y)) +
    geom_point(aes(col=col),size=15) +
    geom_text(aes(label=lab)) +
    scale_color_manual(values=levels(df$col)) +
    guides(color = guide_legend(override.aes = list(size=5)))
  g
  return(g)
}
