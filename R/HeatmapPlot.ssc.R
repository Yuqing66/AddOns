#' Heatmap with single cell level or mean expression value.
#'
#' @description
#' Plot a heatmap with single cell level or mean expression value to show the expression level of genes in the list.
#' @param object Seurat object.
#' @param genes A vector of genes for heatmap.
#' @param group.by Metadata column to aggregate for mean expression value. Null for single-cell.
#' @param split.by Split into columns.
#' @param color User defined color panel for ComplexHeatmap.
#' @param assay Seurat assay to get gene expression values.
#' @param slot Seurat slot of assy for gene expression values.
#' @param title Plot title.
#' @examples
#' HeatmapPlot.ssc(srt.ic, genes = c("CLEC9A","XCR1","IFI35","CLEC10A"), group.by = "DiseaseFiner")
#' @import ggplot2 dplyr ComplexHeatmap
#' @export



HeatmapPlot.ssc <- function(object, genes, group.by=NULL, split.by=NULL, color=NULL,
                        assay="RNA", slot="data", dataScale='row',
                        title=NULL, showRowName=T, showColName=T,
                        clusterRow=FALSE, clusterCol=F, reorderRow=F, reorderCol=F,
                        annRightDf=NULL, annTopShow=F, row_split_str=NULL){
  df <- GetAssayData(object, slot = slot, assay = assay)
  df <- df[genes, ]
  meta <- object@meta.data[,c(group.by, split.by), drop=F] %>%
    unite(condition, c(group.by, split.by), sep=":", remove = F)
  conditions <- meta %>%
    group_by_at(c(group.by, split.by)) %>%
    summarise(condition = unique(condition)) %>% pull(condition)
  df2 <- sapply(conditions, function(x){
    ind <- meta$condition == x
    tmp <- apply(df[,ind], 1, mean)
    return(tmp)
  })

  #### scale data ####
  if (dataScale == 'row'){
    df2 <- t(scale(t(df2)))
  }else if (dataScale == 'column'){
    df2 <- scale(df2)
  }


  #### split columns ####
  if (!is.null(split.by)){
    tmp <- unlist(strsplit(colnames(df2), ":"))
    groupbys <- tmp[c(T,F)]
    splitbys <- tmp[c(F,T)]
    column_split_str <- factor(splitbys)
    if (is.factor(object@meta.data[,split.by])){
      levels(column_split_str) <- levels(object@meta.data[,split.by])
    }
  }else{
    groupbys <- colnames(df2)
    column_split_str <- NULL
  }

  #### show column group color annotation ####
  if (annTopShow){
    topAnn <- HeatmapAnnotation(group=groupbys, which = 'column')
  }else{
    topAnn <- NULL
  }

  #### add gene annotation on the right ####
  if (!is.null(annRightDf)){
    rightAnn <- HeatmapAnnotation(df = annRightDf,
                                  which = 'row',
                                  na_col = "grey",
                                  foo = anno_empty(border = FALSE, width = GOanno.width))
    # show_legend = c(rep(FALSE, showCategory), FALSE)
  }else{
    rightAnn <- NULL
    # HeatmapAnnotation(foo = anno_empty(border = FALSE, width = GOanno.width))
  }

  #### plot heatmap ####
  hmap <- Heatmap(
    df2,
    col = color,
    # column_title = paste0("kmeans k = ",heatmap.k),
    name = "Expression",
    show_row_names = showRowName,
    row_names_side = "left",
    show_column_names = showColName,
    cluster_rows = clusterRow,
    cluster_columns = clusterCol,
    show_column_dend = F,
    show_row_dend = F,
    row_dend_reorder = reorderRow,
    column_dend_reorder = reorderCol,
    # column_order = c(1,7,8,9,10,11,2,3,4,5,6),
    right_annotation =  rightAnn,
    top_annotation = topAnn,
    column_split = column_split_str,
    row_split = row_split_str)

  return(hmap)
}
