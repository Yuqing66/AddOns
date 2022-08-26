#' ToMeta
#' @description A wrapper to change a list of genes to their meta gene names if applicable. 
#'  Return a genelist with new meta gene names and print the part of matagene table describing 
#'  meta genes included in gene list. Can be directly used for ploting with new data
#'  print info for meta genes included.
#' @param genes A vector of genes.
#' @param data The mata gene name table. If NULL, use the default table for FourDisease data.
#' @param to Change to meta or mergedgene name. Or use the name of a custom column.
#' @examples
#' FeaturePlot(srt, ToMeta(c("CD4","FOXP3")))
#' @export



ToMeta <- function(genes, data=NULL, to = "mergedname"){
  
  if (is.null(data)){
    data("metagene")
    data <- metagene
  }
  # are there meta genes in gene list?
  tmp <- sum(genes %in% data$name)
  if (tmp > 0){
    ind <- match(genes, data$name)
    if (!(to %in% colnames(data))) stop("parameter 'to' needs to be one of the columns in data")
    metanames <- data[ind[!is.na(ind)],to]
    genesinmeta <- genes[!is.na(ind)]
    genes[!is.na(ind)] <- metanames
    message(paste0(paste(genesinmeta, collapse = " ")," is/are in metagenes"))
    print(data[data[,to] %in% metanames, ])
  }else{
    message("Not in metagenes")
  }
  return(genes)
}

