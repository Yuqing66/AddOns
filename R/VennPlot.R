
# input is a list of conditions.
# display number or element or error in each category
# shape = ellipse or circle
# elementPerRow, when too many elements, use multiple rows to make it narrower for display
library(eulerr)

VennPlot <- function(genelist, display="number", shape="ellipse", elementPerRow=5){
  geneall <- unique(unlist(genelist))
  mat <- as.data.frame(lapply(genelist, function(x) geneall %in% x))
  rownames(mat) <- geneall
  
  venn <- euler(mat, shape = shape)
  
  if (display == "number"){
    g <- plot(venn, quantities = venn$original.values)
  }else if(display == "element"){
    mat.group <- apply(mat, 1, function(x) paste(colnames(mat)[x], collapse = "&"))
    elements <- c()
    groups <- names(venn$original.values)
    for (group in groups){
      genes <- rownames(mat)[mat.group == group]
      tmp <- paste0(genes, collapse = ",")
      if (length(genes) > elementPerRow){
        tmp <- gsub(paste0("((?:[^,]+,){",elementPerRow-1,"}[^,]+),"), "\\1\n", tmp)
      }
      elements[group] <- tmp
    }
    g <- plot(venn, quantities = elements)
  }else if (display == "error"){
    g <- error_plot(venn)
  }
  return(g)
}


