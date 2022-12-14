% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FeaturePlot.meta.R
\name{FeaturePlot.meta}
\alias{FeaturePlot.meta}
\title{FeaturePlot.meta, FeaturePlot for genelist with metagenes in.}
\usage{
FeaturePlot.meta(
  object,
  features,
  data.ToMeta = NULL,
  to.ToMeta = "mergedname",
  dims = c(1, 2),
  cells = NULL,
  cols = if (blend) {     c("lightgrey", "#ff0000", "#00ff00") } else {    
    c("lightgrey", "blue") },
  pt.size = NULL,
  order = FALSE,
  min.cutoff = NA,
  max.cutoff = NA,
  reduction = NULL,
  split.by = NULL,
  keep.scale = "feature",
  shape.by = NULL,
  slot = "data",
  blend = FALSE,
  blend.threshold = 0.5,
  label = FALSE,
  label.size = 4,
  label.color = "black",
  repel = FALSE,
  ncol = NULL,
  coord.fixed = FALSE,
  by.col = TRUE,
  sort.cell = NULL,
  interactive = FALSE,
  combine = TRUE,
  raster = NULL,
  raster.dpi = c(512, 512)
)
}
\arguments{
\item{data.ToMeta}{meta gene table passed to ToMeta function.}

\item{to.ToMeta}{to parameter in ToMeta function. The gene name format used in this Seurat object.}

\item{...}{same as Seurat FeaturePlot function.}
}
\description{
Although ToMeta can easily transform any gene to its meta name, the saved FeaturePlot can be confusing not
knowing what was the gene name wanted to plot at first. This function embed ToMeta inside, and display the original gene name.
}
\examples{
FeaturePlot.meta(srt, "FOXP3")
}
