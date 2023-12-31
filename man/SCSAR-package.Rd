% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCSAR.R
\docType{package}
\name{SCSAR-package}
\alias{SCSAR}
\alias{SCSAR-package}
\title{SCSAR: Cell-type annotation for single-cell RNA-seq data based on marker genes}
\description{
SCSAR is an R package designed to facilitate cell-type annotation for single-cell RNA-seq data. It takes as input the information on differential gene expression for each cell cluster. Leveraging cell type marker gene information from the CellMarker 2.0 and Cell Taxonomy databases, SCSAR employs scoring algorithms to identify and assign the most accurate and reliable cell-type annotations.
}
\details{
SCSAR is an R package designed to facilitate cell-type annotation for
single-cell RNA-seq data. It takes as input the information on differential
gene expression for each cell cluster. Leveraging cell type marker gene
information from the CellMarker 2.0
(\url{http://bio-bigdata.hrbmu.edu.cn/CellMarker})
and Cell Taxonomy databases (\url{https://ngdc.cncb.ac.cn/celltaxonomy/}),
SCSAR employs scoring algorithms to identify and assign the most accurate and
reliable cell-type annotations.
}
\seealso{
\code{\link{cellMarkerAnno}}

\code{\link{showCellMarkerTissue}}

\code{\link{cellTaxoAnno}}

\code{\link{showCellTaxoSpecies}}
}
\author{
Lanying Wei
}
