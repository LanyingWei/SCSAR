% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cellTypeAnnotator.R
\name{cellTaxoAnno}
\alias{cellTaxoAnno}
\title{Cell Type Annotation based on Cell Taxonomy Database}
\usage{
cellTaxoAnno(
  clusterDEGs,
  species = "Homo sapiens",
  outFile = NA,
  clusters = "all",
  cutoffLFC = 1,
  cutoffPval = 0.05,
  topCell = 5,
  topMarker = 5,
  panelty = -0.1,
  paneltyCutoffLFC = cutoffLFC
)
}
\arguments{
\item{clusterDEGs}{A data frame containing information about differentially expressed genes for each cell cluster. It should have the following columns:
\describe{
\item{cluster}{The cluster identifier or label.}
\item{gene}{The gene symbol.}
\item{p_val}{The p-value representing the significance of differential expression for the gene in the cluster.}
\item{avg_log2FC}{The average log fold change indicating the magnitude of differential expression for the gene in the cluster.}
}}

\item{species}{A character vector specifying the species to consider. Currently it contains data on 34 species. Use "all" to include all species. Default is "Homo sapiens". To explore the available species, use the \code{showCellTaxoSpecies()} function.}

\item{outFile}{The path where the xlsx file containing the annotation results will be saved.}

\item{clusters}{A character vector specifying the clusters to analyze. By default, all clusters are considered (set to "all").}

\item{cutoffLFC}{The log fold change cutoff value for the column "avg_log2FC" in the \code{clusterDEGs} data frame. Only cell type marker genes with an average log fold change greater than or equal to this value are considered as supporting evidence for the corresponding cell type. Default is 1, which corresponds to a two-fold change in average.}

\item{cutoffPval}{The p-value cutoff for the column "p_val" in the \code{clusterDEGs} data frame. Only cell type marker genes with a p-value less than or equal to this value are considered as supporting evidence for the corresponding cell type. Default is 0.05.}

\item{topCell}{The maximum number of top-ranked cell types to report. Default is 5.}

\item{topMarker}{The maximum number of top-ranked marker genes contributing to the cell type annotation to report for each cell type. Default is 5.}

\item{panelty}{The penalty weight assigned to cell type marker genes that are not a detected marker of the cluster. A "detected marker" means that the gene meets both the cutoffPval and paneltyCutoffLFC thresholds. A higher absolute value of the penalty means a more severe punishment for cell types with some of the marker genes that are not overexpressed in the cluster. The default value is -0.1.}

\item{paneltyCutoffLFC}{The log fold change cutoff value for including genes in the penalty calculation. Cell types with marker genes having an average log fold change less than this value are penalized. The default value is the same as "cutoffLFC". The value of "paneltyCutoffLFC" should be less than or equal to "cutoffLFC".}
}
\value{
A list containing the annotation results for each cluster. Each element in the list corresponds to a cluster and includes a data frame. The data frame has annotated cell types as row names and two columns: "score" (representing the final score for the cell type) and "topMarkers" (containing the top marker genes and their scores).
}
\description{
This function assigns cell types to cell clusters of single-cell RNA-seq data based on differential gene expression information and marker gene evidence calculated using the Cell Taxonomy database.
}
\examples{
# Example usage
data(pbmc3kDEGs)
annoRes <- cellTaxoAnno(clusterDEGs = pbmc3kDEGs, 
                        species = "Homo sapiens",
                        cutoffLFC = 0.585,
                        clusters = "all")


}
\seealso{
\code{\link{showCellTaxoSpecies}}, \code{\link{cellMarkerAnno}}
}
