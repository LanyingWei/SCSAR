#' Cell Type Annotation based on CellMarker 2.0 Database
#' 
#' This function assigns cell types to cell clusters of single-cell RNA-seq data based on differential gene expression information and marker gene evidence calculated using the CellMarker 2.0 database.
#'
#' @param clusterDEGs A data frame containing information about differentially expressed genes for each cell cluster. It should have the following columns:
#' \describe{
#'   \item{cluster}{The cluster identifier or label.}
#'   \item{gene}{The gene symbol.}
#'   \item{p_val}{The p-value representing the significance of differential expression for the gene in the cluster.}
#'   \item{avg_log2FC}{The average log fold change indicating the magnitude of differential expression for the gene in the cluster.}
#' } 
#' @param species A character vector specifying the species to consider. Currently it supports "human" and "mouse". Use "all" to include all species.
#' @param tissueClass A character vector specifying the broad tissue classes to consider (case insensitive). Default is "all" that includes all tissue classes. To explore the available tissue classes, use the \code{showCellMarkerTissue()} function.
#' @param tissueType A character vector specifying the specific tissue types to consider within the chosen tissue class (case insensitive). Default is "all" that includes all tissue types. "tissueType" represents the specific subtypes within "tissueClass", For instance, when considering blood as "tissueClass", examples of blood tissue types could be "Peripheral blood", "Serum", "Blood vessel", and so on. To explore the available tissue types, use the \code{showCellMarkerTissue()} function.
#' @param cellState A character vector specifying the cell states to consider. Valid options are "normal" and "cancer". Default is both.
#' @param markerSource A character vector specifying the marker gene sources to consider. Valid options are "experiment", "review", "single-cell sequencing", and "company". Default is all.
#' @param outFile The path where the xlsx file containing the annotation results will be saved. For example, path = "directory/results.xlsx".
#' @param clusters A character vector specifying the clusters to analyze. By default, all clusters are considered (set to "all").
#' @param cutoffLFC The log fold change cutoff value for the column "avg_log2FC" in the \code{clusterDEGs} data frame. Only cell type marker genes with an average log fold change greater than or equal to this value are considered as supporting evidence for the corresponding cell type. Default is 1, which corresponds to a two-fold change in average.
#' @param cutoffPval The p-value cutoff for the column "p_val" in the \code{clusterDEGs} data frame. Only cell type marker genes with a p-value less than or equal to this value are considered as supporting evidence for the corresponding cell type. Default is 0.05.
#' @param topCell The maximum number of top-ranked cell types to report. Default is 5.
#' @param topMarker The maximum number of top-ranked marker genes contributing to the cell type annotation to report for each cell type. Default is 5.
#' @param panelty The penalty weight assigned to cell type marker genes that are not a detected marker of the cluster. A "detected marker" means that the gene meets both the cutoffPval and paneltyCutoffLFC thresholds. A higher absolute value of the penalty means a more severe punishment for cell types with some of the marker genes that are not overexpressed in the cluster. The default value is -0.1.
#' @param paneltyCutoffLFC The log fold change cutoff value for including genes in the penalty calculation. Cell types with marker genes having an average log fold change less than this value are penalized. The default value is the same as "cutoffLFC". The value of "paneltyCutoffLFC" should be less than or equal to "cutoffLFC".
#'
#' @return A list containing the annotation results for each cluster. Each element in the list corresponds to a cluster and includes a data frame. The data frame has annotated cell types as row names and two columns: "score" (representing the final score for the cell type) and "topMarkers" (containing the top marker genes and their scores).
#' @export
#'
#' @examples
#' # Example usage
#' data(pbmc3kDEGs)
#' annoRes <- cellMarkerAnno(clusterDEGs = pbmc3kDEGs, 
#'                           species = "human", 
#'                           tissueClass = "Blood",
#'                           tissueType = "Peripheral blood",
#'                           cutoffLFC = 0.585)
#'
#'
#' @seealso
#' \code{\link{showCellMarkerTissue}}, \code{\link{cellTaxoAnno}}
#'
cellMarkerAnno <- function(clusterDEGs, 
                           species="human",
                           tissueClass="all",
                           tissueType="all",
                           cellState=c("normal", "cancer"),
                           markerSource=c("Experiment", 
                                          "Review",
                                          "Single-cell sequencing",
                                          "Company"),
                           outFile=NA,
                           clusters="all",
                           cutoffLFC=1,
                           cutoffPval=0.05,
                           topCell=5,
                           topMarker=5,
                           panelty=-0.1,
                           paneltyCutoffLFC=cutoffLFC) {
    
    stopifnot(!missing(clusterDEGs))
    markerDB <- cellMarkerDB
    
    markerEviMat <- cellMarkerSupport(markerDB = markerDB, 
                                      species = species, 
                                      tissueClass =  tissueClass,
                                      tissueType = tissueType,
                                      cellState = cellState,
                                      markerSource = markerSource,
                                      logBase = exp(1),
                                      pseudoCount = 1)
    
    cellTypeRes <- cellTypeIdentifier(clusterDEGs = clusterDEGs,
                                      markerEviMat = markerEviMat,
                                      clusters=clusters,
                                      cutoffLFC=cutoffLFC,
                                      cutoffPval=cutoffPval,
                                      topCell=topCell,
                                      topMarker=topMarker,
                                      panelty=panelty,
                                      paneltyCutoffLFC=paneltyCutoffLFC)
    
    if (!is.na(outFile)) {
        saveCellTypeRes(cellTypeRes=cellTypeRes, path=outFile)
    }
    
    return(cellTypeRes)
}
    



#' Cell Type Annotation based on Cell Taxonomy Database
#' 
#' This function assigns cell types to cell clusters of single-cell RNA-seq data based on differential gene expression information and marker gene evidence calculated using the Cell Taxonomy database.
#'
#' @param clusterDEGs A data frame containing information about differentially expressed genes for each cell cluster. It should have the following columns:
#' \describe{
#'   \item{cluster}{The cluster identifier or label.}
#'   \item{gene}{The gene symbol.}
#'   \item{p_val}{The p-value representing the significance of differential expression for the gene in the cluster.}
#'   \item{avg_log2FC}{The average log fold change indicating the magnitude of differential expression for the gene in the cluster.}
#' } 
#' @param species A character vector specifying the species to consider. Currently it contains data on 34 species. Use "all" to include all species. Default is "Homo sapiens". To explore the available species, use the \code{showCellTaxoSpecies()} function.
#' @param outFile The path where the xlsx file containing the annotation results will be saved.
#' @param clusters A character vector specifying the clusters to analyze. By default, all clusters are considered (set to "all").
#' @param cutoffLFC The log fold change cutoff value for the column "avg_log2FC" in the \code{clusterDEGs} data frame. Only cell type marker genes with an average log fold change greater than or equal to this value are considered as supporting evidence for the corresponding cell type. Default is 1, which corresponds to a two-fold change in average.
#' @param cutoffPval The p-value cutoff for the column "p_val" in the \code{clusterDEGs} data frame. Only cell type marker genes with a p-value less than or equal to this value are considered as supporting evidence for the corresponding cell type. Default is 0.05.
#' @param topCell The maximum number of top-ranked cell types to report. Default is 5.
#' @param topMarker The maximum number of top-ranked marker genes contributing to the cell type annotation to report for each cell type. Default is 5.
#' @param panelty The penalty weight assigned to cell type marker genes that are not a detected marker of the cluster. A "detected marker" means that the gene meets both the cutoffPval and paneltyCutoffLFC thresholds. A higher absolute value of the penalty means a more severe punishment for cell types with some of the marker genes that are not overexpressed in the cluster. The default value is -0.1.
#' @param paneltyCutoffLFC The log fold change cutoff value for including genes in the penalty calculation. Cell types with marker genes having an average log fold change less than this value are penalized. The default value is the same as "cutoffLFC". The value of "paneltyCutoffLFC" should be less than or equal to "cutoffLFC".
#'
#' @return A list containing the annotation results for each cluster. Each element in the list corresponds to a cluster and includes a data frame. The data frame has annotated cell types as row names and two columns: "score" (representing the final score for the cell type) and "topMarkers" (containing the top marker genes and their scores).
#' @export
#'
#' @examples
#' # Example usage
#' data(pbmc3kDEGs)
#' annoRes <- cellTaxoAnno(clusterDEGs = pbmc3kDEGs, 
#'                         species = "Homo sapiens",
#'                         cutoffLFC = 0.585,
#'                         clusters = "all")
#'
#'
#' @seealso
#' \code{\link{showCellTaxoSpecies}}, \code{\link{cellMarkerAnno}}
#'

cellTaxoAnno <- function(clusterDEGs,
                         species="Homo sapiens",
                         outFile=NA,
                         clusters="all",
                         cutoffLFC=1,
                         cutoffPval=0.05,
                         topCell=5,
                         topMarker=5,
                         panelty=-0.1,
                         paneltyCutoffLFC=cutoffLFC) {
    
    
    stopifnot(!missing(clusterDEGs))
    markerDB <- cellTaxoDB
    
    markerEviMat <- cellTaxonomySupport(markerDB = markerDB,
                                        species=species,
                                        logBase = exp(1),
                                        pseudoCount = 1)
    
    cellTypeRes <- cellTypeIdentifier(clusterDEGs = clusterDEGs,
                                      markerEviMat = markerEviMat,
                                      clusters=clusters,
                                      cutoffLFC=cutoffLFC,
                                      cutoffPval=cutoffPval,
                                      topCell=topCell,
                                      topMarker=topMarker,
                                      panelty=panelty,
                                      paneltyCutoffLFC=paneltyCutoffLFC)
    
    if (!is.na(outFile)) {
        saveCellTypeRes(cellTypeRes=cellTypeRes, path=outFile)
    }
    
    return(cellTypeRes)
}