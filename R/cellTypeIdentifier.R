#' Cell Type Annotation
#'
#' This is an internal function that assigns cell types to cell clusters of single-cell RNA-seq data based on differential gene expression information and marker gene evidence.
#'
#' @param clusterDEGs A data frame containing information about differentially expressed genes for each cell cluster. It should have the following columns:
#' \describe{
#'   \item{cluster}{The cluster identifier or label.}
#'   \item{gene}{The gene symbol.}
#'   \item{p_val}{The p-value representing the significance of differential expression for the gene in the cluster.}
#'   \item{avg_log2FC}{The average log fold change indicating the magnitude of differential expression for the gene in the cluster.}
#' }
#' @param markerEviMat A matrix containing the evidence levels of marker genes. Rows represent cell types, and columns represent genes. The values in the matrix indicate the strength of evidence that a gene is a marker for a particular cell type. A higher value indicates stronger evidence.
#' @param clusters A character vector specifying the clusters to analyze. By default, all clusters are considered (set to "all").
#' @param cutoffLFC The log fold change cutoff value for the column "avg_log2FC" in the \code{clusterDEGs} data frame. Only cell type marker genes with an average log fold change greater than or equal to this value are considered as supporting evidence for the corresponding cell type. Default is 1, which corresponds to a two-fold change in average.
#' @param cutoffPval The p-value cutoff for the column "p_val" in the \code{clusterDEGs} data frame. Only cell type marker genes with a p-value less than or equal to this value are considered as supporting evidence for the corresponding cell type. Default is 0.05.
#' @param topCell The maximum number of top-ranked cell types to report. Default is 5.
#' @param topMarker The maximum number of top-ranked marker genes contributing to the cell type annotation to report for each cell type. Default is 5.
#' @param panelty The penalty weight assigned to cell type marker genes that are not a detected marker of the cluster. A "detected marker" means that the gene meets both the cutoffPval and paneltyCutoffLFC thresholds. A higher absolute value of the penalty means a more severe punishment for cell types with some of the marker genes that are not overexpressed in the cluster. The default value is -0.1.
#' @param paneltyCutoffLFC The log fold change cutoff value for including genes in the penalty calculation. Cell types with marker genes having an average log fold change less than this value are penalized. The default value is the same as "cutoffLFC". The value of "paneltyCutoffLFC" should be less than or equal to "cutoffLFC".
#'
#' @return A list containing the annotation results for each cluster. Each element in the list corresponds to a cluster and includes a data frame. The data frame has annotated cell types as row names and two columns: "score" (representing the final score for the cell type) and "topMarkers" (containing the top marker genes and their scores).
#'
#' @importFrom stats sd

cellTypeIdentifier <- function(clusterDEGs, 
                               markerEviMat, 
                               clusters="all",
                               cutoffLFC=1,
                               cutoffPval=0.05,
                               topCell=5,
                               topMarker=5,
                               panelty=-0.1,
                               paneltyCutoffLFC=cutoffLFC) {
    
    stopifnot(!missing(clusterDEGs))
    stopifnot(!missing(markerEviMat))
    
    if("all" %in% tolower(clusters)) {
        clusters <- sort(unique(clusterDEGs$cluster))
    }
    
    res <- list()
    for (i in seq_along(clusters)) {
        cluster <- clusters[i]
        thisClusterDEGsAll <- clusterDEGs[clusterDEGs$cluster %in% cluster &
                                              clusterDEGs$gene %in% 
                                              colnames(markerEviMat), 
                                          c("gene", "p_val", "avg_log2FC")]
        thisClusterDEGs <- 
            thisClusterDEGsAll[thisClusterDEGsAll$p_val <= cutoffPval &
                                   thisClusterDEGsAll$avg_log2FC >= cutoffLFC,]
        
        if (nrow(thisClusterDEGs) == 0) {
            stop("For cluster ", clusters[i], ", ",
                 sum(thisClusterDEGsAll$p_val <= cutoffPval), 
                 " gene(s) meet the CutoffPval threshold, and ",   
                 sum(thisClusterDEGsAll$avg_log2FC >= cutoffLFC), 
                 " gene(s) meet the PaneltyCutoffLFC threshold.",
                 " There are no genes left after filtering for both criteria,",
                 " please adjust the thresholds and run again.")
        }
        
        
        clusterMarkers <- thisClusterDEGs$gene
        message("Estimating identity of cell cluster ", cluster,
                " based on ", length(clusterMarkers), 
                " cluster marker genes...")
        
        # Retain only cell types that have at least one marker gene that is DEG
        thisMarkerEviMat <- markerEviMat[, clusterMarkers, drop = FALSE]
        thisMarkerEviMat <- thisMarkerEviMat[rowSums(thisMarkerEviMat) > 0, , 
                                             drop = FALSE]
        
        # Remove genes that are not the marker gene for any remaining cell types
        thisMarkerEviMat <-  markerEviMat[rownames(thisMarkerEviMat), ]
        thisMarkerEviMat <- thisMarkerEviMat[, colSums(thisMarkerEviMat) > 0]
        
        # Initialize weight vector with penalty score
        weightVec <- rep(panelty, ncol(thisMarkerEviMat))
        names(weightVec) <- colnames(thisMarkerEviMat)
        
        # Assign weights for supporting marker genes
        weightVec[thisClusterDEGs$gene] <- thisClusterDEGs$avg_log2FC
        
        # For genes with avg_log2FC above paneltyCutoffLFC, set penalty to 0,
        # so that these genes neither add to nor lower the cell type score
        noPaneltyGene <- 
            thisClusterDEGsAll$gene[thisClusterDEGsAll$gene %in% 
                                        names(weightVec) &
                                        thisClusterDEGsAll$p_val <= cutoffPval &
                                        thisClusterDEGsAll$avg_log2FC >= 
                                        paneltyCutoffLFC &
                                        thisClusterDEGsAll$avg_log2FC < 
                                        cutoffLFC]
        weightVec[noPaneltyGene] <- 0
        
        score <- thisMarkerEviMat %*% weightVec
        colnames(score) <- "score"
        score <- score[order(score, decreasing = TRUE), , drop = FALSE]
        zscore <- (score - mean(score)) / sd(score)
        zscore <- zscore[zscore > 0, , drop = FALSE]
        zscore <- zscore[seq_len(min(topCell, nrow(zscore))), , drop = FALSE]
        zscore <- round(zscore, 2)
        
        markerGeneScore <- 
            t(t(thisMarkerEviMat[rownames(zscore), ]) * weightVec)
        markerGenes <- 
            apply(markerGeneScore, 1, function(x) {
                markerScore <- sort(x[x>0], decreasing = TRUE)
                markerScore <- 
                    markerScore[seq_len(min(topMarker, length(markerScore)))]
                paste(paste(names(markerScore), 
                            round(markerScore, 0), sep = " "), collapse = ", ")
            })
        resInfo <- as.data.frame(zscore)
        resInfo$topMarkers <- markerGenes
        res[[paste("cluster", cluster)]] <- resInfo
    }
    message("Done.")
    return(res)
}
