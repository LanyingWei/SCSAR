#' PBMC 3k Differentially Expressed Genes (DEGs)
#'
#' The `pbmc3kDEGs` data is a data frame containing the differentially expressed genes (DEGs) for each cluster in the peripheral blood mononuclear cell (PBMC) dataset. The dataset consists of 3k PBMCs from a Healthy Donor and is freely available from 10x Genomics. The data can be accessed at: \url{https://www.10xgenomics.com/resources/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-0-0}.
#'
#' The DEGs were identified through normalization, clustering, and differential gene expression analyses using the R package Seurat 4.1. Clustering of cells was performed using the implemented community identification method in the 'FindClusters' function. The specific marker genes of cell clusters (the DEGs for a specific cluster) were identified by comparing cells from that cluster to all other cells using the Seurat 'FindMarkers' function.
#'
#' @format A data frame with the following columns:
#'   \describe{
#'     \item{cluster}{Cluster identifier.}
#'     \item{gene}{Gene symbol.}
#'     \item{p_val}{p-value for differential expression.}
#'     \item{avg_log2FC}{Average log2 fold change.}
#'   }
#'
#' @usage data(pbmc3kDEGs)
#'
"pbmc3kDEGs"