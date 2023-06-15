#' Generate Marker Gene Evidence Matrix (CellMarker 2.0 DB)
#'
#' This internal function generates a marker gene evidence matrix based on the CellMarker 2.0 database and given filtering criteria.
#'
#' @param markerDB A data frame containing the marker gene database.
#' @param species A character vector specifying the species to consider. Currently it supports "human" and "mouse". Use "all" to include all species.
#' @param tissueClass A character vector specifying the broad tissue classes to consider (case insensitive). Default is "all" that includes all tissue classes. To explore the available tissue classes, use the \code{showCellMarkerTissue()} function.
#' @param tissueType A character vector specifying the specific tissue types to consider within the chosen tissue class (case insensitive). Default is "all" that includes all tissue types. "tissueType" represents the specific subtypes within "tissueClass", For instance, when considering blood as "tissueClass", examples of blood tissue types could be "Peripheral blood", "Serum", "Blood vessel", and so on. To explore the available tissue types, use the \code{showCellMarkerTissue()} function.
#' @param cellState A character vector specifying the cell states to consider. Valid options are "normal" and "cancer". Default is both.
#' @param markerSource A character vector specifying the marker gene sources to consider. Valid options are "experiment", "review", "single-cell sequencing", and "company". Default is all.
#' @param logBase The base value for logarithmic transformation. Default is the Euler's number.
#' @param pseudoCount The pseudo count value added to the raw counts before logarithmic transformation. Default is 1.
#'
#' @return A matrix containing marker gene evidence levels. Row names represent cell types, and column names represent genes. The higher the value, the stronger the evidence that the gene is a marker for the cell type.
#' 
#' @importFrom dplyr "%>%" n_distinct

cellMarkerSupport <- function(markerDB,
                              tissueClass="all",
                              tissueType="all",
                              species="human",
                              cellState=c("normal", "cancer"),
                              markerSource=c("Experiment", 
                                             "Review",
                                             "Single-cell sequencing",
                                             "Company"),
                              logBase=exp(1),
                              pseudoCount=1) {
    
    stopifnot(!missing(markerDB))
    stopifnot(!missing(species))
    
    speciesAll <- tolower(sort(unique(markerDB$species)))
    tissueClassAll <- tolower(sort(unique(markerDB$tissue_class)))
    tissueTypeAll <- tolower(sort(unique(markerDB$tissue_type)))
    
    markerDB <- markerDB[!is.na(markerDB$species) &
                             !is.na(markerDB$tissue_class) &
                             !is.na(markerDB$tissue_type) &
                             !is.na(markerDB$cell_type) &
                             !is.na(markerDB$Symbol) &
                             !is.na(markerDB$PMID), ]
    
    markerDB$species <- tolower(markerDB$species)
    markerDB$tissue_class <- tolower(markerDB$tissue_class)
    markerDB$tissue_type <- tolower(markerDB$tissue_type)
    markerDB$marker_source <- tolower(markerDB$marker_source)
    species <- tolower(species)
    tissueClass <- tolower(tissueClass)
    tissueType <- tolower(tissueType)
    markerSource <- tolower(markerSource)
    
    stopifnot(markerSource %in% c("experiment", "review",
                                  "single-cell sequencing","company"))
    
    stopifnot(cellState %in% c("normal", "cancer"))
    
    if("all" %in% tolower(species)) {
        species <- speciesAll
    }
    
    if("all" %in% tolower(tissueClass)) {
        tissueClass <- tissueClassAll
    }
    
    if("all" %in% tolower(tissueType)) {
        tissueType <- tissueTypeAll
    }
    
    if (!all(species %in% speciesAll)) {
        stop("Species", 
             paste(species[!species %in% speciesAll], collapse = ", "), 
             " not in the marker database.",
             " To explore the available species,",
             " use the showCellMarkerTissue() function.")
    }
    
    if (!all(tissueClass %in% tissueClassAll)) {
        stop("Tissue Class", paste(tissueClass[!tissueClass %in% tissueClassAll], collapse = ", "), 
             " not in the marker database.",
             " To explore the available tissue classes,",
             " use the showCellMarkerTissue() function.")
    }
    
    if (!all(tissueType %in% tissueTypeAll)) {
        stop("Tissue Type", paste(tissueType[!tissueType %in% tissueTypeAll], collapse = ", "), 
             " not in the marker database.",
             " To explore the available tissue types,",
             " use the showCellMarkerTissue() function.")
    }
    
    cellState[cellState %in% "normal"] <- "Normal cell"
    cellState[cellState %in% "cancer"] <- "Cancer cell"
    
    markerDB2Use <- markerDB[markerDB$marker_source %in% markerSource &
                                 markerDB$species %in% species &
                                 markerDB$tissue_class %in% tissueClass &
                                 markerDB$tissue_type %in% tissueType &
                                 markerDB$cell_type %in% cellState, 
                             c("species", "tissue_class", "tissue_type", "cell_type",
                               "cell_name", "Symbol", "PMID")]
    
    markerDB2Use$cellInfo <- paste(markerDB2Use$species,
                                   markerDB2Use$tissue_class,
                                   markerDB2Use$tissue_type,
                                   markerDB2Use$cell_type,
                                   markerDB2Use$cell_name, sep = ":")
    
    pmidCounts <- markerDB2Use %>% 
        group_by(cellInfo, Symbol) %>% 
        summarise(count = n_distinct(PMID))
    
    # Spread the data into a matrix format
    pmidCounts <- pmidCounts  %>% 
        tidyr::spread(Symbol, count)
    cellInfo <- pmidCounts$cellInfo
    
    pmidCounts <- as.matrix(pmidCounts[, -1])
    row.names(pmidCounts) <- cellInfo
    pmidCounts[is.na(pmidCounts)] <- 0
    
    markerEvi <- log(pmidCounts + pseudoCount, base = logBase)
    
    pmidCountsCellInfo <- markerDB2Use %>% 
        group_by(cellInfo) %>% 
        summarise(count = n_distinct(PMID))
    
    pmidCountsCellInfo <- 
        pmidCountsCellInfo[match(rownames(pmidCounts), 
                                 pmidCountsCellInfo$cellInfo), ]
    
    markerCountsCellInfo <- markerDB2Use %>% 
        group_by(cellInfo) %>% 
        summarise(count = n_distinct(Symbol))
    
    markerCountsCellInfo <- 
        markerCountsCellInfo[match(rownames(pmidCounts), 
                                   markerCountsCellInfo$cellInfo), ]
    
    cellInfoEvi <- log(pmidCountsCellInfo$count * markerCountsCellInfo$count +
                           pseudoCount, base = logBase)
    names(cellInfoEvi) <- rownames(pmidCounts)
    
    weight <- markerEvi * cellInfoEvi
    
    return(weight)
}