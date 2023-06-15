#' Generate Marker Gene Evidence Matrix (Cell Taxonomy DB)
#'
#' This internal function generates a marker gene evidence matrix based on the Cell Taxonomy database and given filtering criteria.
#'
#' @param markerDB A data frame containing the marker gene database.
#' @param species A character vector specifying the species to consider. Currently it contains data on 34 species. Use "all" to include all species. Default is "Homo sapiens". To explore the available species, use the \code{showCellTaxoSpecies()} function.
#' @param logBase The base value for logarithmic transformation. Default is the Euler's number.
#' @param pseudoCount The pseudo count value added to the raw counts before logarithmic transformation. Default is 1.
#'
#' @return A matrix containing marker gene evidence levels. Row names represent cell types, and column names represent genes. The higher the value, the stronger the evidence that the gene is a marker for the cell type.
#' 
#' @importFrom dplyr "%>%" n_distinct

cellTaxonomySupport <- function(markerDB,
                                species="Homo sapiens",
                                logBase=exp(1),
                                pseudoCount=1) {
    
    stopifnot(!missing(markerDB))
    
    speciesAll <- tolower(sort(unique(markerDB$species)))
    
    markerDB <- markerDB[!is.na(markerDB$species) &
                             !is.na(markerDB$Symbol) &
                             !is.na(markerDB$PMID), ]
    
    markerDB$species <- tolower(markerDB$species)
    species <- tolower(species)
    
    if("all" %in% tolower(species)) {
        species <- speciesAll
    }
    
    
    if (!all(species %in% speciesAll)) {
        stop("Species", 
             paste(species[!species %in% speciesAll], collapse = ", "), 
             " not in the marker database.",
             " To explore the available species,",
             " use the showCellTaxoSpecies() function.")
    }
    
    
    markerDB2Use <- markerDB[markerDB$species %in% species, 
                             c("species", "cell_name", 
                               "Symbol", "PMID")]
    
    markerDB2Use$cellInfo <- paste(markerDB2Use$species, 
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
