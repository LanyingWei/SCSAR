#' Show CellMarker Database Tissue Information
#'
#' This function displays a searchable table with the tissue information and the corresponding species from the CellMarker 2.0 database.
#'
#' @return an HTML widget to display tabular data.
#' @export
#'
#' @examples
#' showCellMarkerTissue()
#'
#' @importFrom DT datatable
#'
showCellMarkerTissue <- function() {
    DB <- cellMarkerDB
    DB <- DB[, c("species", "tissue_class", "tissue_type")]
    DB <- unique(DB)
    rownames(DB) <- seq(length = nrow(DB))
    resTable <- DT::datatable(DB)
    return(resTable)
}

#' Show Cell Taxonomy Database Species Information
#'
#' This function displays a searchable table with the species information and the number of corresponding cell types from the Cell Taxonomy database.
#' 
#' @return an HTML widget to display tabular data.
#' @export
#'
#' @examples
#' showCellTaxoSpecies()
#'
#' @importFrom DT datatable
#' @importFrom dplyr group_by summarise
#'
showCellTaxoSpecies <- function() {
    DB <- cellTaxoDB
    DB <- DB[, c("species", "cell_name")]
    DB <- DB %>% 
        group_by(species) %>% 
        summarise(count = n_distinct(cell_name))
    colnames(DB)[2] <- "number_of_cell_types"
    DB <- DB[order(DB$number_of_cell_types, decreasing = TRUE), ]
    DB <- as.data.frame(DB)
    resTable <- DT::datatable(DB)
    return(resTable)
}
