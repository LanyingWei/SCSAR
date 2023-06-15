#' Save Cell Type Annotation Results
#'
#' This function saves the cell type annotation results for each cluster to an Excel file.
#'
#' @param cellTypeRes A list containing the annotation results for each cell cluster.
#' @param path The path where the Excel file will be saved.
#'
#' @return None
#'
#' @import openxlsx

saveCellTypeRes <- function(cellTypeRes, path) {
    outFile <- openxlsx::createWorkbook(path)
    for (i in seq_along(cellTypeRes)) {
        openxlsx::addWorksheet(outFile, tolower(names(cellTypeRes)[i]))
        openxlsx::write.xlsx(outFile, 
                             sheet = i,
                             x = cellTypeRes[[i]], 
                             rowNames = TRUE) 
    }
    openxlsx::saveWorkbook(outFile)
}
