#' Save Cell Type Results
#'
#' This function saves the cell type results to an Excel file.
#'
#' @param cellTypeRes A list containing the annotation results for each cell cluster.
#' @param path The path where the Excel file will be saved. For example, path = "directory/results.xlsx".
#'
#' @return None
#'
#' @import openxlsx

saveCellTypeRes <- function(cellTypeRes, path) {
    openxlsx::write.xlsx(x = cellTypeRes,
                         file = path,
                         rowNames = TRUE)
}
