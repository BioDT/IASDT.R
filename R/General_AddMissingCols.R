## |------------------------------------------------------------------------| #
# AddMissingCols ----
## |------------------------------------------------------------------------| #

#' Add missing columns to data frame
#'
#' Add missing columns to data frame
#' @name AddMissingCols
#' @author Ahmed El-Gabbas
#' @return NULL
#' @param DT data frame
#' @param FillVal value to be used: default: NA_character
#' @param ... list of column names to add
#' @export
#' @examples
#' mtcars %>%
#'  dplyr::select(1:3) %>%
#'  AddMissingCols(FillVal = NA_character_, A, B, C) %>%
#'  AddMissingCols(FillVal = as.integer(10), D)
#'
#' # -------------------------------------------
#'
#' AddCols <- c("Add1", "Add2")
#' mtcars %>%
#'  dplyr::select(1:3) %>%
#'  AddMissingCols(FillVal = NA_real_, AddCols)

AddMissingCols <- function(DT, FillVal = NA_character_, ...) {
  Cols <- as.character(rlang::ensyms(...))

  if (any(Cols %in% ls(envir = parent.env(rlang::caller_env())))) {
    Cols <- get(Cols, envir = parent.env(rlang::caller_env()))
  }

  Cols2Add <- setdiff(Cols, names(DT))

  Add_DF <- rep(FillVal, length(Cols2Add)) %>%
    matrix(nrow = 1) %>%
    as.data.frame() %>%
    stats::setNames(Cols2Add) %>%
    tibble::as_tibble()

  if (length(Cols2Add) != 0) {
    DT <- tibble::add_column(DT, !!!Add_DF)
  }
  return(tibble::tibble(DT))
}
