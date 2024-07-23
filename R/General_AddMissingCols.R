## |------------------------------------------------------------------------| #
# AddMissingCols ----
## |------------------------------------------------------------------------| #

#' Add missing columns to a data frame with specified fill values
#'
#' This function checks a data frame for missing columns specified by the user. If any are missing, it adds these columns to the data frame, filling them with a specified value.
#'
#' @name AddMissingCols
#' @author Ahmed El-Gabbas
#' @return NULL
#' @param DT A data frame to which missing columns will be added. This parameter cannot be NULL.
#' @param FillVal The value to fill the missing columns with. This parameter defaults to `NA_character_`, but can be changed to any scalar value as required.
#' @param ... Column names as character strings.
#' @return a data frame with the missing columns added, if any were missing.
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

  if (is.null(DT) || is.null(FillVal)) {
    stop("DT can not be NULL")
  }

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
    DT <- tibble::add_column(DT, !!!Add_DF) %>%
      tibble::tibble()
  }
  return(DT)
}
