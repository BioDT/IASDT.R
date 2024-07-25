## |------------------------------------------------------------------------| #
# AllObjSizes ----
## |------------------------------------------------------------------------| #

#' Size of objects in memory
#'
#' This function calculates the size of objects in the global environment of R using [pryr::object_size] and prints a summary of objects that are greater than a specified size threshold. It is useful for memory management and identifying large objects in the workspace.
#' @param GreaterThan  numeric value specifying the size threshold in MB. Only objects larger than this value will be shown. Default is 0, which means all objects will be shown. `GreaterThan` must be a non-negative number.
#' @return The function prints a tibble containing the variables' names, their sizes in MB, and their percentage of the total size of all variables. If no objects meet the criteria, a message is printed instead. Output is sorted in descending order of the size of the objects. The function also prints the total size of all variables and the number of objects that were examined.
#' @author Ahmed El-Gabbas
#' @importFrom rlang .data
#' @export
#' @name AllObjSizes
#' @examples
#' AA1 <<- rep(1:1000, 10000)
#' AA2 <<- rep(1:1000, 100)
#'
#' AllObjSizes()
#'
#' AllObjSizes(GreaterThan = 1)
#'
#' AllObjSizes(GreaterThan = 50)

AllObjSizes <- function(GreaterThan = 0) {

  if (!is.numeric(GreaterThan) || is.na(GreaterThan) || GreaterThan < 0) {
    stop("GreaterThan must be a non-negative number")
  }

  AllVars <- ls(envir = .GlobalEnv)

  if (length(AllVars) == 0) {
    cat("No Objects are available in the global environment!\n")
  } else {
    AllVarsSize <- AllVars %>%
      sapply(
        FUN = function(x) {
          pryr::object_size(get(x)) / (1024 * 1024)
        }) %>%
      dplyr::tibble() %>%
      stats::setNames("Size") %>%
      dplyr::mutate(
        Vars = AllVars,
        Size = as.numeric(.data$Size),
        Size = round(.data$Size, 4),
        Percent = 100 * (.data$Size / sum(.data$Size)),
        Percent = round(.data$Percent, 2),
        Percent = paste0(.data$Percent, " %")) %>%
      dplyr::arrange(dplyr::desc(.data$Size)) %>%
      dplyr::select(tidyselect::all_of(c("Vars", "Size", "Percent")))


    AllVarsSize <- dplyr::filter(AllVarsSize, .data$Size >= GreaterThan)

    if (nrow(AllVarsSize) > 0) {
      cat(crayon::blue(
        "---------------------------------------------\n",
        crayon::bold(nrow(AllVarsSize)),
        " Object(s) fulfill the criteria.\n",
        "---------------------------------------------\n",
        sep = ""),
        sep = "")
      print(AllVarsSize, row.names = FALSE, n = Inf)
      cat(crayon::blue(
        "Object sizes are in MB.\n",
        "---------------------------------------------\n", sep = ""),
        sep = "")
    } else {
      cat(crayon::red(
        paste0("No variables have Size > ", GreaterThan, " MB\n")), sep = "")
    }
  }
}
