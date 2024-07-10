## |------------------------------------------------------------------------| #
# AllObjSizes ----
## |------------------------------------------------------------------------| #

#' Size of objects in memory
#'
#' Size of objects in memory. This function print the size allocated by any of loaded objects into `R`
#'
#' @param GreaterThan Only show objects with size > specified value in MB. Default: 0, which means show all variables size
#' @author Ahmed El-Gabbas
#' @importFrom rlang .data
#' @export
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

    if (!is.na(GreaterThan)) {
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
          paste0("No variables have Size > ",
                 GreaterThan, " MB\n")), sep = "")
      }
    } else {
      print(AllVarsSize, row.names = FALSE, n = Inf)
      cat(crayon::green("Object sizes are in MB"), sep = "")
    }
  }
}
