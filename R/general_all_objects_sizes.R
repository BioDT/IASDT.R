## |------------------------------------------------------------------------| #
# all_objects_sizes ----
## |------------------------------------------------------------------------| #

#' Size of objects in memory
#'
#' This function calculates the size of objects in the global environment of R
#' using [lobstr::obj_size] and prints a summary of objects that are greater
#' than a specified size threshold. It is useful for memory management and
#' identifying large objects in the workspace.
#' @param greater_than Numeric. Size threshold in MB. Only objects larger than
#'   this value will be shown. Default is 0, which means all objects will be
#'   shown. `greater_than` must be a non-negative number.
#' @param in_function Logical. This controls the scope of the function. It
#'   indicates whether the execution is done inside or outside of a function.
#'   Defaults to `FALSE` to show sizes of objects in the global environment. If
#'   set to `TRUE`, sizes of objects in the function are returned.
#' @param SizeDecimals Integer; representing the number of decimal places to
#'   show in the `Size` column. Defaults to 2.
#' @param n_objects Number of objects to show. Defaults to `Inf` meaning show
#'   all available objects.
#' @return The function prints a tibble containing the variables' names, their
#'   sizes in MB, and their percentage of the total size of all variables. If no
#'   objects meet the criteria, a message is printed instead. Output is sorted
#'   in descending order of the size of the objects. The function also prints
#'   the total size of all variables and the number of objects that were
#'   examined.
#' @author Ahmed El-Gabbas
#' @importFrom rlang .data
#' @export
#' @name all_objects_sizes
#' @examples
#' AA1 <<- rep(seq_len(1000), 10000)
#' AA2 <<- rep(seq_len(1000), 100)
#'
#' # All objects in memory
#' all_objects_sizes()
#'
#' # Objects larger than 1 MB
#' all_objects_sizes(greater_than = 1)
#'
#' # Objects larger than 50 MB
#' all_objects_sizes(greater_than = 50)
#'
#' # When called with another function, it shows the objects only available
#' # within the function
#' TestFun <- function(XX = 10) {
#'   Y <- 20
#'   C <- matrix(data = seq_len(10000), nrow = 100, ncol = 100)
#'   all_objects_sizes(in_function = TRUE)
#' }
#'
#' TestFun()
#'
#' TestFun(XX = "TEST")

all_objects_sizes <- function(
    greater_than = 0, in_function = FALSE, SizeDecimals = 2, n_objects = Inf) {

  if (in_function) {
    Environment <- parent.frame()
  } else {
    Environment <- .GlobalEnv
  }

  if (!is.numeric(greater_than) || is.na(greater_than) || greater_than < 0) {
    stop("`greater_than` must be a non-negative number", call. = FALSE)
  }

  AllVars <- ls(envir = Environment, all.names = TRUE)

  if (length(AllVars) == 0) {
    cat("No Objects are available in the global environment!\n")
  } else {

    AllVarsSize <- purrr::map_dfr(
      .x = AllVars,
      .f = ~{
        Obj <- get(.x, envir = Environment)
        Class <- paste(class(Obj), collapse = "_")

        tryCatch({
          Size <- lobstr::obj_size(Obj) / (1024 * 1024)
          Size <- round(as.numeric(Size), SizeDecimals)
          return(tibble::tibble(Object = .x, Class = Class, Size = Size))

        }, error = function(e) {
          return(tibble::tibble(Object = .x, Class = Class, Size = NA_real_))
        })
      }) %>%
      dplyr::mutate(
        Percent = round(.data$Size / sum(.data$Size, na.rm = TRUE), 2),
        Percent = scales::percent(.data$Percent, accuracy = 0.01)) %>%
      dplyr::arrange(dplyr::desc(.data$Size)) %>%
      dplyr::filter(.data$Size >= greater_than | is.na(.data$Size))

    if (nrow(AllVarsSize) > 0) {
      cat(crayon::blue(
        "---------------------------------------------------\n\t",
        crayon::bold(sum(!is.na(AllVarsSize$Size))),
        " Object(s) fulfill the criteria\n",
        "---------------------------------------------------\n",
        sep = ""),
        sep = "")

      withr::local_options(list(pillar.sigfig = 4))
      print(AllVarsSize, n = n_objects)

      if (sum(is.na(AllVarsSize$Size)) > 0) {
        NA_Var <- dplyr::filter(AllVarsSize, is.na(.data$Size)) %>%
          dplyr::pull(.data$Object) %>%
          paste(collapse = " | ")

        cat(crayon::blue(
          paste0(
            "`lobstr::obj_size` was not able to get the object ",
            "size of the following object(s): ", NA_Var, "\n"), sep = ""),
          sep = "")
      }

      cat(crayon::blue(
        "Object sizes are in MB.\n",
        "---------------------------------------------------\n", sep = ""),
        sep = "")
    } else {
      cat(crayon::red(
        paste0("No object has Size > ", greater_than, " MB\n")), sep = "")
    }
  }
}
