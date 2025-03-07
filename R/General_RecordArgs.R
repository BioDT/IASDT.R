#' Capture and record function arguments
#'
#' `RecordArgs()` is a utility function designed to capture and record both the
#' evaluated and unevaluated forms of arguments passed to a parent function. It
#' returns a `tibble` where each argument is represented by two columns: one for
#' the unevaluated expression (suffix `_orig`) and one for the evaluated value
#' (suffix `_eval`). The function dynamically handles scalars, call objects, and
#' complex objects (e.g., `lm` models, `SpatRaster` objects), preserving their
#' structure appropriately.
#' @return A tibble containing the unevaluated and evaluated forms of the parent
#'   function’s arguments, with columns named using the argument names followed
#'   by `_orig` and `_eval` suffixes. if `ExportPath` is `NULL` (default), the
#'   tibble is returned; otherwise, it is saved as an `.RData` file and returns
#'   `NULL`.
#' @author Ahmed El-Gabbas
#' @name RecordArgs
#' @export
#' @param ExportPath Character. The path of an `RData` file to export the tibble
#'   to . If `NULL` (default), the tibble will not be saved to disk, but will be
#'   returned as an object.
#' @details This function must be called from within another function. It uses
#'   `sys.call(-1)` to capture the parent function’s call, evaluates arguments
#'   in the parent environment, and combines them with default values from the
#'   parent function’s formal arguments. Unevaluated expressions (e.g., `a + b`)
#'   are converted to character strings, while evaluated values are kept as-is
#'   for scalars or wrapped in lists for complex objects. Columns are ordered
#'   based on the original argument sequence in the parent function’s
#'   definition, with unevaluated variants preceding evaluated variants (e.g.,
#'   `x_orig`, `x_eval`).
#' @examples
#' # Basic usage with scalars and expressions
#' a <- 5
#' b <- 3
#' Fun1 <- function(w = 5, x, y, z = 10) {
#'   Args <- RecordArgs()
#'   return(Args)
#' }
#' AA <- Fun1(x = a + b, y = 2)
#' AA$w_orig    # 5
#' AA$w_eval    # 5
#' AA$x_orig    # "a + b"
#' AA$x_eval    # 8
#' AA$y_orig    # 2
#' AA$y_eval    # 2
#' AA$z_orig    # 10
#' AA$z_eval    # 10
#'
#' # Usage with complex objects (lm and SpatRaster)
#' AA <- Fun1(
#'   w = 10,
#'   x = a + b,
#'   y = stats::lm(mpg ~ disp + hp, data = mtcars),
#'   z = terra::rast(system.file("ex/logo.tif", package = "terra")))
#' AA$w_orig       # 10
#' AA$w_eval       # 10
#' AA$x_orig       # "a + b"
#' AA$x_eval       # 8
#' AA$y_orig       # "lm(mpg ~ disp + hp, data = mtcars)"
#' AA$y_eval[[1]]  # lm object
#' AA$z_orig    # "terra::rast(system.file("ex/logo.tif", package = "terra"))"
#' AA$z_eval[[1]]  # SpatRaster object

RecordArgs <- function(ExportPath = NULL) {

  # Get the call to the parent function (one level up)
  call_info <- sys.call(-1)

  if (is.null(call_info)) {
    stop(
      "RecordArgs() must be called from within another function",
      call. = FALSE)
  }

  # Extract the arguments, excluding the function name
  args_list <- as.list(call_info)[-1]

  # Extract the name of the calling function
  calling_func <- deparse(call_info[[1]])

  # Get the parent function's environment and formal arguments
  parent_env <- parent.frame()
  parent_func <- sys.function(-1)
  formals_full <- formals(parent_func)

  # Initialize results
  Evaluated <- NULL
  Unevaluated <- NULL

  # Evaluate the arguments in the parent environment
  args_values <- lapply(args_list, eval, envir = parent_env)
  recorded_values <- stats::setNames(args_values, names(args_list))
  # Combine with default values
  Evaluated <- modifyList(formals_full, recorded_values)
  # Store unevaluated expressions
  Unevaluated <- modifyList(formals_full, args_list)

  # Get argument names in their original order - use formal argument order
  arg_names <- names(formals_full)

  # Construct column names
  eval_cols <- paste0(arg_names, "_eval")
  uneval_cols <- paste0(arg_names, "_orig")

  # Prepare Evaluated values: keep scalars as-is, wrap complex objects in lists
  eval_values <- purrr::map(
    .x = as.list(Evaluated),
    .f = function(x) {
      if (is.vector(x) && length(x) == 1 && !is.list(x)) {
        # Scalars (numeric, character) stay as-is
        return(x)
      } else {
        # Complex objects (e.g., lm, SpatRaster) wrapped in list

        if (inherits(x, "SpatRaster")) {
          x <- terra::wrap(x)
        }

        return(list(x))
      }
    })

  # Prepare Unevaluated values: keep calls as-is, wrap in list, scalars stay
  # as-is
  uneval_values <- purrr::map(
    .x = as.list(Unevaluated),
    .f = function(x) {
      if (is.call(x)) {
        # Keep call objects as-is, wrapped in list
        return(list(x))
      } else if (is.vector(x) && length(x) == 1 && !is.list(x)) {
        # Scalars stay as-is
        return(x)
      } else {
        # Other objects wrapped in list
        if (inherits(x, "SpatRaster")) {
          x <- terra::wrap(x)
        }
        return(list(x))
      }
    })

  # Combine into a tibble
  tibble_data <- c(
    # Unevaluated first
    stats::setNames(uneval_values, uneval_cols),
    # Evaluated second
    stats::setNames(eval_values, eval_cols))

  # Create tibble and reorder columns by argument order
  result <- tibble::as_tibble(tibble_data)

  # Define the desired column order: unevaluated before evaluated for each arg
  desired_order <- lapply(
    arg_names,
    function(n) {
      c(paste0(n, "_orig"), paste0(n, "_eval"))
    }) %>%
    unlist()

  result <- dplyr::select(result, tidyselect::all_of(desired_order))

  # Export the tibble to a file if a path is provided
  if (is.null(ExportPath)) {
    return(result)
  } else {
    IASDT.R::SaveAs(
      InObj = result, OutObj = paste0(calling_func, "_args"),
      OutPath = ExportPath)
    return(invisible(NULL))
  }


}
