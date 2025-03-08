## |------------------------------------------------------------------------| #
# RecordArgs ----
## |------------------------------------------------------------------------| #

#' Capture and record function arguments
#'
#' `RecordArgs()` is a utility function that captures and records both the
#' unevaluated and evaluated forms of arguments passed to a parent function. It
#' returns a tibble with columns reflecting argument states: when unevaluated
#' and evaluated values differ, columns are named with `_orig` and `_eval`
#' suffixes; when they are the same, a single column is used with the argument
#' name alone. The function dynamically handles scalars, call objects, and
#' complex objects (e.g., `lm` models, `SpatRaster` objects), preserving their
#' structure appropriately.
#'
#' @param ExportPath Character. The path to an `.RData` file where the tibble
#'   will be exported. If `NULL` (default), the tibble is returned without
#'   saving. If provided, the tibble is saved to the specified file and `NULL`
#'   is returned invisibly.
#'
#' @details This function must be called from within another function. It uses
#'   `sys.call(-1)` to capture the parent function’s call, evaluates arguments
#'   in the parent environment, and combines them with default values from the
#'   parent function’s formal arguments. Unevaluated expressions (e.g., `a + b`)
#'   are preserved as `call` objects, while evaluated values are kept as-is for
#'   scalars or wrapped in lists for complex objects. Columns are ordered based
#'   on the original argument sequence: single columns (for matching values)
#'   appear first, followed by `_orig` and `_eval` pairs (for differing values)
#'   in that order.
#'
#' @return A `tibble` containing the unevaluated and evaluated forms of the
#'   parent function’s arguments. Column naming depends on whether unevaluated
#'   and evaluated values match:
#'   - **Single columns** (e.g., `y`): Used when unevaluated and evaluated
#'   values are identical (e.g., scalars like `2` or defaults like `10`),
#'   containing the evaluated value as-is.
#'   - **Paired columns** (e.g., `x_orig`, `x_eval`):
#'     - `*_orig`: Unevaluated expressions as list columns with `call` objects
#'   (e.g., `a + b`) or scalars as-is.
#'     - `*_eval`: Evaluated values, either scalars (e.g., `8`) or list columns
#'   for complex objects (e.g., `lm`, `SpatRaster`).
#'
#'   If `ExportPath` is `NULL` (default), the tibble is returned. If provided,
#'   the tibble is saved to the specified `.RData` file and `NULL` is returned
#'   invisibly.
#'
#' @author Ahmed El-Gabbas
#' @export
#' @examples
#' a <- 5
#' b <- 3
#' Function1 <- function(w = 5, x, y, z = 10) {
#'   Args <- RecordArgs()
#'   return(Args)
#' }
#'
#' # --------------------------------------------------------------
#'
#' # Basic usage with scalars and expressions
#' Out1 <- Function1(x = a + b, y = 2)
#' Out1$w              # 5 (single column, same as orig and eval)
#' Out1$x_orig         # call object: a + b
#' Out1$x_eval         # 8
#' Out1$y              # 2 (single column)
#' Out1$z              # 10 (single column)
#'
#' # --------------------------------------------------------------
#'
#' #' # Usage with complex objects (lm and SpatRaster)
#' Out2 <- Function1(
#'   w = 10,
#'   x = a + b,
#'   y = stats::lm(mpg ~ disp + hp, data = mtcars),
#'   z = terra::rast(system.file("ex/logo.tif", package = "terra")))
#' Out2$w              # 10 (single column)
#' Out2$x_orig         # call object: a + b
#' Out2$x_eval         # 8
#' Out2$y_orig         # call object: lm(mpg ~ disp + hp, data = mtcars)
#' Out2$y_eval[[1]]    # lm object
#' Out2$z_orig         # call object: terra::rast(system.file(...))
#' Out2$z_eval[[1]]    # SpatRaster object
#'

RecordArgs <- function(ExportPath = NULL) {
  # Get the call to the parent function (one level up)
  call_info <- sys.call(-1)

  if (is.null(call_info)) {
    stop(
      "RecordArgs() must be called from within another function", call. = FALSE)
  }

  # Extract the arguments, excluding the function name
  args_list <- as.list(call_info)[-1]

  # Extract the name of the calling function
  calling_func <- deparse(call_info[[1]])

  # Get the parent function's environment and formal arguments
  parent_env <- parent.frame()
  parent_func <- sys.function(-1)
  formals_full <- formals(parent_func)

  # Evaluate the arguments in the parent environment
  args_values <- lapply(args_list, eval, envir = parent_env)
  recorded_values <- stats::setNames(args_values, names(args_list))
  # Combine with default values
  Evaluated <- utils::modifyList(formals_full, recorded_values)
  # Store unevaluated expressions
  Unevaluated <- utils::modifyList(formals_full, args_list)

  # Get argument names in their original order
  arg_names <- names(formals_full)

  # Determine which arguments have identical unevaluated and evaluated values
  same_values <- purrr::map2_lgl(
    .x = Unevaluated,
    .y = Evaluated,
    .f = function(u, e) {
      if (is.call(u)) return(FALSE)  # Calls always differ from evaluated
      identical(u, e)
    })

  # Names for single columns
  single_cols <- arg_names[same_values]
  # Names for orig/eval pairs
  diff_cols <- arg_names[!same_values]

  # Construct column names
  single_cols <- single_cols                      # e.g., "y", "z"
  eval_cols <- paste0(diff_cols, "_eval")         # e.g., "x_eval"
  uneval_cols <- paste0(diff_cols, "_orig")       # e.g., "x_orig"

  # Prepare Evaluated values: keep scalars as-is, wrap complex objects in lists
  eval_values <- purrr::map(
    .x = as.list(Evaluated),
    .f = function(x) {
      if (is.vector(x) && length(x) == 1 && !is.list(x)) {
        return(x)
      } else {
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
        return(list(x))
      } else if (is.vector(x) && length(x) == 1 && !is.list(x)) {
        return(x)
      } else {
        if (inherits(x, "SpatRaster")) {
          x <- terra::wrap(x)
        }
        return(list(x))
      }
    })

  # Combine into a tibble
  tibble_data <- c(
    stats::setNames(uneval_values[diff_cols], uneval_cols),
    stats::setNames(eval_values[diff_cols], eval_cols),
    stats::setNames(eval_values[single_cols], single_cols)
  )

  # Create tibble and reorder columns by argument order
  result <- tibble::as_tibble(tibble_data)

  # Define the desired column order: single columns first, then orig/eval pairs
  desired_order <- c(
    single_cols,
    unlist(lapply(
      diff_cols,
      function(n) {
        c(paste0(n, "_orig"), paste0(n, "_eval"))
      }
    ))
  )

  result <- dplyr::select(result, tidyselect::all_of(desired_order))

  # Export the tibble to a file if a path is provided
  if (is.null(ExportPath)) {
    return(result)
  } else {
    IASDT.R::SaveAs(
      InObj = result, OutObj = paste0("Args_", calling_func),
      OutPath = ExportPath)
    return(invisible(NULL))
  }
}
