## |------------------------------------------------------------------------| #
# RecordArgs ----
## |------------------------------------------------------------------------| #

#' Capture and record function arguments
#'
#' `RecordArgs()` is a utility function that captures and records both the
#' unevaluated and evaluated forms of arguments passed to a parent function. It
#' returns a tibble with columns reflecting argument states: when unevaluated
#' and evaluated values differ, columns are named with `_orig` and `_eval`
#' suffixes; when they are the same (including symbols evaluated to scalars), a
#' single column is used with the argument name alone. The function dynamically
#' handles scalars, call objects, and complex objects (e.g., `lm` models,
#' `SpatRaster` objects), preserving their structure appropriately.
#'
#' @param ExportPath Character. The path to an `.RData` file where the output
#'   tibble will be exported. If `NULL` (default), the tibble is returned
#'   without saving. If provided, the tibble is saved to the specified file and
#'   `NULL` is returned invisibly.
#'
#' @param call Language object (optional). The call to the parent function, as
#'   provided by `match.call()` from the caller. If `NULL` (default), the
#'   function falls back to `sys.call(-1)` to capture the parent call. Used to
#'   ensure accurate argument capture in iterative contexts (e.g., `lapply`,
#'   `purrr::map`).
#' @param env Environment (optional). The environment in which to evaluate the
#'   arguments, typically provided by `parent.frame()` from the caller. If
#'   `NULL` (default), the function uses `parent.frame()` to determine the
#'   evaluation environment. Used to resolve variables in iterative contexts.
#'
#' @details This function must be called within another function to capture and
#'   record the parent function’s arguments. It uses either a provided `call`
#'   (e.g., from `match.call()`) or `sys.call(-1)` to retrieve the call,
#'   evaluates arguments in the specified or default parent environment, and
#'   merges them with the parent function’s formal argument defaults. The
#'   handling of argument types is as follows:
#'   - `NULL` values (e.g., from defaults or arguments) are converted to the
#'   string `"NULL"` in the output tibble for compatibility.
#'   - Symbols (e.g., `i` in `lapply`) match their evaluated scalar values,
#'   appearing as single columns.
#'   - Calls (e.g., `a + b`) are preserved as strings via `deparse()` in
#'   `_orig` columns and evaluated in `_eval` columns.
#'   - Multi-element vectors (e.g., `c(1, 2)`) appear as single columns if
#'   unevaluated and evaluated forms match, otherwise as `_orig`/`_eval` pairs.
#'   - Complex objects (e.g., `lm`, `SpatRaster`) are wrapped in list columns
#'   in `_eval`.
#'
#'   The output tibble columns are ordered by the original argument sequence,
#'   with single columns (matching values) first, followed by `_orig` and
#'   `_eval` pairs.
#'
#' @return A `tibble` containing the unevaluated and evaluated forms of the
#'   parent function’s arguments. Column naming depends on whether unevaluated
#'   and evaluated values match:
#'   - **Single columns** (e.g., `y`): Used when unevaluated and evaluated
#'   values are identical or effectively equivalent (e.g., scalars like `2`,
#'   defaults like `10`, or symbols like `i` evaluating to `1` in loops),
#'   containing the evaluated value as-is; `NULL` is recorded as `"NULL"`.
#'   - **Paired columns** (e.g., `x_orig`, `x_eval`):
#'     - `*_orig`: Unevaluated expressions as character strings (e.g.,
#'   `"a + b"`) or scalars as-is for non-call objects; `NULL` is recorded as
#'   `"NULL"`.
#'     - `*_eval`: Evaluated values, either scalars (e.g., `8`) or list columns
#'   for complex objects (e.g., `lm`, `SpatRaster`); `NULL` is recorded as
#'   `"NULL"`.
#'
#'   If `ExportPath` is `NULL` (default), the tibble is returned. If provided,
#'   the tibble is saved to the specified `.RData` file and `NULL` is returned
#'   invisibly.
#'
#' @author Ahmed El-Gabbas
#' @export
#' @name RecordArgs
#' @examples
#' a <- 5
#' b <- 3
#'
#' Function1 <- function(w = 5, x, y, z = c(1, 2)) {
#'   Args <- IASDT.R::RecordArgs(call = match.call(), env = parent.frame())
#'   return(Args)
#' }
#'
#' # Basic usage with scalars and expressions
#' Out1 <- Function1(x = a + b, y = 2)
#'
#' Out1
#'
#' Out1$w              # 5 (single column, default matches evaluated)
#' Out1$x_orig         # "a + b" (unevaluated expression)
#' Out1$x_eval         # 8 (evaluated result)
#' Out1$y              # 2 (single column, scalar matches evaluated)
#' Out1$z_orig         # "c(1, 2)"
#' Out1$z_eval         # c(1, 2) (single column, default matches evaluated)
#'
#' # --------------------------------------------------------------
#'
#' # Usage with NULL, scalars, and complex objects (lm and Raster)
#' Out2 <- Function1(
#'   w = NULL,
#'   x = a + b,
#'   y = stats::lm(mpg ~ disp + hp, data = mtcars),
#'   z = raster::raster())
#'
#' Out2
#'
#' Out2$w              # "NULL" (recorded as character string)
#' Out2$x_orig         # "a + b" (unevaluated call)
#' Out2$x_eval         # 8 (evaluated scalar)
#' Out2$y_orig         # "stats::lm(mpg ~ disp + hp, data = mtcars)" (call)
#' Out2$y_eval[[1]]    # <lm object> (evaluated complex object)
#' Out2$z_orig         # "raster::raster()" (call)
#' Out2$z_eval[[1]]    # <RasterLayer object> (evaluated complex object)
#'
#' # --------------------------------------------------------------
#'
#' # Usage with purrr::pmap for multiple inputs
#' w_values <- 1:3
#' x_values <- c(a + b, 10, 15)
#' y_values <- c("ABCD", "XYZ123", "TEST")
#' Out3 <- purrr::pmap(
#'   .l = list(w = w_values, x = x_values, y = y_values),
#'   .f = function(w, x, y) {
#'     Function1(
#'       w = w,
#'       x = x,
#'       y = stringr::str_extract(y, "B.+$"),
#'       z = terra::rast(system.file("ex/elev.tif", package="terra")))
#'   }) %>%
#'   dplyr::bind_rows()
#'
#' Out3
#'
#' Out3$w        # 1, 2, 3
#' Out3$x        # 8, 10, 15
#' Out3$y_orig   # 'stringr::str_extract(y, "B.+$")', repeated for each row
#' Out3$y_eval   # "BCD", NA, NA
#' Out3$z_orig   # "terra::rast(...))", repeated for each row
#' Out3$z_eval   # Packed SpatRaster, repeated for each row

RecordArgs <- function(ExportPath = NULL, call = NULL, env = NULL) {

  # Capture the parent function's call: use provided call (e.g., from
  # match.call()) or fall back to sys.call(-1) for direct calls
  call_info <- if (!is.null(call)) {
    call
  }  else {
    sys.call(-1)
  }

  # Check if call_info is valid; stop if not called within a function
  if (is.null(call_info)) {
    stop(
      "RecordArgs() must be called from within another function", call. = FALSE)
  }

  # Get the name of the calling function as a character string
  calling_func <- deparse(call_info[[1]])

  # Determine the environment for evaluation: use provided env or default to
  # parent.frame()
  parent_env <- env %||% parent.frame()

  # Retrieve the parent function and its formal arguments (including defaults)
  parent_func <- sys.function(-1)
  formals_full <- formals(parent_func)

  # Replace NULL values with "NULL" for tibble compatibility
  ReplaceNULL <- function(x) {
    purrr::map(x, ~ if (is.null(.x)) "NULL" else .x) %>%
      stats::setNames(names(x))
  }


  # Extract the arguments from the call, excluding the function name (first
  # element), and preserve defaults for missing arguments
  args_list <- as.list(call_info)[-1]
  args_list <- utils::modifyList(formals_full, args_list) %>%
    # Replace NULL values with "NULL" to avoid errors
    ReplaceNULL()

  # Evaluate the captured arguments in the parent environment, handling all
  # cases safely
  args_values <- lapply(args_list, function(x) {
    tryCatch(eval(x, envir = parent_env), error = function(e) {
      # If evaluation fails in parent_env, try globalenv() to resolve globals
      # (e.g., 'a' and 'b' in 'a + b')
      tryCatch(eval(x, envir = globalenv()), error = function(e2) NA)
    })
  })

  # Name the evaluated values with their corresponding argument names
  recorded_values <- stats::setNames(args_values, names(args_list))

  # Merge evaluated values with defaults, overriding defaults with provided
  # values
  Evaluated <- utils::modifyList(formals_full, recorded_values) %>%
    # Replace NULL values with "NULL" to avoid errors
    ReplaceNULL()

  # Merge unevaluated expressions with defaults, keeping unevaluated forms
  Unevaluated <- utils::modifyList(formals_full, args_list) %>%
    # Replace NULL values with "NULL" to avoid errors
    ReplaceNULL()

  # Get the argument names in their original order from the function definition
  arg_names <- names(formals_full)

  # Identify which arguments have identical unevaluated and evaluated values
  # - Calls (e.g., a + b) are always different
  # - Symbols (e.g., i in loops) match their evaluated scalar
  # - Scalars and defaults are compared directly
  same_values <- purrr::map2_lgl(
    # Coerce pairlist to list for purrr compatibility
    .x = as.list(Unevaluated),
    .y = as.list(Evaluated),
    .f = function(u, e) {
      # Calls always differ from evaluated
      if (is.call(u)) return(FALSE)
      # Symbols match their evaluated value
      if (is.symbol(u)) return(identical(u, e))
      # Direct comparison for scalars and defaults
      identical(u, e)
    })

  # Define column names: single columns for matching values, pairs for differing
  # ones
  #
  # Names for arguments with identical values
  single_cols <- arg_names[same_values]
  # Names for arguments needing _orig/_eval pairs
  diff_cols <- arg_names[!same_values]

  # Construct column names for the tibble

  # e.g., "y", "z" (unchanged)
  single_cols <- single_cols
  # e.g., "x_eval" (evaluated values)
  eval_cols <- paste0(diff_cols, "_eval")
  # e.g., "x_orig" (unevaluated forms)
  uneval_cols <- paste0(diff_cols, "_orig")

  # Format evaluated values: scalars as-is (e.g., 5, 8, NA), multi-element
  # vectors (e.g., c(1, 2)) and complex objects in lists
  eval_values <- purrr::map(
    .x = as.list(Evaluated),
    .f = function(x) {
      if (is.vector(x) && !is.list(x)) {
        if (length(x) == 1) {
          # Scalars (e.g., 5, 8, NA) remain as-is
          return(x)
        } else {
          # Multi-element vectors (e.g., c(1, 2)) wrapped in a list
          return(list(x))
        }
      } else if (is.language(x)) {
        # Evaluate language objects and handle based on length
        eval_result <- eval(x, envir = parent_env)
        if (is.vector(eval_result) && !is.list(eval_result) &&
            length(eval_result) == 1) {
          return(eval_result)
        } else {
          return(list(eval_result))
        }
      } else {
        if (inherits(x, "SpatRaster")) {
          # Wrap SpatRaster objects for storage
          x <- terra::wrap(x)
        }
        # Wrap complex objects (e.g., lm, RasterLayer) in a list
        return(list(x))
      }
    })

  # Format unevaluated values: calls as strings (e.g., "a + b"), scalars as-is
  # (e.g., 5), multi-element vectors (e.g., c(1, 2)) in lists
  uneval_values <- purrr::map(
    .x = as.list(Unevaluated),
    .f = function(x) {
      if (is.call(x)) {
        # Convert calls (e.g., a + b) to strings without quotes
        return(noquote(deparse(x)))
      } else if (is.vector(x) && !is.list(x)) {
        if (length(x) == 1) {
          # Scalars (e.g., 5) remain as-is
          return(x)
        } else {
          # Multi-element vectors (e.g., c(1, 2)) wrapped in a list
          return(list(x))
        }
      } else {
        if (inherits(x, "SpatRaster")) {
          # Wrap SpatRaster objects
          x <- terra::wrap(x)
        }
        # Wrap complex objects in a list
        return(list(x))
      }
    })

  # Combine all values into a named list for tibble construction
  tibble_data <- c(
    # Unevaluated differing values, only if diff_cols is non-empty
    if (length(diff_cols) > 0) {
      stats::setNames(uneval_values[diff_cols], uneval_cols)
    },
    # Evaluated differing values, only if diff_cols is non-empty
    if (length(diff_cols) > 0) {
      stats::setNames(eval_values[diff_cols], eval_cols)
    },
    # Single columns for matching values
    stats::setNames(eval_values[single_cols], single_cols))

  # Create the tibble from the combined data
  result <- tibble::as_tibble(tibble_data)

  # Reorder the tibble columns according to the desired order
  # Define the desired column order
  desired_order <- purrr::map(
    arg_names, ~ c(.x, paste0(.x, "_orig"), paste0(.x, "_eval"))) %>%
    unlist()
  result <- dplyr::select(result, tidyselect::any_of(desired_order))

  # Return the tibble or save it to a file based on ExportPath
  if (is.null(ExportPath)) {

    # Return the tibble if no export path is provided
    return(result)

  } else {

    # Extract the calling function name without the namespace
    calling_func2 <- stringr::str_remove(calling_func, "^.+::")

    # Save to .RData file if ExportPath is specified
    IASDT.R::SaveAs(
      InObj = result, OutObj = paste0("Args_", calling_func2),
      OutPath = ExportPath)
    # Return NULL invisibly after saving
    return(invisible(NULL))
  }
}
