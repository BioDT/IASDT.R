## |------------------------------------------------------------------------| #
# trim_hmsc ----
## |------------------------------------------------------------------------| #

#' Trim an Hmsc Model Object by Removing Specified Components
#'
#' Removes specified components from an `Hmsc` model object one at a time. This
#' is used to keep only a smaller `Hmsc` object, containing only needed list
#' items. This is useful in situations in which the fitted Hmsc model is large
#' (e.g. in GBs) and downstream computations are implemented on parallel.
#'
#' @param model An object of class `Hmsc`, containing fitted Hmsc model. Must
#'   not be `NULL`.
#' @param names_to_remove A character vector specifying the names of components
#'   to remove from the model (e.g., `"postList"`, `"Y"`). If `NULL`, no
#'   trimming is implemented. Must be non-empty and match names in the model.
#' @return An `Hmsc` object with the specified components removed.
#' @details This function is used to reduce the memory footprint of `Hmsc`
#'   models by removing unnecessary components before further processing. It
#'   converts the model to a plain list to avoid S3 method overhead, removes
#'   components iteratively, and restores the `Hmsc` class. The iterative
#'   approach is slower than vectorized subsetting but may be preferred for
#'   specific use cases requiring step-by-step removal. The simple trimming of
#'   list items; e.g. `model[c("postList", X")] <- NULL` took comparably too
#'   much time for trimming done inside functions for jobs submitted via SLURM
#'
#' @author Ahmed El-Gabbas
#' @export
#' @examples
#' require(Hmsc)
#'
#' (model <- Hmsc::TD$m)
#'
#' (trimmed_model <- trim_hmsc(
#'   model, names_to_remove = c("postList", "rL", "ranLevels")))
#'
#' setdiff(names(model), names(trimmed_model))
#'
#' lobstr::obj_size(model)
#' lobstr::obj_size(trimmed_model)
#'

trim_hmsc <- function(model, names_to_remove = NULL) {

  # Input validation
  if (is.null(model)) {
    ecokit::stop_ctx("model cannot be NULL", include_backtrace = TRUE)
  }

  if (!inherits(model, "Hmsc")) {
    ecokit::stop_ctx(
      "model must be an object of class Hmsc",
      class_model = class(model), include_backtrace = TRUE)
  }

  if (!is.character(names_to_remove) || length(names_to_remove) == 0L) {
    ecokit::stop_ctx(
      "names_to_remove must be a non-empty character vector",
      class_names_to_remove = class(names_to_remove),
      include_backtrace = TRUE)
  }

  names_to_remove <- unique(names_to_remove)

  names_difference <- setdiff(names_to_remove, names(model))
  if (length(names_difference) > 0L) {
    ecokit::stop_ctx(
      "Some names_to_remove are not present in the model object",
      names_difference = names_difference, include_backtrace = TRUE)
  }

  # Convert to plain list to avoid S3 method overhead
  model <- as.list(model)

  # Remove one item at a time with for loop
  for (name in names_to_remove) {
    model[[name]] <- NULL
  }

  # Restore Hmsc class
  class(model) <- "Hmsc"

  invisible(gc())
  return(model)
}
