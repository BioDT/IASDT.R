## |------------------------------------------------------------------------| #
# save_as ----
## |------------------------------------------------------------------------| #

#' Save an object to a file with a new name
#'
#' This function saves an R object to a specified file path with a potentially
#' new name. It is useful for renaming objects during the save process. The
#' function supports saving objects in `RData`, `qs2`, `feather`, and `rds`
#' formats. The format is determined by the extension of the file path.
#'
#' @param object The input object to be saved. This can be an actual R object or
#'   a character string representing the name of an object.
#' @param object_name Character. The new name for the saved `RData` object. This
#'   name is used when the object is loaded back into R. Default is `NULL`. This
#'   is required when saving `RData` files.
#' @param out_path Character. File path (ends with either `*.RData`, `*.qs2`,
#'   `feather`, and `rds`) where the object be saved. This includes the
#'   directory and the file name.
#' @param n_threads Character. Number of threads to use when compressing data.
#'   See [qs2::qs_save].
#' @param feather_compression Character. The compression algorithm to use when
#'   saving the object in the `feather` format. The default is "zstd". See
#'   [arrow::write_feather].
#' @param ... Additional arguments to be passed to the respective save
#'   functions. [base::save] for `RData` files; [qs2::qs_save] for `qs2` files;
#'   [arrow::write_feather] for `feather` files; and [base::saveRDS] for `rds`
#'   files.
#' @name save_as
#' @author Ahmed El-Gabbas
#' @return The function does not return a value but saves an object to the
#'   specified file path.
#' @export
#' @examples
#' TMP_Folder <- IASDT.R::path(tempdir(), stringi::stri_rand_strings(1, 5))
#' fs::dir_create(TMP_Folder)
#' list.files(TMP_Folder)
#'
#' # save iris data in `iris2.RData` with `iris2` object name
#' save_as(iris, "iris2", IASDT.R::path(TMP_Folder, "iris2.RData"))
#' list.files(TMP_Folder, pattern = "^.+.RData")
#'
#' (load(IASDT.R::path(TMP_Folder, "iris2.RData")))
#'
#' tibble::tibble(iris2)

save_as <- function(
    object = NULL, object_name = NULL, out_path = NULL, n_threads = 5L,
    feather_compression = "zstd", ...) {

  if (is.null(object) || is.null(out_path)) {
    stop("`object` and `out_path` cannot be NULL", call. = FALSE)
  }

  if (inherits(object, "character")) {
    object <- get(object)
  }

  Extension <- stringr::str_to_lower(tools::file_ext(out_path))

  if (!Extension %in% c("qs2", "rdata", "feather", "rds")) {
    stop(
      "Extension of `out_path` must be either 'qs2', ",
      "'rdata', 'feather', or 'rds' (case-insensitive).", call. = FALSE)
  }

  # Create directory if not available
  fs::dir_create(dirname(out_path))

  switch(
    Extension,
    qs2 = {
      qs2::qs_save(object = object, file = out_path, nthreads = n_threads, ...)
    },
    rdata = {
      if (is.null(object_name)) {
        stop(
          "`object_name` cannot be `NULL` for saving RData files",
          call. = FALSE)
      }
      object_name <- eval(object_name)
      assign(object_name, object)
      save(list = object_name, file = out_path, ...)
    },
    feather = {
      arrow::write_feather(
        x = object, sink = out_path, compression = feather_compression, ...)
    },
    rds = {
      saveRDS(object = object, file = out_path, ...)
    },
    stop("Invalid file extension", call. = FALSE)
  )

  return(invisible())
}
