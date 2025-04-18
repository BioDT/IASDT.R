## |------------------------------------------------------------------------| #
# check_data ----
## |------------------------------------------------------------------------| #

#' Check the integrity of `RData` / `qs2` / `rds` / `feather` files
#'
#' These functions validate a given file by checking its extension and
#' attempting to load its contents. A file is considered valid if it loads
#' successfully and contains a non-null object, returning `TRUE`. Otherwise, it
#' returns `FALSE`.
#' @param file Character. The path to the file to be checked. Cannot be empty.
#' @param warning Logical. If `TRUE` (default), warnings will be printed if the
#'   file does not exist.
#' @param n_threads Integer. The number of threads to use when reading qs2
#'   files. Default is 5.
#' @return Logical. `TRUE` if the file is valid, otherwise `FALSE`.
#' @name check_data
#' @rdname check_data
#' @order 1
#' @author Ahmed El-Gabbas
#' @export
#' @details The `check_data()` function determines the file type based on its
#'   extension. If the extension is unrecognized, it returns `FALSE`. Supported
#'   file types:
#' - **RData**: Checked with `check_RData()`, read using [load_as]
#' - **qs2**: Checked with `check_qs()`, read using [qs2::qs_read]
#' - **rds**: Checked with `check_rds()`, read using [readRDS]
#' - **feather**: Checked with `check_feather()`, read using
#'   [arrow::read_feather]

check_data <- function(file, warning = TRUE, n_threads = 5) {

  if (!file.exists(file) || is.null(file) || !nzchar(file)) {
    if (warning) {
      warning("The provided file does not exist: `", file, "`", call. = FALSE)
    }
    return(FALSE)
  }

  Extension <- stringr::str_to_lower(tools::file_ext(file))

  OutFile <- switch(
    Extension,
    qs2 = IASDT.R::check_qs(file, n_threads = n_threads, warning = warning),
    rdata = IASDT.R::check_RData(file, warning = warning),
    rds = IASDT.R::check_rds(file, warning = warning),
    feather = IASDT.R::check_feather(file, warning = warning),
    FALSE)

  return(OutFile)
}


## |------------------------------------------------------------------------| #
# check_RData ----
## |------------------------------------------------------------------------| #

#' @export
#' @name check_data
#' @rdname check_data
#' @order 2
#' @author Ahmed El-Gabbas

check_RData <- function(file, warning = TRUE) {

  if (!file.exists(file) || is.null(file) || !nzchar(file)) {
    if (warning) {
      warning("The provided file does not exist: `", file, "`", call. = FALSE)
    }
    return(FALSE)
  }

  Extension <- stringr::str_to_lower(tools::file_ext(file))

  if (Extension == "rdata") {

    Obj <- try(IASDT.R::load_as(file), silent = TRUE)

    if (inherits(Obj, "try-error")) {
      return(FALSE)
    }

    if (exists("Obj") && !is.null(Obj)) {
      return(TRUE)
    } else {
      return(FALSE)
    }

  } else {
    if (warning) {
      warning("The provided file is not an RData file", call. = FALSE)
    }
    return(FALSE)
  }
}

## |------------------------------------------------------------------------| #
# check_qs ----
## |------------------------------------------------------------------------| #

#' @export
#' @name check_data
#' @rdname check_data
#' @order 3
#' @author Ahmed El-Gabbas

check_qs <- function(file, warning = TRUE, n_threads = 5) {

  if (!file.exists(file) || is.null(file) || !nzchar(file)) {
    if (warning) {
      warning("The provided file does not exist: `", file, "`", call. = FALSE)
    }
    return(FALSE)
  }

  Extension <- stringr::str_to_lower(tools::file_ext(file))

  if (Extension == "qs2") {

    Obj <- try(qs2::qs_read(file = file, nthreads = n_threads), silent = TRUE)

    if (inherits(Obj, "try-error")) {
      return(FALSE)
    }

    if (exists("Obj") && !is.null(Obj)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    if (warning) {
      warning(
        "Unsupported file type. Please provide a `qs2` file.",
        call. = FALSE)
    }
    return(FALSE)

  }
}

## |------------------------------------------------------------------------| #
# check_rds ----
## |------------------------------------------------------------------------| #

#' @export
#' @name check_data
#' @rdname check_data
#' @order 4
#' @author Ahmed El-Gabbas

check_rds <- function(file, warning = TRUE) {

  if (!file.exists(file) || is.null(file) || !nzchar(file)) {
    if (warning) {
      warning("The provided file does not exist: `", file, "`", call. = FALSE)
    }
    return(FALSE)
  }

  Extension <- stringr::str_to_lower(tools::file_ext(file))

  if (Extension == "rds") {

    Obj <- try(readRDS(file), silent = TRUE)

    if (inherits(Obj, "try-error")) {
      return(FALSE)
    }

    if (exists("Obj") && !is.null(Obj)) {
      return(TRUE)
    } else {
      return(FALSE)
    }

  } else {
    if (warning) {
      warning("The provided file is not an rds file", call. = FALSE)
    }
    return(FALSE)
  }
}

## |------------------------------------------------------------------------| #
# check_feather ----
## |------------------------------------------------------------------------| #

#' @export
#' @name check_data
#' @rdname check_data
#' @order 5
#' @author Ahmed El-Gabbas

check_feather <- function(file, warning = TRUE) {

  if (!file.exists(file) || is.null(file) || !nzchar(file)) {
    if (warning) {
      warning("The provided file does not exist: `", file, "`", call. = FALSE)
    }
    return(FALSE)
  }

  Extension <- stringr::str_to_lower(tools::file_ext(file))

  if (Extension == "feather") {

    Obj <- try(arrow::read_feather(file), silent = TRUE)

    if (inherits(Obj, "try-error")) {
      return(FALSE)
    }

    if (exists("Obj") && !is.null(Obj)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    if (warning) {
      warning(
        "Unsupported file type. Please provide a feather file.", call. = FALSE)
    }
    return(FALSE)

  }
}
