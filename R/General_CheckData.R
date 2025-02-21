## |------------------------------------------------------------------------| #
# CheckData ----
## |------------------------------------------------------------------------| #

#' Check the integrity of `RData` / `qs2` / `rds` / `feather` files
#'
#' These functions validate a given file by checking its extension and
#' attempting to load its contents. A file is considered valid if it loads
#' successfully and contains a non-null object, returning `TRUE`. Otherwise, it
#' returns `FALSE`.
#' @param File Character. The path to the file to be checked. Cannot be empty.
#' @param warning logical. If `TRUE` (default), warnings will be printed if the
#'   file does not exist.
#' @param nthreads Integer. The number of threads to use when reading qs2 files.
#'   Default is 5.
#' @return Logical. `TRUE` if the file is valid, otherwise `FALSE`.
#' @name Check_data
#' @rdname Check_data
#' @order 1
#' @author Ahmed El-Gabbas
#' @export
#' @details The `CheckData()` function determines the file type based on its
#'   extension. If the extension is unrecognized, it returns `FALSE`. Supported
#'   file types:
#' - **RData**: Checked with `CheckRData()`, read using [LoadAs]
#' - **qs2**: Checked with `CheckQs()`, read using [qs2::qs_read]
#' - **rds**: Checked with `CheckRDS()`, read using [readRDS]
#' - **feather**: Checked with `CheckFeather()`, read using
#' [arrow::read_feather]

CheckData <- function(File, warning = TRUE, nthreads = 5) {

  if (!file.exists(File) || is.null(File) || !nzchar(File)) {
    if (warning) {
      warning(paste0("The provided file does not exist: `", File, "`"))
    }
    return(FALSE)
  }

  Extension <- stringr::str_to_lower(tools::file_ext(File))

  OutFile <- switch(
    Extension,
    qs2 = IASDT.R::CheckQs(File, nthreads = nthreads, warning = warning),
    rdata = IASDT.R::CheckRData(File, warning = warning),
    rds = IASDT.R::CheckRDS(File, warning = warning),
    feather = IASDT.R::CheckFeather(File, warning = warning),
    FALSE)

  return(OutFile)
}


## |------------------------------------------------------------------------| #
# CheckRData ----
## |------------------------------------------------------------------------| #

#' @export
#' @name Check_data
#' @rdname Check_data
#' @order 2
#' @author Ahmed El-Gabbas

CheckRData <- function(File, warning = TRUE) {

  if (!file.exists(File) || is.null(File) || !nzchar(File)) {
    if (warning) {
      warning(paste0("The provided file does not exist: `", File, "`"))
    }
    return(FALSE)
  }

  Extension <- stringr::str_to_lower(tools::file_ext(File))

  if (Extension == "rdata") {

    Obj <- try(IASDT.R::LoadAs(File), silent = TRUE)

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
# CheckQs ----
## |------------------------------------------------------------------------| #

#' @export
#' @name Check_data
#' @rdname Check_data
#' @order 3
#' @author Ahmed El-Gabbas

CheckQs <- function(File, warning = TRUE, nthreads = 5) {

  if (!file.exists(File) || is.null(File) || !nzchar(File)) {
    if (warning) {
      warning(paste0("The provided file does not exist: `", File, "`"))
    }
    return(FALSE)
  }

  Extension <- stringr::str_to_lower(tools::file_ext(File))

  if (Extension == "qs2") {

    Obj <- try(qs2::qs_read(file = File, nthreads = nthreads), silent = TRUE)

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
# CheckRDS ----
## |------------------------------------------------------------------------| #

#' @export
#' @name Check_data
#' @rdname Check_data
#' @order 4
#' @author Ahmed El-Gabbas

CheckRDS <- function(File, warning = TRUE) {

  if (!file.exists(File) || is.null(File) || !nzchar(File)) {
    if (warning) {
      warning(paste0("The provided file does not exist: `", File, "`"))
    }
    return(FALSE)
  }

  Extension <- stringr::str_to_lower(tools::file_ext(File))

  if (Extension == "rds") {

    Obj <- try(readRDS(File), silent = TRUE)

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
      warning("The provided file is not an rds file")
    }
    return(FALSE)
  }
}

## |------------------------------------------------------------------------| #
# CheckFeather ----
## |------------------------------------------------------------------------| #

#' @export
#' @name Check_data
#' @rdname Check_data
#' @order 5
#' @author Ahmed El-Gabbas

CheckFeather <- function(File, warning = TRUE) {

  if (!file.exists(File) || is.null(File) || !nzchar(File)) {
    if (warning) {
      warning(paste0("The provided file does not exist: `", File, "`"))
    }
    return(FALSE)
  }

  Extension <- stringr::str_to_lower(tools::file_ext(File))

  if (Extension == "feather") {

    Obj <- try(arrow::read_feather(File), silent = TRUE)

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
