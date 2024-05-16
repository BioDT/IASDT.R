# |---------------------------------------------------| #
# LoadMultiple ----
# |---------------------------------------------------| #
#
#' Load multiple RData files together
#'
#' Load multiple RData files together
#' @param Files vector containing the list of `.RData` files to be loaded
#' @param Verbose Verbose? Only valid when `OneObject = FALSE`. Default = `TRUE`
#' @param OneObject Should the objects loaded as one object, each as a list item? If `FALSE`, then all original object names will be loaded to the global environment. Default = `TRUE`.
#' @param ReturnNames Should the name of the objects loaded be returned? Only valid when `OneObject = FALSE`. Default = `TRUE`.
#' @author Ahmed El-Gabbas
#' @return NULL
#' @examples
#' (Files <- system.file("testdata", c("culcita_dat.RData", "gopherdat2.RData"), package = "lme4"))
#'
#' # ---------------------------------------------------
#' # Load multiple *.rd files to one list object
#' # ---------------------------------------------------
#' file.exists(Files)
#' ls()
#'
#' MultiObj <- LoadMultiple(Files = Files, OneObject = TRUE)
#' ls()
#'
#' str(MultiObj, 1)
#'
#' # ---------------------------------------------------
#' # Load multiple *.rd files current environment
#' # ---------------------------------------------------
#' rm("MultiObj")
#' ls()
#'
#' LoadMultiple(Files = Files, OneObject = FALSE)
#'
#' print(c("culcita_dat", "Gdat") %in% ls())
#'
#' str(Gdat, 1)
#'
#' str(culcita_dat, 1)
#'
#' # ---------------------------------------------------
#' # Load multiple *.rd files, one object already exists
#' # ---------------------------------------------------
#' rm("culcita_dat")
#' ls()
#'
#' print(c("culcita_dat", "Gdat") %in% ls())
#'
#' try(LoadMultiple(Files = Files, OneObject = FALSE))
#'
#' @export

LoadMultiple <- function(
    Files = NULL, Verbose = TRUE, OneObject = TRUE, ReturnNames = TRUE) {

  if (!inherits(Files, "character")) {
    cat(paste0(crayon::red("Files should be character object"), "\n"))
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  }

  if (all(file.exists(Files))) {

    if (OneObject) {
      ListNames <- basename(Files) %>%
        stringr::str_replace_all(".RData", "")

      ObjectsLoaded00 <- lapply(
        X = seq_along(Files),
        FUN = function(y) {
          assign(x = ListNames[y], value = LoadAs(Files[y]))
        }) %>%
        stats::setNames(ListNames)

      return(ObjectsLoaded00)

    } else {

      Env1 <- new.env()
      ObjectsLoaded <- lapply(
        X = seq_along(Files),
        FUN = function(y) {
          load(file = Files[y], envir = Env1)
        })

      if (any(ls(envir = Env1) %in% ls(envir = parent.frame()))) {
        cat(paste0(crayon::red("ERROR: Some of the new object names already exists in the current environment. No files were loaded!"), "\n"), sep = "")
        opt <- options(show.error.messages = FALSE)
        on.exit(options(opt))
        stop()
      } else {

        ObjectsLoaded <- lapply(
          X = seq_along(Files),
          FUN = function(y) {
            load(file = Files[y], envir = parent.frame(3))
          }) %>%
          do.call(what = c)

        if (Verbose) {
          cat(
            paste0(crayon::blue(
              "Object:", crayon::red(ObjectsLoaded), "was loaded successfully.")),
            sep = "\n")
        }
        if (ReturnNames) {
          return(ObjectsLoaded)
        } else {
          return(invisible(NULL))
        }
      }

    }
  } else {
    cat(paste0(crayon::red("Some of these files do not exist. No objects were loaded!"), "\n"))
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  }
  # return(invisible(NULL))
}
