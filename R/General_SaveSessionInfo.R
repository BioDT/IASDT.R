## |------------------------------------------------------------------------| #
# SaveSessionInfo ----
## |------------------------------------------------------------------------| #
#
#' Save Session Information to a File
#'
#' This function saves the current R session information, including installed packages, session details, and optionally, information about specific objects in the session, to a text file.
#'
#' @param Path  string specifying the directory path where the output file should be saved. The default is the current working directory ([base::getwd]).
#' @param SessionObj An optional list of objects to include in the session information output. This is typically the result of a session management function like [IASDT.R::SaveSession]. If provided, details of these objects (excluding functions and pre-selected objects, with sizes in megabytes) are appended to the session information file.
#' @param Prefix A string to be used as a prefix for the output file name. The default prefix is `S`.
#' @author Ahmed El-Gabbas
#' @return The primary effect of this function is the side effect of writing  session information to a file.
#' @export
#' @name SaveSessionInfo

SaveSessionInfo <- function(Path = getwd(), SessionObj = NULL, Prefix = "S") {

  FileName <- paste0(
    Path, "/", Prefix, "_", format(lubridate::now(tzone = "CET"), "%Y%m%d_%H%M"),
    ".txt")

  utils::capture.output(
    IASDT.R::InfoChunk("Session Info"), file = FileName, append = TRUE)
  utils::capture.output(
    sessioninfo::session_info(), file = FileName, append = TRUE)

  if (magrittr::not(is.null(SessionObj))) {
    utils::capture.output(
      IASDT.R::InfoChunk("Objects in the current session (except functions and\npre-selected objects; Size in megabytes)"),
      file = FileName, append = TRUE)

    sink(FileName, append = TRUE)
    print.data.frame(tibble::tibble(SessionObj), row.names = FALSE)
    sink()
  }

  return(invisible(NULL))
}
