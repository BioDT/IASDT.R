# |---------------------------------------------------| #
# ScriptLocation ----
# |---------------------------------------------------| #
#
#' The location of current script
#'
#' The location of current script
#'
#' @name ScriptLocation
#' @references [Click here](https://stackoverflow.com/questions/47044068/)
#' @importFrom rlang .data
#' @export
#' @examples
#'  \dontrun{
#' ScriptLocation()
#' }

ScriptLocation <-  function() {
  this_file <- commandArgs() %>%
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col = .data$value, into = c("key", "value"), sep = "=", fill = "right") %>%
    dplyr::filter(.data$key == "--file") %>%
    dplyr::pull(.data$value)
  if (length(this_file) == 0) {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(this_file)
}
