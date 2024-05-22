# |---------------------------------------------------| #
# SaveSession ----
# |---------------------------------------------------| #
#
#' Save all objects (except functions) of the global environment as list items
#'
#' Save all objects (except functions) of the global environment as list items
#'
#' @param Path String. Path of where to save the output RData file
#' @param ExcludeObj Objects not to save
#' @param Prefix String. File prefix
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

SaveSession <- function(Path = getwd(), ExcludeObj = NULL, Prefix = "S") {

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Object <- Class <- Size <- Obj <- NULL

  fs::dir_create(Path)
  ExcludeObj <- c(ExcludeObj, "Grid_10_sf_s", "Grid_10_Raster", "Bound_sf_Eur_s", "Bound_sf_Eur")

  AllObjs <- ls(envir = .GlobalEnv) %>%
    tibble::tibble(Object = .) %>%
    dplyr::mutate(
      Class = purrr::map_chr(
        .x = Object,
        .f = ~{
          get(.x, envir = .GlobalEnv) %>%
            class() %>%
            stringr::str_c(collapse = "_")
        }
      )) %>%
    dplyr::filter(
      Class != "function",
      magrittr::not(Object %in% ExcludeObj)) %>%
    dplyr::pull(Object)

  AllObjs <- AllObjs %>%
    purrr::map(
      .f = ~{
        Obj <- get(.x, envir = .GlobalEnv)
        if (class(Obj)[1] == "SpatRaster") {
          suppressWarnings(terra::wrap(Obj))
        } else {
          Obj
        }
      }) %>%
    stats::setNames(AllObjs)

  FF2 <- lubridate::now(tzone = "CET") %>%
    purrr::map_chr(
      .f = ~{
        c(lubridate::year(.x), lubridate::month(.x),
          lubridate::day(.x), "__",
          lubridate::hour(.x), lubridate::minute(.x)) %>%
          sapply(stringr::str_pad, width = 2, pad = "0") %>%
          stringr::str_c(collapse = "") %>%
          stringr::str_replace_all("__", "_") %>%
          stringr::str_c(Prefix, "_", ., collapse = "_")
      })

  IASDT.R::SaveAs(
    InObj = AllObjs, OutObj = FF2,
    OutPath = file.path(Path, paste0(FF2, ".RData")))

  AllObjs %>%
    lapply(pryr::object_size) %>%
    tibble::tibble(Obj = names(.), Size = as.numeric(.)) %>%
    dplyr::mutate(Size = Size / (1024 * 1024), Size = round(Size, 1)) %>%
    dplyr::select(Obj, Size) %>%
    dplyr::arrange(dplyr::desc(Size)) %>%
    return()
}
