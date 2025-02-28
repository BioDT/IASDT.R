## |------------------------------------------------------------------------| #
# SaveSession ----
## |------------------------------------------------------------------------| #

#' Saves all non-function objects from the global environment to an RData file
#'
#' This function saves all objects (except functions and specified exclusions)
#' from the global environment as list items in an `RData` file. It also creates
#' a summary of these objects' sizes in memory.
#'
#' @param Path Character. Directory path where the output `RData` file should be
#'   saved. Defaults to the current working directory [base::getwd()].
#' @param ExcludeObj Character vector. Object names (as strings) to exclude from
#'   saving.
#' @param Prefix Character. Prefix the saved file name with. Defaults to `S`.
#' @author Ahmed El-Gabbas
#' @return A tibble containing the names and sizes (in MB, rounded to 1 decimal
#'   place) of the saved objects.
#' @export
#' @name SaveSession

SaveSession <- function(Path = getwd(), ExcludeObj = NULL, Prefix = "S") {

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Object <- Class <- Size <- Obj <- NULL

  fs::dir_create(Path)

  ExcludeObj <- c(ExcludeObj, "Grid_10_sf_s", "Grid_10_Raster",
                  "Bound_sf_Eur_s", "Bound_sf_Eur")

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
    dplyr::filter(Class != "function", !(Object %in% ExcludeObj)) %>%
    dplyr::pull(Object)

  AllObjs <- purrr::map(
      .x = AllObjs,
      .f = ~{
        Obj <- get(.x, envir = .GlobalEnv)
        if (class(Obj)[1] == "SpatRaster") {
          suppressWarnings(terra::wrap(Obj))
        } else {
          Obj
        }
      }) %>%
    stats::setNames(AllObjs)

  FF2 <- purrr::map_chr(
      .x = lubridate::now(tzone = "CET"),
      .f = ~{
        c(lubridate::year(.x), lubridate::month(.x),
          lubridate::day(.x), "__",
          lubridate::hour(.x), lubridate::minute(.x)) %>%
          purrr::map_chr(stringr::str_pad, width = 2, pad = "0") %>%
          stringr::str_c(collapse = "") %>%
          stringr::str_replace_all("__", "_") %>%
          stringr::str_c(Prefix, "_", ., collapse = "_")
      })

  IASDT.R::SaveAs(
    InObj = AllObjs, OutObj = FF2,
    OutPath = IASDT.R::Path(Path, paste0(FF2, ".RData")))

  Out <- AllObjs %>%
    lapply(lobstr::obj_size) %>%
    tibble::tibble(Obj = names(.), Size = as.numeric(.)) %>%
    dplyr::mutate(Size = Size / (1024 * 1024), Size = round(Size, 1)) %>%
    dplyr::select(Obj, Size) %>%
    dplyr::arrange(dplyr::desc(Size))

  return(Out)
}
