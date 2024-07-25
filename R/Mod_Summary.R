## |------------------------------------------------------------------------| #
# Mod_Summary ----
## |------------------------------------------------------------------------| #

#' Summary of Hmsc model parameters
#'
#' This function provides a comprehensive summary of Hmsc model parameters, including Alpha, Beta, Rho, and Omega. It processes the model's output, performs statistical summaries, and optionally returns the summarized data.
#'
#' @param Path_Coda String. Path to the `.RData` file containing a coda object.
#' @param EnvFile String. Path to the environment variables file. This file is used to read necessary environmental variables for processing. The default value is `.env`.
#' @param ReturnData Logical. Indicates whether the summarized data should be returned as an R object. If `TRUE`, the function returns a list containing summaries of Alpha, Beta, Rho, and Omega parameters. The default value is `FALSE`, which means the function will not return any data but will save the summaries to a specified directory.
#' @author Ahmed El-Gabbas
#' @return If `ReturnData` is `FALSE` (default), the function does not return anything and saves the summaries to a directory. If `ReturnData` is `TRUE`, it also returns the data as R object.
#' @export

Mod_Summary <- function(Path_Coda = NULL, EnvFile = ".env", ReturnData = FALSE) {

  if (is.null(Path_Coda)) {
    stop("Path_Coda cannot be empty")
  }

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  VarSp <- Variable <- Species <- Sp1_abb <- Sp2_abb <- IAS_ID <- Val <-
    taxon_name <- Species_name <- Sp1_Species_name <- Sp2_Species_name <-
    Q2_5 <- Q97_5 <- Sp1_Sp_abb <- Sp2_Sp_abb <- NULL

  # DataPrep helper function -------
  DataPrep <- function(DT) {
    DT %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "VarSp") %>%
      tibble::as_tibble() %>%
      dplyr::rename(
        dplyr::any_of(c(
          Naive_SE = "Naive SE", TimeSeries_SE = "Time-series SE",
          Q2_5 = "2.5%", Q25 = "25%", Q50 = "50%", Q75 = "75%",
          Q97_5 = "97.5%")))
  }

  # Loading coda object ------
  IASDT.R::CatTime("Loading coda object")
  Coda <- IASDT.R::LoadAs(Path_Coda)

  # Beta ------
  IASDT.R::CatTime("Beta")
  Beta_Summary <- summary(Coda$Beta)
  Beta_Summary <- DataPrep(Beta_Summary$statistics) %>%
    dplyr::full_join(DataPrep(Beta_Summary$quantiles), by = "VarSp") %>%
    dplyr::mutate(
      VarSp = purrr::map(
        .x = VarSp,
        .f = ~{
          stringr::str_split(string = .x, pattern = ",", simplify = TRUE) %>%
            as.vector() %>%
            stringr::str_trim() %>%
            stats::setNames(c("Variable", "Species"))
        })) %>%
    tidyr::unnest_wider(col = VarSp) %>%
    dplyr::mutate(
      Variable = stringr::str_remove_all(Variable, "B\\[| \\(.+\\)|\\(|\\)"),
      Species = stringr::str_remove_all(Species, "^Species_| \\(S.+\\)\\]|\\]"),
      CI_Overlap_0 = purrr::map2_lgl(
        .x = Q2_5, .y = Q97_5, ~dplyr::between(x = 0, left = .x, right = .y)))

  invisible(gc())

  # Alpha -----

  IASDT.R::CatTime("Alpha")
  Alpha_Summary <- summary(Coda$Alpha[[1]])
  Alpha_Summary <- DataPrep(Alpha_Summary$statistics) %>%
    dplyr::full_join(DataPrep(Alpha_Summary$quantiles), by = "VarSp") %>%
    dplyr::mutate(
      VarSp = purrr::map(
        .x = VarSp,
        .f = ~{
          stringr::str_split(string = .x, pattern = "\\[", simplify = TRUE) %>%
            stringr::str_remove(pattern = "\\]") %>%
            as.vector() %>%
            stringr::str_trim() %>%
            stats::setNames(c("Alpha", "Factor"))
        })) %>%
    tidyr::unnest_wider(col = VarSp) %>%
    dplyr::mutate(
      CI_Overlap_0 = purrr::map2_lgl(
        .x = Q2_5, .y = Q97_5, ~dplyr::between(x = 0, left = .x, right = .y)))

  invisible(gc())

  # Rho ----
  IASDT.R::CatTime("Rho")
  Rho_Summary <- summary(Coda$Rho)
  Rho_Summary <- dplyr::bind_rows(
    as.data.frame(as.matrix(Rho_Summary$statistics)),
    as.data.frame(as.matrix(Rho_Summary$quantiles))) %>%
    stats::setNames("Val") %>%
    tibble::rownames_to_column(var = "Var") %>%
    tibble::as_tibble() %>%
    tidyr::pivot_wider(names_from = "Var", values_from = Val) %>%
    dplyr::rename(
      dplyr::any_of(c(
        Naive_SE = "Naive SE", TimeSeries_SE = "Time-series SE",
        Q2_5 = "2.5%", Q25 = "25%", Q50 = "50%", Q75 = "75%", Q97_5 = "97.5%")
      )) %>%
    dplyr::mutate(Rho = "Taxonomy", .before = 1) %>%
    dplyr::mutate(
      CI_Overlap_0 = purrr::map2_lgl(
        .x = Q2_5, .y = Q97_5, ~dplyr::between(x = 0, left = .x, right = .y)))

  # Omega ------

  IASDT.R::CatTime("Omega")

  # Prepare Species list
  if (file.exists(EnvFile)) {
    readRenviron(EnvFile)
  } else {
    MSG <- paste0("Path for environment variables: ", EnvFile, " was not found")
    stop(MSG)
  }
  ListSp <- Sys.getenv("DP_R_Mod_Path_TaxaList") %>%
    file.path("Species_List_ID.txt") %>%
    utils::read.delim(sep = "\t") %>%
    tibble::tibble() %>%
    dplyr::mutate(
      IAS_ID = stringr::str_pad(IAS_ID, pad = "0", width = 4),
      IAS_ID = paste0("Sp_", IAS_ID)) %>%
    dplyr::select(
      Sp_abb = IAS_ID, taxon_name, Species_name, tidyselect::everything())

  # Prepare summary for the omega parameter
  Omega_Summary <- summary(Coda$Omega[[1]])
  Omega_Summary <- DataPrep(Omega_Summary$statistics) %>%
    dplyr::full_join(DataPrep(Omega_Summary$quantiles), by = "VarSp") %>%
    dplyr::mutate(
      VarSp = purrr::map(
        .x = VarSp,
        .f = ~{
          stringr::str_remove_all(.x, pattern = "\\]|Omega1|\\[") %>%
            stringr::str_split(pattern = ",", simplify = TRUE) %>%
            as.data.frame() %>%
            tibble::tibble() %>%
            stats::setNames(c("Sp1_abb", "Sp2_abb")) %>%
            dplyr::mutate_all(stringr::str_trim)
        })) %>%
    tidyr::unnest_wider(col = VarSp) %>%
    dplyr::left_join(
      dplyr::rename_all(ListSp, ~paste0("Sp1_", .x)),
      by = dplyr::join_by(Sp1_abb == Sp1_Sp_abb)) %>%
    dplyr::left_join(
      dplyr::rename_all(ListSp, ~paste0("Sp2_", .x)),
      by = dplyr::join_by(Sp2_abb == Sp2_Sp_abb)) %>%
    dplyr::select(
      Sp1_abb, Sp2_abb, tidyselect::everything(),
      tidyselect::starts_with("Sp1"), tidyselect::starts_with("Sp2"),
      -tidyselect::ends_with("_name2"), -tidyselect::ends_with("_File")) %>%
    dplyr::relocate(Sp1 = Sp1_Species_name, .after = "Sp1_abb") %>%
    dplyr::relocate(Sp2 = Sp2_Species_name, .after = "Sp2_abb") %>%
    dplyr::mutate(
      CI_Overlap_0 = purrr::map2_lgl(
        .x = Q2_5, .y = Q97_5, ~dplyr::between(x = 0, left = .x, right = .y)))

  # Saving ------
  IASDT.R::CatTime("Saving")
  Path_Out <- Path_Coda %>%
    dirname() %>%
    dirname() %>%
    file.path("Model_Postprocessing")
  fs::dir_create(Path_Out)

  IASDT.R::SaveAs(
    InObj = list(
      Alpha = Alpha_Summary, Beta = Beta_Summary,
      Rho = Rho_Summary, Omega = Omega_Summary),
    OutObj = "Parameters_Summary",
    OutPath = file.path(Path_Out, "Parameters_Summary.RData"))

  return(invisible(NULL))
  if (ReturnData) {
    return(InObj = list(
      Alpha = Alpha_Summary, Beta = Beta_Summary,
      Rho = Rho_Summary, Omega = Omega_Summary))
  } else {
    return(invisible(NULL))
  }
}
