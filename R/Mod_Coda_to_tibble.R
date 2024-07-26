## |------------------------------------------------------------------------| #
# Coda_to_tibble ----
## |------------------------------------------------------------------------| #

#' Convert a Coda object to a tibble with specified parameter transformations
#'
#' This function converts a Coda object (`mcmc.list` or `mcmc`) into a tibble format, facilitating further analysis and visualization. It supports transformation for specific parameter types: rho, alpha, omega, and beta.
#'
#' @param CodaObj An object of class `mcmc.list` or `mcmc`, representing the MCMC output.
#' @param Type A character string specifying the parameter type to transform and extract. Must be one of `rho`, `alpha`, `omega`, or `beta`.
#' @param EnvFile A character string specifying the path to the environment file that contains necessary variables for beta parameter transformation. Defaults to `.env`.
#' @param NOmega An integer specifying the number of species to be sampled for the Omega parameter transformation. Defaults to 1000.
#' @name Coda_to_tibble
#' @author Ahmed El-Gabbas
#' @return A tibble containing the transformed parameters based on the specified `Type`. The structure of the returned tibble varies depending on the `Type` parameter.
#' @export
#' @examples
#' library(Hmsc)
#' library(coda)
#' library(dplyr)
#' Coda  <- Hmsc::convertToCodaObject(Hmsc::TD$m)
#' IASDT.R::Coda_to_tibble(CodaObj = Coda$Alpha[[1]], Type = "Alpha")

Coda_to_tibble <- function(
    CodaObj = NULL, Type = NULL, EnvFile = ".env", NOmega = 1000) {

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Chain <- Iter <- Alpha <- AlphaNum <- Factor <- CHAIN <- ITER <-
    Value <- IAS_ID <- Species_name <- Var_Sp <- Variable <- Species <-
    SpComb <- Sp1 <- IAS1 <- Sp2 <- IAS2 <- NULL

  # # |||||||||||||||||||||||||||||||||||||||
  # Check missing arguments ----
  # # |||||||||||||||||||||||||||||||||||||||

  if (is.null(CodaObj) || is.null(Type)) {
    stop("None of `CodaObj` or `Type` cannot be empty")
  }

  # # |||||||||||||||||||||||||||||||||||||||
  # Check class of input coda object -----
  # # |||||||||||||||||||||||||||||||||||||||

  Type <- stringr::str_to_lower(Type)

  if (!Type %in% c("rho", "alpha", "omega", "beta")) {
    stop("Type has to be one of rho, alpha, omega, or beta")
  }

  if (!(inherits(CodaObj, "mcmc.list") && !inherits(CodaObj, "mcmc"))) {
    stop("Input Coda object has to be of class mcmc.list or mcmc")
  }


  # # |||||||||||||||||||||||||||||||||||||||
  # Convert to tibble ----
  # # |||||||||||||||||||||||||||||||||||||||

  Coda <- as.matrix(CodaObj, iter = TRUE, chain = TRUE) %>%
    tibble::as_tibble() %>%
    dplyr::rename(Chain = CHAIN, Iter = ITER) %>%
    dplyr::mutate(Chain = factor(Chain), Iter = as.integer(Iter)) %>%
    dplyr::arrange(Chain, Iter)

  # # |||||||||||||||||||||||||||||||||||||||
  # Rho parameter ----
  # # |||||||||||||||||||||||||||||||||||||||

  if (Type == "rho") {
    Coda <- dplyr::rename(Coda, dplyr::any_of(c(Value = "var1")))
  }

  # # |||||||||||||||||||||||||||||||||||||||
  # Alpha parameter ----
  # # |||||||||||||||||||||||||||||||||||||||

  if (Type == "alpha") {

    Coda <- tidyr::pivot_longer(
      data = Coda,
      cols = -c(Chain, Iter), names_to = "Alpha", values_to = "Value")

    Coda <- Coda %>%
      dplyr::distinct(Alpha) %>%
      dplyr::mutate(
        Alpha2 = purrr::map(
          .x = Alpha,
          .f = ~{
            stringr::str_split(.x, "\\[", simplify = TRUE) %>%
              as.character() %>%
              stringr::str_remove("\\]")
          })) %>%
      tidyr::unnest_wider("Alpha2", names_sep = "_") %>%
      purrr::set_names(c("Alpha", "AlphaNum", "Factor")) %>%
      dplyr::right_join(Coda, by = "Alpha") %>%
      dplyr::mutate(
        Alpha = factor(Alpha), AlphaNum = factor(AlphaNum),
        Factor = factor(Factor)) %>%
      dplyr::select(
        Alpha, AlphaNum, Factor, Chain, Iter, Value, dplyr::everything()) %>%
      dplyr::arrange(Alpha, Chain, Iter)
  }

  # # |||||||||||||||||||||||||||||||||||||||
  # Beta parameters ------
  # # |||||||||||||||||||||||||||||||||||||||

  if (Type == "beta") {

    if (file.exists(EnvFile)) {

      readRenviron(EnvFile)

      SpeciesNames <- Sys.getenv("DP_R_Mod_Path_TaxaList") %>%
        file.path("Species_List_ID.txt") %>%
        readr::read_tsv(show_col_types = FALSE) %>%
        dplyr::select(IAS_ID, Species = Species_name) %>%
        dplyr::mutate(
          IAS_ID = stringr::str_pad(IAS_ID, pad = "0", width = 4),
          IAS_ID = paste0("Sp_", IAS_ID))

    } else {
      stop("Path for environment variables: ", EnvFile, " was not found")
    }

    Coda <- tidyr::pivot_longer(
      Coda, -c(Chain, Iter), names_to = "Var_Sp", values_to = "Value")

    Coda <- Coda %>%
      dplyr::distinct(Var_Sp) %>%
      dplyr::mutate(
        Betas = purrr::map(
          .x = Var_Sp,
          .f = ~{
            .x %>%
              stringr::str_remove_all("B\\[|\\(|\\)|\\]") %>%
              stringr::str_split(", ", simplify = TRUE) %>%
              as.vector() %>%
              as.list() %>%
              purrr::set_names(c("Variable", "IAS_ID"))
          })) %>%
      tidyr::unnest_wider("Betas") %>%
      dplyr::right_join(Coda, by = "Var_Sp") %>%
      dplyr::left_join(SpeciesNames, by = "IAS_ID") %>%
      dplyr::select(
        Var_Sp, Variable, Species, IAS_ID, Chain, Iter, Value,
        dplyr::everything()) %>%
      tidyr::nest(DT = -c(Variable, IAS_ID, Species, Var_Sp)) %>%
      dplyr::arrange(Variable, IAS_ID)
  }

  # # |||||||||||||||||||||||||||||||||||||||
  # Omega parameters ------
  # # |||||||||||||||||||||||||||||||||||||||

  if (Type == "omega") {

    SampleOmega <- setdiff(names(Coda), c("Chain", "Iter")) %>%
      sample(size = min((ncol(Coda) - 2), NOmega))

    Coda <- Coda %>%
      dplyr::select(Chain, Iter, dplyr::all_of(SampleOmega)) %>%
      tidyr::pivot_longer(
        -c(Chain, Iter), names_to = "SpComb", values_to = "Value")

    Coda <- tibble::tibble(SpComb = SampleOmega) %>%
      dplyr::mutate(
        SP = purrr::map(
          .x = SpComb,
          .f = ~{
            IAS_N <- .x %>%
              stringr::str_remove_all("Omega1\\[|\\]") %>%
              stringr::str_split(", ", simplify = TRUE) %>%
              as.character() %>%
              purrr::set_names(c("Sp1", "Sp2"))

            IAS1 <- IASDT.R::GetSpeciesName(
              EnvFile = EnvFile, SpID = IAS_N[1]) %>%
              dplyr::pull(Species_name)
            IAS2 <- IASDT.R::GetSpeciesName(
              EnvFile = EnvFile, SpID = IAS_N[2]) %>%
              dplyr::pull(Species_name)
            return(c(IAS_N, IAS1 = IAS1, IAS2 = IAS2))
          })) %>%
      tidyr::unnest_wider("SP") %>%
      dplyr::right_join(Coda, by = "SpComb") %>%
      dplyr::select(
        SpComb, Sp1, IAS1, Sp2, IAS2, Chain, Iter, Value,
        dplyr::everything()) %>%
      tidyr::nest(DT = -c(SpComb, Sp1, IAS1, Sp2, IAS2)) %>%
      dplyr::arrange(SpComb)

  }
  return(Coda)
}
