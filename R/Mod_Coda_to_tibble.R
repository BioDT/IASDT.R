## |------------------------------------------------------------------------| #
# Coda_to_tibble ----
## |------------------------------------------------------------------------| #

#' Convert `mcmc.list` to tibble
#'
#' Convert `mcmc.list` to tibble
#'
#' @param CodaObj `mcmc.list` object
#' @param Type String. Parameter type. It has to be one of `rho`, `alpha`, `omega`, or `beta`.
#' @param EnvFile String. Path to read the environment variables. Default value: `.env`
#' @param NOmega Integer. Number of species to be sampled for the Omega parameter
#' @name Coda_to_tibble
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

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

  if (any(missing(CodaObj), missing(Type))) {
    MSG <- "None of `CodaObj` or `Type` can be empty"
    stop(MSG)
  }

  # # |||||||||||||||||||||||||||||||||||||||
  # Check class of input coda object -----
  # # |||||||||||||||||||||||||||||||||||||||

  Type <- stringr::str_to_lower(Type)

  if (sum(stringr::str_detect(Type, c("rho", "alpha", "omega", "beta"))) == 0) {
    MSG <- "Type has to be one of rho, alpha, omega, or beta"
    stop(MSG)
  }

  CodaClass <- c(
    inherits(CodaObj, "mcmc.list"), inherits(CodaObj, "mcmc")) %>%
    sum() %>%
    magrittr::equals(0)

  if (CodaClass) {
    stop("Input Coda object has to be of class mcmc.list or mcmc")
  }

  # # |||||||||||||||||||||||||||||||||||||||
  # Convert to tibble ----
  # # |||||||||||||||||||||||||||||||||||||||

  Coda <- CodaObj %>%
    as.matrix(iter = TRUE, chain = TRUE) %>%
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
    Coda <- Coda %>%
      tidyr::pivot_longer(
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
      MSG <- paste0(
        "Path for environment variables: ", EnvFile, " was not found")
      stop(MSG)
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

    SampleOmega <- names(Coda) %>%
      setdiff(c("Chain", "Iter")) %>%
      sample(size = NOmega)

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

            IAS1 <- dplyr::pull(GetSpeciesName(SpID = IAS_N[1]), Species_name)
            IAS2 <- dplyr::pull(GetSpeciesName(SpID = IAS_N[2]), Species_name)
            return(c(IAS_N, IAS1 = IAS1, IAS2 = IAS2))
          })) %>%
      tidyr::unnest_wider("SP") %>%
      dplyr::right_join(Coda, by = "SpComb") %>%
      dplyr::select(
        SpComb, Sp1, IAS1, Sp2, IAS2, Chain, Iter, Value, dplyr::everything()) %>%
      tidyr::nest(DT = -c(SpComb, Sp1, IAS1, Sp2, IAS2)) %>%
      dplyr::arrange(SpComb)
  }
  return(Coda)
}
