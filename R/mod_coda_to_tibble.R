## |------------------------------------------------------------------------| #
# coda_to_tibble ----
## |------------------------------------------------------------------------| #

#' Convert a Coda object to a tibble with specified parameter transformations
#'
#' This function converts a Coda object (`mcmc.list` or `mcmc`) into a tibble
#' format, facilitating further analysis and visualization. It supports
#' transformation for specific parameter types: `rho`, `alpha`, `omega`, and
#' `beta`.
#' @param coda_object An object of class `mcmc.list` or `mcmc`, representing the
#'   MCMC output.
#' @param posterior_type Character. The parameter type to transform and extract.
#'   Must be one of `rho`, `alpha`, `omega`, or `beta`.
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param n_omega Integer. The number of species to be sampled for the `Omega`
#'   parameter transformation. Defaults to 100.
#' @name coda_to_tibble
#' @author Ahmed El-Gabbas
#' @return A tibble containing the transformed parameters based on the specified
#'   `posterior_type`. The structure of the returned tibble varies depending on
#'   the `posterior_type` parameter.
#' @export
#' @examples
#' library(Hmsc)
#' library(coda)
#' library(dplyr)
#' Coda <- Hmsc::convertToCodaObject(Hmsc::TD$m)
#' IASDT.R::coda_to_tibble(
#'    coda_object = Coda$Alpha[[1]], posterior_type = "Alpha")

coda_to_tibble <- function(
    coda_object = NULL, posterior_type = NULL,
    env_file = ".env", n_omega = 100) {

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Chain <- Iter <- Alpha <- AlphaNum <- Factor <- CHAIN <- ITER <-
    Value <- IAS_ID <- Species_name <- Var_Sp <- Variable <- Species <-
    SpComb <- Sp1 <- IAS1 <- Sp2 <- IAS2 <- TaxaInfoFile <- NULL

  # # |||||||||||||||||||||||||||||||||||||||
  # Check missing arguments ----
  # # |||||||||||||||||||||||||||||||||||||||

  if (is.null(coda_object) || is.null(posterior_type)) {
    IASDT.R::stop_ctx(
      "None of `coda_object` or `posterior_type` cannot be empty",
      coda_object = coda_object, posterior_type = posterior_type)
  }

  # # |||||||||||||||||||||||||||||||||||||||
  # Check class of input coda object -----
  # # |||||||||||||||||||||||||||||||||||||||

  posterior_type <- tolower(posterior_type)

  if (!posterior_type %in% c("rho", "alpha", "omega", "beta")) {
    IASDT.R::stop_ctx(
      "posterior_type has to be one of rho, alpha, omega, or beta",
      posterior_type = posterior_type)
  }

  if (!(inherits(coda_object, "mcmc.list") || inherits(coda_object, "mcmc"))) {
    IASDT.R::stop_ctx(
      "Input Coda object has to be of class mcmc.list or mcmc",
      class_coda_obj = class(coda_object))
  }

  # # |||||||||||||||||||||||||||||||||||||||
  # Convert to tibble ----
  # # |||||||||||||||||||||||||||||||||||||||

  if (posterior_type == "omega") {
    CombSample <- sample.int(n = dim(coda_object[[1]])[2], size = n_omega)

    Coda <- purrr::map(coda_object, ~.x[, CombSample]) %>%
      coda::as.mcmc.list() %>%
      as.matrix(iter = TRUE, chain = TRUE) %>%
      tibble::as_tibble() %>%
      dplyr::rename(Chain = CHAIN, Iter = ITER) %>%
      dplyr::mutate(Chain = factor(Chain), Iter = as.integer(Iter)) %>%
      dplyr::arrange(Chain, Iter) %>%
      tidyr::pivot_longer(
        -c(Chain, Iter), names_to = "SpComb", values_to = "Value")
  } else {
    Coda <- as.matrix(coda_object, iter = TRUE, chain = TRUE) %>%
      tibble::as_tibble() %>%
      dplyr::rename(Chain = CHAIN, Iter = ITER) %>%
      dplyr::mutate(Chain = factor(Chain), Iter = as.integer(Iter)) %>%
      dplyr::arrange(Chain, Iter)
  }


  # # |||||||||||||||||||||||||||||||||||||||
  # Rho parameter ----
  # # |||||||||||||||||||||||||||||||||||||||

  if (posterior_type == "rho") {
    Coda <- dplyr::rename(Coda, dplyr::any_of(c(Value = "var1")))
  }

  # # |||||||||||||||||||||||||||||||||||||||
  # Alpha parameter ----
  # # |||||||||||||||||||||||||||||||||||||||

  if (posterior_type == "alpha") {

    Coda <- tidyr::pivot_longer(
      data = Coda,
      cols = -c(Chain, Iter), names_to = "Alpha", values_to = "Value")

    Coda <- dplyr::distinct(Coda, Alpha) %>%
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

  if (posterior_type == "beta") {

    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "TaxaInfoFile", "DP_R_Taxa_info", FALSE, TRUE)
    # Assign environment variables and check file and paths
    IASDT.R::assign_env_vars(
      env_file = env_file, env_variables_data = EnvVars2Read)
    rm(EnvVars2Read, envir = environment())

    SpeciesNames <- readr::read_tsv(
      file = TaxaInfoFile, show_col_types = FALSE, progress = FALSE) %>%
      dplyr::select(IAS_ID, Species = Species_name) %>%
      dplyr::mutate(
        IAS_ID = stringr::str_pad(IAS_ID, pad = "0", width = 4),
        IAS_ID = paste0("Sp_", IAS_ID))

    Coda <- tidyr::pivot_longer(
      Coda, -c(Chain, Iter), names_to = "Var_Sp", values_to = "Value")


    VarSp <- Coda %>%
      dplyr::distinct(Var_Sp) %>%
      dplyr::mutate(
        Species = purrr::map(
          .x = Var_Sp,
          .f = ~{
            .x %>%
              stringr::str_remove_all("B\\[|B\\[|\\(|\\)|\\]|stats::poly") %>%
              stringr::str_replace_all(", degree = 2, raw = TRUE", "_") %>%
              stringr::str_split(",", simplify = TRUE) %>%
              as.data.frame() %>%
              tibble::as_tibble() %>%
              stats::setNames(c("Variable", "IAS_ID")) %>%
              dplyr::mutate_all(stringr::str_trim) %>%
              dplyr::mutate(
                Variable = purrr::map_chr(
                  .x = Variable,
                  .f = function(V) {
                    V %>%
                      stringr::str_replace_all("_1$", "_L") %>%
                      stringr::str_replace_all("_2$", "_Q")
                  }))
          })) %>%
      tidyr::unnest_wider("Species") %>%
      dplyr::left_join(SpeciesNames, by = "IAS_ID")



    Coda <- VarSp %>%
      dplyr::right_join(Coda, by = "Var_Sp") %>%
      dplyr::select(
        Var_Sp, Variable, Species, IAS_ID, Chain, Iter, Value,
        dplyr::everything()) %>%
      tidyr::nest(DT = -c(Variable, IAS_ID, Species, Var_Sp)) %>%
      dplyr::arrange(Variable, IAS_ID)
  }

  # # |||||||||||||||||||||||||||||||||||||||
  # Omega parameters ------
  # # |||||||||||||||||||||||||||||||||||||||

  if (posterior_type == "omega") {

    IAS <- IASDT.R::get_species_name(env_file = env_file) %>%
      dplyr::select(IAS_ID, Species_name)

    Coda <- tibble::tibble(SpComb = unique(Coda$SpComb)) %>%
      dplyr::mutate(
        SP = purrr::map(
          .x = SpComb,
          .f = ~{
            IAS_N <- stringr::str_remove_all(.x, "Omega1\\[|\\]") %>%
              stringr::str_split(",| ", simplify = TRUE) %>%
              as.character() %>%
              stringr::str_subset("^Sp_") %>%
              purrr::set_names(c("Sp1", "Sp2"))
            IAS1 <- dplyr::filter(IAS, IAS_ID == IAS_N[1])$Species_name
            IAS2 <- dplyr::filter(IAS, IAS_ID == IAS_N[2])$Species_name

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
