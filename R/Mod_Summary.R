## |------------------------------------------------------------------------| #
# BetaSummary ----
## |------------------------------------------------------------------------| #

#' BetaSummary
#'
#' BetaSummary
#'
#' @param Beta Beta
#' @name BetaSummary
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

BetaSummary <- function(Beta) {

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  "Naive SE" <- "Time-series SE" <- `2.5%` <- `25%` <- `50%` <-
    `75%` <- `97.5%` <- rowname <- Var <- Sp <- `Naive SE` <-
    `Time-series SE` <- NULL

  Stats <- summary(Beta)$statistics %>%
    round(3) %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    tibble::as_tibble()

  Quant <- summary(Beta)$quantiles %>%
    round(3) %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    tibble::as_tibble()

  Summary <- dplyr::full_join(Stats, Quant, by = "rowname") %>%
    dplyr::rename(
      Naive_SE = `Naive SE`, TimeSeries_SE = `Time-series SE`,
      Q2_5 = `2.5%`, Q25 = `25%`, Q50 = `50%`,
      Q75 = `75%`, Q97_5 = `97.5%`) %>%
    dplyr::mutate(
      rowname = purrr::map(
        .x = rowname,
        .f = ~{
          stringr::str_split(
            string = .x, pattern = ",", simplify = TRUE) %>%
            as.vector() %>%
            stringr::str_trim() %>%
            setNames(c("Var", "Sp"))
        })) %>%
    tidyr::unnest_wider(col = rowname) %>%
    dplyr::mutate(
      Var = stringr::str_remove_all(Var, "B\\[| \\(.+\\)|\\(|\\)"),
      Sp = stringr::str_remove_all(Sp, "^Sp_| \\(S.+\\)\\]|\\]"))
  return(Summary)
}
