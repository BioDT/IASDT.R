## -----------------------------------------------------------------------------
library(kableExtra)
library(knitr)
library(tibble)
knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>", eval = FALSE, warning = FALSE,
  message = FALSE, dev = "ragg_png", dpi = 300, tidy = "styler",
  out.width = "100%", fig.show = "hold", echo = FALSE)

## -----------------------------------------------------------------------------
c("<style>", readLines("style.css"), "</style>") %>% 
  cat(sep = "\n")

## -----------------------------------------------------------------------------
DT <- tibble::tribble(
  ~`Climate model`, ~Institution,
  "mpi-esm1-2-hr", "Max Planck Institute for Meteorology, Germany",
  "ipsl-cm6a-lr", "Institut Pierre Simon Laplace, France",
  "ukesm1-0-ll", "Met Office Hadley Centre, UK",
  "gfdl-esm4", "National Oceanic and Atmospheric Administration, USA",
  "mri-esm2-0", "Meteorological Research Institute, Japan") %>% 
  dplyr::mutate(`Climate model` = paste0(`Climate model`, "<tab0>"))

DT %>% 
  knitr::kable(format = "html", escape = FALSE, align = "l")  %>%
  kableExtra::kable_styling(
    c("basic", "condensed", "hover"), font_size = 16, full_width = TRUE) %>%
  kableExtra::add_indent(
    positions = seq_len(nrow(DT)), level_of_indent = 2) %>% 
  kableExtra::column_spec(
    column = 1, extra_css = "white-space: nowrap;", width_min = "3in")

## -----------------------------------------------------------------------------
DT <- tibble::tribble(
  ~"Shared Socioeconomic Pathway", ~Description,
  "ssp126", "SSP1-RCP2.6 climate as simulated by the GCMs",
  "ssp370", "SSP3-RCP7 climate as simulated by the GCMs",
  "ssp585", "SSP5-RCP8.5 climate as simulated by the GCMs") %>% 
  dplyr::mutate(
    `Shared Socioeconomic Pathway` = paste0(
      `Shared Socioeconomic Pathway`, "<tab0>"))

DT %>% 
  knitr::kable(format = "html", escape = FALSE, align = "l")  %>%
  kableExtra::kable_styling(
    c("basic", "condensed", "hover"), font_size = 16, full_width = TRUE) %>%
  kableExtra::add_indent(
    positions = seq_len(nrow(DT)), level_of_indent = 2) %>% 
  kableExtra::column_spec(
    column = 1, extra_css = "white-space: nowrap;", width_min = "3in")

