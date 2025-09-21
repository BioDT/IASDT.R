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

