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
  ~"class",         ~"order",            ~`# species`,
  "Liliopsida",     "Acorales",           2,
  "Liliopsida",     "Alismatales",        8,
  "Liliopsida",     "Arecales",           5,
  "Liliopsida",     "Asparagales",       74,
  "Liliopsida",     "Commelinales",       4,
  "Liliopsida",     "Liliales",           6,
  "Liliopsida",     "Poales",           163,
  "Liliopsida",     "Zingiberales",       7,
  "Lycopodiopsida", "Selaginellales",     1,
  "Magnoliopsida",  "Apiales",           15,
  "Magnoliopsida",  "Asterales",        148,
  "Magnoliopsida",  "Boraginales",       10,
  "Magnoliopsida",  "Brassicales",       13,
  "Magnoliopsida",  "Caryophyllales",   141,
  "Magnoliopsida",  "Celastrales",        4,
  "Magnoliopsida",  "Cornales",          10,
  "Magnoliopsida",  "Cucurbitales",       6,
  "Magnoliopsida",  "Dipsacales",        19,
  "Magnoliopsida",  "Ericales",          40,
  "Magnoliopsida",  "Escalloniales",      1,
  "Magnoliopsida",  "Fabales",           52,
  "Magnoliopsida",  "Fagales",           20,
  "Magnoliopsida",  "Garryales",          1,
  "Magnoliopsida",  "Gentianales",        9,
  "Magnoliopsida",  "Geraniales",         9,
  "Magnoliopsida",  "Gunnerales",         1,
  "Magnoliopsida",  "Lamiales",          63,
  "Magnoliopsida",  "Laurales",           1,
  "Magnoliopsida",  "Magnoliales",        1,
  "Magnoliopsida",  "Malpighiales",      41,
  "Magnoliopsida",  "Malvales",          10,
  "Magnoliopsida",  "Myrtales",          74,
  "Magnoliopsida",  "Oxalidales",        19,
  "Magnoliopsida",  "Piperales",          3,
  "Magnoliopsida",  "Proteales",          4,
  "Magnoliopsida",  "Ranunculales",      27,
  "Magnoliopsida",  "Rosales",          125,
  "Magnoliopsida",  "Sapindales",        23,
  "Magnoliopsida",  "Saxifragales",      38,
  "Magnoliopsida",  "Solanales",         51,
  "Magnoliopsida",  "Vitales",            9,
  "Pinopsida",      "Pinales",           46,
  "Polypodiopsida", "Cyatheales",         1,
  "Polypodiopsida", "Equisetales",        1,
  "Polypodiopsida", "Polypodiales",      12,
  "Polypodiopsida", "Schizaeales",        1)

knitr::kable(
  cbind(DT[1:23,], DT[24:46,]), format = "html", 
  align = "l", escape = FALSE) %>% 
  kableExtra::column_spec(c(1, 2, 4, 5), width_min = "1.4in")  %>% 
  kableExtra::column_spec(c(3, 6), width_min = "1.2in")  %>%
  kableExtra::kable_styling(
    c("basic", "condensed", "hover"), font_size = 16,
    position = "center", full_width = TRUE)

