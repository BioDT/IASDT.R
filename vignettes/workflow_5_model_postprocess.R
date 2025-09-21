## -----------------------------------------------------------------------------
cat('
<script>
function toggleScript(id) {
  var content = document.getElementById(id);
  if (content.style.display === "none") {
    content.style.display = "block";
  } else {
    content.style.display = "none";
  }
}
</script>
')

## -----------------------------------------------------------------------------
library(kableExtra)
library(knitr)
library(tibble)
knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>", eval = FALSE, warning = FALSE,
  message = FALSE, dev = "ragg_png", dpi = 300, tidy = "styler",
  out.width = "100%", fig.show = "hold", echo = FALSE)

collapsible_script_section <- function(
    file_path, button_label, section_id, 
    button_class = "ShortBlockButton", 
    div_class = "ShortBlock") {
  if (!file.exists(file_path)) {
    script_content <- "Error: File not found."
  } else {
    script_content <- readLines(file_path) %>% 
      paste(collapse = "\n") %>% 
      htmltools::htmlEscape()
  }
  html_output <- paste0(
    '<button onclick="toggleScript(\'', section_id, '\')" class="', button_class, '">', button_label, '</button>\n',
    '<div id="', section_id, '" class="', div_class, '" style="display: none;">\n',
    '  <pre><code class="bash">', script_content, '</code></pre>\n',
    '</div>\n'
  )
  cat(html_output)
}

## -----------------------------------------------------------------------------
c("<style>", readLines("style.css"), "</style>") %>% 
  cat(sep = "\n")

## -----------------------------------------------------------------------------
DT <- tibble::tribble(
  ~R_function, ~Description,
  "`mod_slurm_refit()`", "checks for unsuccessful model fits.",
  "`mod_merge_chains()`", "merges MCMC chains and saves the fitted model and coda objects to `.qs2` or `.RData` files.",
  "`convergence_plot_all()`", "visualises the convergence of `rho`, `alpha`, `omega`, and `beta` parameters across all model variants. This function is particularly useful for comparing convergence across models with different thinning values, with and without phylogenetic relationships, or varying GPP knot distances. It is unnecessary when only a single model variant is considered.",
  "`convergence_plot()`", "displays the convergence of `rho`, `alpha`, `omega`, and `beta` parameters for a selected model, providing a more detailed view compared to `convergence_plot_all()`.",
  "`plot_gelman()`", "visualises the Gelman-Rubin-Brooks diagnostics for the selected model.",
  "`mod_summary()`", "extracts and saves a summary of the model.",
  "`mod_heatmap_beta()`", "generates heatmaps of the `beta` parameters.",
  "`mod_heatmap_omega()`", "generates heatmaps of the `omega` parameter, which represents residual species associations.",
  "`mod_cv_fit()`", "prepares the necessary data for fitting cross-validated models. 
  <ul><li>output files are saved in the `model_fitting_cv` subdirectory.</li><li>the type of cross-validation strategy is controlled by the `cv_name` argument, which defaults to both `cv_dist` and `cv_large`.</li><li>unfitted model objects are saved in the `model_init` subdirectory.</li><li>commands for model fitting are saved as text files, with a separate file for each cross-validation strategy (e.g., `commands_to_fit_cv_dist.txt`, `commands_to_fit_cv_large.txt`).</li><li>model fitting commands are submitted as batch jobs using SLURM scripts, with a separate script for each strategy (e.g., `cv_bash_fit_dist.slurm`, `cv_bash_fit_large.slurm`).</li></ul>") %>% 
  dplyr::mutate(R_function = paste0("<cc>", R_function, "</cc><tab0>"))

DT %>% 
  knitr::kable(col.names = NULL, format = "html", escape = FALSE)  %>%
  kableExtra::kable_styling(
    c("basic", "condensed", "hover"), 
    font_size = 16, full_width = TRUE) %>%
  kableExtra::add_indent(
    positions = seq_len(nrow(DT)), level_of_indent = 2) %>% 
  kableExtra::column_spec(
    column = 1, extra_css = "white-space: nowrap;")

## -----------------------------------------------------------------------------
DT <- tibble::tribble(
  ~R_function, ~Description,
  "`rc_prepare_data()`", "prepares data for predicting latent factors for response curves",
  "`predict_maps()`", "prepares data for predicting latent factors at new sampling units",
  "`variance_partitioning_compute()`", "prepares data for computing variance partitioning") %>% 
  dplyr::mutate(R_function = paste0("<cc>", R_function, "</cc><tab0>"))

DT %>% 
  knitr::kable(col.names = NULL, format = "html", escape = FALSE)  %>%
  kableExtra::kable_styling(
    c("basic", "condensed", "hover"),
    font_size = 16, full_width = TRUE) %>%
  kableExtra::add_indent(
    positions = seq_len(nrow(DT)), level_of_indent = 2) %>% 
  kableExtra::column_spec(
    column = 1, extra_css = "white-space: nowrap;")

## -----------------------------------------------------------------------------
collapsible_script_section(
  "Example_lf_rc_commands.txt", 
  "&#187;&#187; Example lf_rc_commands.txt file", 
  "scriptContentNS1")

## -----------------------------------------------------------------------------
collapsible_script_section(
  "Example_lf_new_sites_Commands.txt", 
  "&#187;&#187; Example lf_new_sites_Commands.txt file", 
  "scriptContentNS2")

## -----------------------------------------------------------------------------
collapsible_script_section(
  "Example_TF_Chunk.txt", 
  "&#187;&#187; Example tf_chunk_*.txt file", 
  "scriptContentNS3")

## -----------------------------------------------------------------------------
collapsible_script_section(
  "Example_lf_SLURM.slurm", 
  "&#187;&#187; Example lf_SLURM.slurm file", 
  "scriptContentNS4")

## -----------------------------------------------------------------------------
collapsible_script_section(
  "Example_VP_Commands.txt", 
  "&#187;&#187; Example VP_Commands.txt file", 
  "scriptContentNS5")

## -----------------------------------------------------------------------------
collapsible_script_section(
  "Example_VP_SLURM.slurm", 
  "&#187;&#187; Example VP_SLURM.slurm file", 
  "scriptContentNS6")

## -----------------------------------------------------------------------------
DT <- tibble::tribble(
  ~R_function, ~Description,
  "`rc_prepare_data()`,<br/>`rc_plot_sr()`,<br/>`rc_plot_species()`,<br/>`rc_plot_species_all()`", "continues the processing and visualisation of response curves.",
  "`predict_maps()`", "<ul><li>predicts habitat suitability across various climate scenarios.</li><li>computes the model's explanatory power (internal evaluation without cross-validation) using four metrics: AUC (area under the ROC curve), RMSE (root mean square error), continuous Boyce index, and Tjur R<sup>2</sup>.</li><li>prepares GeoTIFF maps for free access via the `IASDT` <a href='http://opendap.biodt.eu/ias-pdt/' target='_blank' class='ll'>OPeNDAP server</a> and the `IASDT` <a href='https://app.biodt.eu/app/biodtshiny' target='_blank' class='ll'>Shiny app</a>.</li></ul>",
  "`plot_prediction()`", "visualises predictions of species and species richness as JPEG images.",
  "`plot_latent_factor()`", "visualises the spatial variation in site loadings of HMSC models as JPEG images.",
  "`variance_partitioning_compute()`,<br/>`variance_partitioning_plot()`", "continues the processing and visualisation of variance partitioning.",
  "`plot_evaluation()`", "visualises the explanatory power of the model (internal evaluation).") %>% 
  dplyr::mutate(R_function = paste0("<cc>", R_function, "</cc><tab0>"))

DT %>% 
  knitr::kable(col.names = NULL, format = "html", escape = FALSE)  %>%
  kableExtra::kable_styling(
    c("basic", "condensed", "hover"), 
    font_size = 16, full_width = TRUE) %>%
  kableExtra::add_indent(
    positions = seq_len(nrow(DT)), level_of_indent = 1) %>% 
  kableExtra::column_spec(
    column = 1, extra_css = "white-space: nowrap;")

## -----------------------------------------------------------------------------
DT <- tibble::tribble(
  ~R_function, ~Description,
  "`mod_merge_chains_cv()`", "merges fitted cross-validated model chains into `Hmsc` model objects and saves them to disk.",
  "`predict_maps_cv()`", "prepares scripts for predicting latent factors for each cross-validation strategy at new sampling units (evaluation folds). The arguments `lf_only` and `lf_commands_only` are set to `TRUE` to prepare only the necessary script files.") %>% 
  dplyr::mutate(R_function = paste0("<cc>", R_function, "</cc><tab0>"))

DT %>% 
  knitr::kable(col.names = NULL, format = "html", escape = FALSE)  %>%
  kableExtra::kable_styling(
    c("basic", "condensed", "hover"), 
    font_size = 16, full_width = TRUE) %>%
  kableExtra::add_indent(
    positions = seq_len(nrow(DT)), level_of_indent = 1) %>% 
  kableExtra::column_spec(
    column = 1, extra_css = "white-space: nowrap;")

