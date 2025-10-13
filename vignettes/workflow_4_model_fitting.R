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
  ~Argument, ~Description,
  "hab_abb", "abbreviation of a single habitat type to be modelled",
  "directory_name", "directory path for storing all model files",
  "min_efforts_n_species", 'minimum number of vascular plant species per grid cell required for inclusion in model fitting. This reflects the total count of vascular plant species (including native species) recorded in GBIF across Europe, as computed during the <a href="workflow_2_abiotic_data.html#sampling-efforts" class="ll">sampling effort preparation</a> step (`efforts_process()`). This argument filters out grid cells with insufficient sampling effort',
  "exclude_cultivated", "whether to exclude countries with cultivated or casual observations for each species",
  "exclude_0_habitat", 'whether to exclude grid cells with zero <a href="workflow_2_abiotic_data.html#corine-land-cover-habitat-data" class="ll">habitat coverage</a> of the respective habitat type',
  "n_pres_per_species", "minimum number of presence grid cells required for a species to be included in the models, calculated after excluding grid cells with low sampling effort (<cc>min_efforts_n_species</cc>), zero habitat coverage (<cc>exclude_0_habitat</cc>), and countries with cultivated or casual observations (<cc>exclude_cultivated</cc>)") %>% 
  dplyr::mutate(Argument = paste0("<cc>", Argument, "</cc><tab0>"))

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
  ~Argument, ~Description,
  "cv_n_folds", "number of cross-validation folds",
  "cv_n_grids", "number of grid cells in each directions for the <i>cv_dist</i> cross-validation strategy (default: 20, yielding 20 &times; 20 grid cell blocks).",
  
  "cv_n_rows / cv_n_columns", "number of rows and columns defining in the <i>cv_large</i> cross-validation strategy, partitioning the study area into large blocks (default: *cv_n_rows = cv_n_columns = 2*, resulting in four blocks divided at median coordinates).") %>% 
  dplyr::mutate(Argument = paste0("<cc>", Argument, "</cc><tab0>"))
  
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
  ~Argument, ~Description,
  "gpp", "whether to incorporate spatial random effects using the Gaussian Predictive Process (GPP)",
  "gpp_dists", "distance (in kilometres; controlled by the `min_distance` argument of `prepare_knots()`) specifying both the spacing between knots and the minimum distance between a knot and the nearest sampling point",
  "gpp_plot", "whether to plot the coordinates of sampling units and knots",
  "min_lf / max_lf", "minimum and maximum number of latent factors to be include",
  "alphapw", "prior specification for the alpha parameter") %>% 
  dplyr::mutate(Argument = paste0("<cc>", Argument, "</cc><tab0>"))

DT %>% 
  knitr::kable(col.names = NULL, format = "html", escape = FALSE)  %>%
  kableExtra::kable_styling(
    c("basic", "condensed", "hover"), font_size = 16, full_width = TRUE) %>%
  add_indent(positions = seq_len(nrow(DT)), level_of_indent = 1) %>% 
  kableExtra::column_spec(
    column = 1, extra_css = "white-space: nowrap;")

## -----------------------------------------------------------------------------
DT <- tibble::tribble(
  ~Argument, ~Description,
  "job_name", "name assigned to the SLURM job",
  "ntasks", "number of tasks to execute",
  "cpus_per_task / gpus_per_node", "Number of CPUs and GPUs allocated per node",
  "memory_per_cpu", "memory allocation per CPU",
  "job_runtime", "maximum duration for job execution",
  "hpc_partition", "name of the HPC partition",
  "n_array_jobs", "number of jobs within each SLURM script") %>% 
  dplyr::mutate(Argument = paste0("<cc>", Argument, "</cc><tab0>"))

DT %>% 
  knitr::kable(col.names = NULL, format = "html", escape = FALSE)  %>%
  kableExtra::kable_styling(
    c("basic", "condensed", "hover"), font_size = 16, full_width = TRUE) %>%
  kableExtra::add_indent(
    positions = seq_len(nrow(DT)), level_of_indent = 1) %>% 
  kableExtra::column_spec(
    column = 1, extra_css = "white-space: nowrap;")

## -----------------------------------------------------------------------------
DT <- tibble::tribble(
  ~Argument, ~Description,
  "bio_variables", 'names of <a href="workflow_2_abiotic_data.html#chelsa-climate-data" class="ll">CHELSA</a> variables to include in the model',
  "quadratic_variables", "names of variables for which quadratic terms are incorporated",
  "efforts_as_predictor", 'whether to include the (log<sub>10</sub>-transformed) <a href="workflow_2_abiotic_data.html#sampling-efforts" class="ll">sampling effort</a> as a predictor',
  "road_rail_as_predictor", 'whether to include the (log<sub>10</sub>-transformed) <a href="workflow_2_abiotic_data.html#railways-and-roads-intensity" class="ll">summed road and railway intensity</a> as a predictor',
  "habitat_as_predictor", 'whether to include the (log<sub>10</sub>-transformed) <a href="workflow_2_abiotic_data.html#corine-land-cover-habitat-data" class="ll">percentage coverage of the respective habitat type</a> per grid cell as a predictor',
  "river_as_predictor", 'whether to include the (log<sub>10</sub>-transformed) total <a href="workflow_2_abiotic_data.html#river-length" class="ll">river length</a> per grid cell as a predictor') %>% 
  dplyr::mutate(Argument = paste0("<cc>", Argument, "</cc><tab0>"))

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
  ~Argument, ~Description,
  "mcmc_n_chains", "number of MCMC chains",
  "mcmc_thin", "thinning value(s) in MCMC sampling",
  "mcmc_samples", "number of MCMC samples per chain",
  "mcmc_transient_factor", "transient multiplication factor. The value of transient will equal the multiplication of <cc>transientFactor</cc> and <cc>thin</cc>",
  "mcmc_verbose", "interval at which MCMC sampling progress is reported",
  "precision", "floating-point precision mode for Hmsc-HPC sampling") %>% 
  dplyr::mutate(Argument = paste0("<cc>", Argument, "</cc><tab0>"))

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
collapsible_script_section(
  "Example_python.txt", 
  "&#187;&#187; Example model fitting commands", 
  "scriptContentNS1")

## -----------------------------------------------------------------------------
collapsible_script_section(
  "Example_SLURM.txt", 
  "&#187;&#187; Example SLURM script", 
  "scriptContentNS2")

