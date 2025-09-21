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
  ~Abbreviation, ~`Habitat Type`, ~Description,
  "1", "Forests", "closed vegetation dominated by deciduous or evergreen trees",
  "2", "Open forests", "woodlands with canopy openings created by environmental stress or disturbance, including forest edges",
  "3", "Scrub", "shrublands maintained by environmental stress (aridity) or disturbance",
  "4a", "Natural grasslands", "grasslands maintained by climate (aridity, unevenly distributed precipitation), herbivores or environmental stress (aridity, instability or toxicity of substrate)",
  "4b", "Human-maintained grasslands", "grasslands dependent on regular human-induced management (mowing, grazing by livestock, artificial burning)",
  "10", "Wetlands", "sites with the permanent or seasonal influence of moisture, ranging from oligotrophic to eutrophic",
  "12a", "Ruderal habitats", "anthropogenically disturbed or eutrophicated sites, where the anthropogenic disturbance or fertilization is typically a side-product and not the aim of the management",
  "12b", "Agricultural habitats", "synanthropic habitats directly associated with growing of agricultural products, thus dependent on specific type of management (ploughing, fertilization)")

DT %>% 
  knitr::kable(format = "html")  %>%
  kableExtra::kable_styling(
    c("basic", "condensed", "hover"), font_size = 16, full_width = TRUE)

## -----------------------------------------------------------------------------
DT <- tibble::tribble(
  ~Variable, ~Category, ~Description, ~"Default value",
  "DP_R_bioreg_interim", "Biogeographical regions", "Directory path to biogeographical regions interim data", "datasets/interim/biogeoregions",
  "DP_R_bioreg_processed", "Biogeographical regions", "Directory path to biogeographical regions processed data", "datasets/processed/biogeoregions",
  "DP_R_bioreg_raw", "Biogeographical regions", "Directory path to biogeographical regions raw data", "datasets/raw/biogeoregions",
  "DP_R_bioreg_url", "Biogeographical regions", "URL for downloading biogeographical regions data", "https://www.eea.europa.eu/en/datahub/datahubitem-view/11db8d14-f167-4cd5-9205-95638dfd9618",
  "DP_R_chelsa_links", "CHELSA", "Directory path containing CHELSA download links", "references/chelsa/dwnload_links",
  "DP_R_chelsa_processed", "CHELSA", "Directory path to processed CHELSA data", "datasets/processed/chelsa",
  "DP_R_chelsa_raw", "CHELSA", "Directory path to raw CHELSA data", "datasets/raw/chelsa",
  "DP_R_chelsa_url", "CHELSA", "Base URL for CHELSA data", "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL",
  "DP_R_clc_crosswalk", "CLC", "Path to a text file containing custom cross-walk between CLC values at level 3 and their corresponding values for *EUNIS* and *SynHab* habitat types", "references/CrossWalk.txt",
  "DP_R_clc_processed", "CLC", "Directory path to processed CLC data", "datasets/processed/corine",
  "DP_R_clc_tif", "CLC", "Path to the input CLC tif file", "datasets/raw/corine/u2018_clc2018_v2020_20u1_raster100m/DATA/U2018_CLC2018_V2020_20u1.tif",
  "DP_R_lumi_cpu", "HPC", "LUMI project number for CPU computations", "project_465001857",
  "DP_R_lumi_gpu", "HPC", "LUMI project number for GPU computations", "project_465001857",
  "DP_R_lumi_gpu_check", "HPC", "File path to a python script for reporting if the GPU was used in the running SLURM job", "references/LUMI_Check_GPU.py",
  "DP_R_country_codes", "Misc", "Path to a file containing countries ISO codes", "references/country_codes.csv",
  "DP_R_country_boundaries", "Misc", "Path to `RData` file containing country boundaries", "references/EU_boundaries_sf.RData",
  "DP_R_model_root_path", "Models", "Directory path for model fitting", "datasets/processed/model_fitting",
  "DP_R_railway_interim", "Railways", "Directory path to interim railways data", "datasets/interim/railway",
  "DP_R_railway_processed", "Railways", "Directory path to processed railways data", "datasets/processed/railway",
  "DP_R_railway_raw", "Railways", "Directory path to raw railways data", "datasets/raw/railway",
  "DP_R_railway_url", "Railways", "URL for railways data", "https://download.geofabrik.de/",
  "DP_R_grid_processed", "Reference grid", "Directory path for reference grid (resulted from processing CLC data)", "datasets/processed/grid",
  "DP_R_grid_raw", "Reference grid", "Directory path for reference grid (original)", "references/grid",
  "DP_R_rivers_interim", "Rivers", "Directory path to interim rivers data", "datasets/interim/rivers",
  "DP_R_rivers_processed", "Rivers", "Directory path to processed rivers data", "datasets/processed/rivers",
  "DP_R_rivers_raw", "Rivers", "Directory path to raw rivers data", "datasets/raw/rivers",
  "DP_R_rivers_zip", "Rivers", "Path to zip file containing river data", "datasets/raw/rivers/EU_hydro_gpkg_eu.zip",
  "DP_R_roads_interim", "Roads", "Directory path to interim roads data", "datasets/interim/roads",
  "DP_R_roads_processed", "Roads", "Directory path to processed roads data", "datasets/processed/roads",
  "DP_R_roads_raw", "Roads", "Directory path to raw roads data", "datasets/raw/roads",
  "DP_R_roads_url", "Roads", "URL for the Global Roads Inventory Project (GRIP) data", "https://dataportaal.pbl.nl/downloads/GRIP4/GRIP4_Region4_vector_fgdb.zip",
  "DP_R_efforts_interim", "Sampling efforts", "Directory path to interim sampling efforts data", "datasets/interim/sampling_efforts",
  "DP_R_efforts_processed", "Sampling efforts", "Directory path to processed sampling efforts data", "datasets/processed/sampling_efforts",
  "DP_R_efforts_raw", "Sampling efforts", "Directory path to raw sampling efforts data", "datasets/raw/sampling_efforts",
  "DP_R_easin_interim", "Species distribution", "Directory path to EASIN data", "datasets/interim/EASIN",
  "DP_R_easin_processed", "Species distribution", "Directory path to processed EASIN data", "datasets/processed/EASIN",
  "DP_R_easin_summary", "Species distribution", "Directory path to summary of processed EASIN data", "datasets/processed/EASIN/Summary",
  "DP_R_easin_taxa_url", "Species distribution", "URL for EASIN API for downloading taxa list", "https://easin.jrc.ec.europa.eu/apixg/catxg",
  "DP_R_easin_data_url", "Species distribution", "URL for EASIN API for downloading species data", "https://easin.jrc.ec.europa.eu/apixg/geoxg",
  "DP_R_elter_processed", "Species distribution", "Path to `RData` containing processed eLTER presence-absence data", "datasets/processed/naps_pa/elter_naps.RData",
  "DP_R_elter_raw", "Species distribution", "Path to `rds` file containing processed and standardized eLTER data", "references/elter_data_gbif_2024-02-07.rds",
  "DP_R_gbif_interim", "Species distribution", "Directory path to interim GBIF data", "datasets/interim/GBIF",
  "DP_R_gbif_processed", "Species distribution", "Directory path to processed GBIF data", "datasets/processed/GBIF",
  "DP_R_gbif_raw", "Species distribution", "Directory path to raw GBIF data", "datasets/raw/gbif",
  "DP_R_pa", "Species distribution", "Directory path to species-specific presence-absence data", "datasets/processed/naps_pa",
  "DP_R_hab_affinity", "Taxa information", "Path to `rds` file containing species affinity to habitat types", "references/taxon-habitat-combined_2024-02-07.rds",
  "DP_R_taxa_country", "Taxa information", "Path to excel file containing the number of grid cells per species and country", "references/cell_count_per_species_and_country_2024-01-30-IA-VK-MG.xlsx",
  "DP_R_taxa_easin", "Taxa information", "Path to `rds` file containing EASIN standardized taxonomy", "references/easin_taxon-list-gbif_2024-02-07.rds",
  "DP_R_opendap_url", "", "Path to OpENDAP server for the `IASDT`", "http://opendap.biodt.eu/ias-pdt/") %>% 
  dplyr::mutate(
    Variable = paste0("<cc>", Variable, "</cc>"),
    Category = NULL,
    `Default value` = paste0("<i>", `Default value`, "</i>"))

DT %>% 
  knitr::kable(format = "html", escape = FALSE)  %>%
  kableExtra::kable_styling(
    c("basic", "condensed", "hover"), font_size = 16, full_width = TRUE) %>% 
  kableExtra::column_spec(column = 1, extra_css = "white-space: nowrap;") %>% 
  kableExtra::column_spec(column = 2, width = "25em")


