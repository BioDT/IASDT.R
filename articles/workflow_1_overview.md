# IASDT modelling workflow — 1. Overview

## Overview

This document delineates the workflow for modelling the distribution and
level of invasion (species richness) of naturalised alien plant species
(NAPS) across Europe. The invasive alien species (IAS) Digital Twin
(`IASDT`) is a component of the European [BioDT](https://www.biodt.eu/)
project, which seeks to establish a Digital Twin framework for
biodiversity in Europe. For a detailed exposition of the `IASDT`, refer
to Khan, El-Gabbas, *et al.*
([2024](https://doi.org/10.3897/rio.10.e124579)). The complete `IASDT`
workflow is documented at
[Zenodo](https://doi.org/10.5281/zenodo.14756907).

The `IASDT` leverages the
[`IASDT.R`](https://doi.org/10.5281/zenodo.14834384) and
[`ecokit`](https://doi.org/10.5281/zenodo.15477683) R packages to
execute a comprehensive workflow encompassing data processing, model
fitting, post-processing, and preparation maps for the IASDT
[Shiny](https://app.biodt.eu/app/biodtshiny) dashboard. The `IASDT.R`
package facilitates the preparation of abiotic data (e.g., climate and
land cover) and biotic data (i.e., species distribution). Model outputs
from the `IASDT` are made publicly available to end-users and
stakeholders through an [OPeNDAP cloud
server](http://opendap.biodt.eu/ias-pdt/).

------------------------------------------------------------------------

## Models

Species distribution models are constructed using the Hierarchical
Modelling of Species Communities
([HMSC](https://www.helsinki.fi/en/researchgroups/statistical-ecology/software/hmsc))
R package, a hierarchical Bayesian framework that incorporates spatial
autocorrelation and species associations. Spatial autocorrelation is
modelled via the Gaussian Predictive Process (GPP; Tikhonov *et al.*,
[2019](https://doi.org/10.1002/ecy.2929)), offering a flexible and
computationally efficient approach to capturing spatial dependencies.
Given the substantial computational demands of fitting these spatial
models at a European scale, we utilise the `HMSC-HPC` extension (Rahman
*et al.*, [2024](https://doi.org/10.1371/journal.pcbi.1011914)) to
leverage GPU-based processing for enhanced efficiency.

Models are fitted at the habitat level, with a distinct model fitted for
each habitat type, incorporating only those invasive alien species (IAS)
associated with the respective habitat type. We employ the habitat
classification delineated by Pyšek *et al.*
([2022](https://doi.org/10.23855/preslia.2022.447)). We fitted the
models at eight habitat types (see table below). For each habitat type,
predictions of individual species distributions and species richness are
generated across multiple climate scenarios; further details are
provided in the [abiotic data processing
section](https://biodt.github.io/IASDT.R/articles/workflow_2_abiotic_data.html#chelsa-climate-data).
Model performance is assessed using spatial block cross-validation to
ensure spatial independence between training and testing datasets.

| Abbreviation | Habitat Type                | Description                                                                                                                                                            |
|:-------------|:----------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 1            | Forests                     | closed vegetation dominated by deciduous or evergreen trees                                                                                                            |
| 2            | Open forests                | woodlands with canopy openings created by environmental stress or disturbance, including forest edges                                                                  |
| 3            | Scrub                       | shrublands maintained by environmental stress (aridity) or disturbance                                                                                                 |
| 4a           | Natural grasslands          | grasslands maintained by climate (aridity, unevenly distributed precipitation), herbivores or environmental stress (aridity, instability or toxicity of substrate)     |
| 4b           | Human-maintained grasslands | grasslands dependent on regular human-induced management (mowing, grazing by livestock, artificial burning)                                                            |
| 10           | Wetlands                    | sites with the permanent or seasonal influence of moisture, ranging from oligotrophic to eutrophic                                                                     |
| 12a          | Ruderal habitats            | anthropogenically disturbed or eutrophicated sites, where the anthropogenic disturbance or fertilization is typically a side-product and not the aim of the management |
| 12b          | Agricultural habitats       | synanthropic habitats directly associated with growing of agricultural products, thus dependent on specific type of management (ploughing, fertilization)              |

------------------------------------------------------------------------

## Environment variables

The workflow necessitates the configuration of multiple environment
variables to ensure proper execution. Certain functions within the
`IASDT.R` package include an `env_file` argument, which defaults to
`.env`. The table below enumerates the environment variables essential
to the workflow, accompanied by their descriptions and default values.

  

| Variable                | Description                                                                                                                                          | Default value                                                                                  |
|:------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------|
| DP_R_bioreg_interim     | Directory path to biogeographical regions interim data                                                                                               | *datasets/interim/biogeoregions*                                                               |
| DP_R_bioreg_processed   | Directory path to biogeographical regions processed data                                                                                             | *datasets/processed/biogeoregions*                                                             |
| DP_R_bioreg_raw         | Directory path to biogeographical regions raw data                                                                                                   | *datasets/raw/biogeoregions*                                                                   |
| DP_R_bioreg_url         | URL for downloading biogeographical regions data                                                                                                     | *<https://www.eea.europa.eu/en/datahub/datahubitem-view/11db8d14-f167-4cd5-9205-95638dfd9618>* |
| DP_R_chelsa_links       | Directory path containing CHELSA download links                                                                                                      | *references/chelsa/dwnload_links*                                                              |
| DP_R_chelsa_processed   | Directory path to processed CHELSA data                                                                                                              | *datasets/processed/chelsa*                                                                    |
| DP_R_chelsa_raw         | Directory path to raw CHELSA data                                                                                                                    | *datasets/raw/chelsa*                                                                          |
| DP_R_chelsa_url         | Base URL for CHELSA data                                                                                                                             | *<https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL>*                                            |
| DP_R_clc_crosswalk      | Path to a text file containing custom cross-walk between CLC values at level 3 and their corresponding values for *EUNIS* and *SynHab* habitat types | *references/clc_crosswalk.txt*                                                                 |
| DP_R_clc_processed      | Directory path to processed CLC data                                                                                                                 | *datasets/processed/corine*                                                                    |
| DP_R_clc_tif            | Path to the input CLC tif file                                                                                                                       | *datasets/raw/corine/u2018_clc2018_v2020_20u1_raster100m/DATA/U2018_CLC2018_V2020_20u1.tif*    |
| DP_R_lumi_cpu           | LUMI project number for CPU computations                                                                                                             | *project_XXXXXXXXX*                                                                            |
| DP_R_lumi_gpu           | LUMI project number for GPU computations                                                                                                             | *project_XXXXXXXXX*                                                                            |
| DP_R_lumi_gpu_check     | File path to a python script for reporting if the GPU was used in the running SLURM job                                                              | *references/lumi_check_gpu.py*                                                                 |
| DP_R_country_codes      | Path to a file containing countries ISO codes                                                                                                        | *references/country_codes.csv*                                                                 |
| DP_R_country_boundaries | Path to `RData` file containing country boundaries                                                                                                   | *references/EU_boundaries_sf.RData*                                                            |
| DP_R_model_root_path    | Directory path for model fitting                                                                                                                     | *datasets/processed/model_fitting*                                                             |
| DP_R_railway_interim    | Directory path to interim railways data                                                                                                              | *datasets/interim/railways*                                                                    |
| DP_R_railway_processed  | Directory path to processed railways data                                                                                                            | *datasets/processed/railways*                                                                  |
| DP_R_railway_raw        | Directory path to raw railways data                                                                                                                  | *datasets/raw/railways*                                                                        |
| DP_R_railway_url        | URL for railways data                                                                                                                                | *<https://download.geofabrik.de/>*                                                             |
| DP_R_grid_processed     | Directory path for reference grid (resulted from processing CLC data)                                                                                | *datasets/processed/grid*                                                                      |
| DP_R_grid_raw           | Directory path for reference grid (original)                                                                                                         | *references/grid*                                                                              |
| DP_R_rivers_interim     | Directory path to interim rivers data                                                                                                                | *datasets/interim/rivers*                                                                      |
| DP_R_rivers_processed   | Directory path to processed rivers data                                                                                                              | *datasets/processed/rivers*                                                                    |
| DP_R_rivers_raw         | Directory path to raw rivers data                                                                                                                    | *datasets/raw/rivers*                                                                          |
| DP_R_rivers_zip         | Path to zip file containing river data                                                                                                               | *datasets/raw/rivers/EU_hydro_gpkg_eu.zip*                                                     |
| DP_R_roads_interim      | Directory path to interim roads data                                                                                                                 | *datasets/interim/roads*                                                                       |
| DP_R_roads_processed    | Directory path to processed roads data                                                                                                               | *datasets/processed/roads*                                                                     |
| DP_R_roads_raw          | Directory path to raw roads data                                                                                                                     | *datasets/raw/roads*                                                                           |
| DP_R_roads_url          | URL for the Global Roads Inventory Project (GRIP) data                                                                                               | *<https://dataportaal.pbl.nl/downloads/GRIP4/GRIP4_Region4_vector_fgdb.zip>*                   |
| DP_R_efforts_interim    | Directory path to interim sampling efforts data                                                                                                      | *datasets/interim/sampling_efforts*                                                            |
| DP_R_efforts_processed  | Directory path to processed sampling efforts data                                                                                                    | *datasets/processed/sampling_efforts*                                                          |
| DP_R_efforts_raw        | Directory path to raw sampling efforts data                                                                                                          | *datasets/raw/sampling_efforts*                                                                |
| DP_R_easin_interim      | Directory path to EASIN data                                                                                                                         | *datasets/interim/easin*                                                                       |
| DP_R_easin_processed    | Directory path to processed EASIN data                                                                                                               | *datasets/processed/easin*                                                                     |
| DP_R_easin_summary      | Directory path to summary of processed EASIN data                                                                                                    | *datasets/processed/easin/Summary*                                                             |
| DP_R_easin_taxa_url     | URL for EASIN API for downloading taxa list                                                                                                          | *<https://easin.jrc.ec.europa.eu/apixg/catxg>*                                                 |
| DP_R_easin_data_url     | URL for EASIN API for downloading species data                                                                                                       | *<https://easin.jrc.ec.europa.eu/apixg/geoxg>*                                                 |
| DP_R_elter_processed    | Path to `RData` containing processed eLTER presence-absence data                                                                                     | *datasets/processed/naps_pa/elter_naps.RData*                                                  |
| DP_R_elter_raw          | Path to `rds` file containing processed and standardized eLTER data                                                                                  | *references/elter_data_gbif_202X-XX-XX7.rds*                                                   |
| DP_R_gbif_interim       | Directory path to interim GBIF data                                                                                                                  | *datasets/interim/gbif*                                                                        |
| DP_R_gbif_processed     | Directory path to processed GBIF data                                                                                                                | *datasets/processed/gbif*                                                                      |
| DP_R_gbif_raw           | Directory path to raw GBIF data                                                                                                                      | *datasets/raw/gbif*                                                                            |
| DP_R_pa                 | Directory path to species-specific presence-absence data                                                                                             | *datasets/processed/naps_pa*                                                                   |
| DP_R_hab_affinity       | Path to `rds` file containing species affinity to habitat types                                                                                      | *references/taxon-habitat-combined_202X-XX-XX.rds*                                             |
| DP_R_taxa_country       | Path to excel file containing the number of grid cells per species and country                                                                       | *references/2_taxon-status_202X-XX-XX.xlsx*                                                    |
| DP_R_taxa_easin         | Path to `rds` file containing EASIN standardized taxonomy                                                                                            | *references/easin_taxon-list-gbif_202X-XX-XX.rds*                                              |
| DP_R_opendap_url        | Path to OpENDAP server for the `IASDT`                                                                                                               | *<http://opendap.biodt.eu/ias-pdt/>*                                                           |

------------------------------------------------------------------------

**Next articles:**  
↠[2. Processing abiotic
data](https://biodt.github.io/IASDT.R/articles/workflow_2_abiotic_data.md)  
↠[3. Processing biotic
data](https://biodt.github.io/IASDT.R/articles/workflow_3_biotic_data.md)  
↠[4. Model
fitting](https://biodt.github.io/IASDT.R/articles/workflow_4_model_fitting.md)  
↠[5. Model
post-processing](https://biodt.github.io/IASDT.R/articles/workflow_5_model_postprocess.md)  
