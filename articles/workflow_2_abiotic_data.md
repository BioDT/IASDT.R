# IASDT modelling workflow — 2. abiotic data

  
  

This article details the processing of abiotic data within the `IASDT`
modelling workflow. These processed data serve as predictor variables in
the species distribution models. All data preparation adheres to the
[FAIR](https://doi.org/10.1038/sdata.2016.18) principles (Findable,
Accessible, Interoperable, and Reusable) to ensure scientific integrity
and reproducibility.

  
  

## Reference grid

The workflow employs the European Environment Agency (EEA) [reference
grid](https://www.eea.europa.eu/en/datahub/datahubitem-view/3c362237-daa4-45e2-8c16-aaadfb1a003b),
standardized at a 10×10 km resolution across the study area. This grid
utilises the ETRS89-LAEA Europe coordinate reference system (CRS;
`EPSG:3035`).

  
  

## Corine land cover and habitat data

The
[`clc_process()`](https://biodt.github.io/IASDT.R/reference/CLC_Process.md)
function manages the processing of Corine Land Cover
([CLC](https://land.copernicus.eu/en/products/corine-land-cover/clc2018))
data within the `IASDT` workflow. It computes the percentage coverage
and predominant classes per grid cell across all three CLC levels,
alongside `EUNIS` and `SynHab` habitat classifications. A custom
crosswalk was used to transform the CLC Level 3 data into ecologically
meaningful `EUNIS` and `SynHab` habitat classes. The resulting data
serve the following purposes:

- **model grid selection**: Grid cells with at least 15% land cover
  (default threshold) are retained for modelling.
- **habitat-specific modelling**: The percentage coverage of `SynHab`
  habitat types informs model fitting by:
  1.  optionally excluding grid cells with zero coverage of the relevant
      habitat type during model fitting.
  2.  serving as a potential predictor variable in the models.

  
  

## Biogeographical regions

The
[`bioreg_process()`](https://biodt.github.io/IASDT.R/reference/BioReg_Process.md)
function retrieves and processes the biogeographical regions dataset
from the
[EEA](https://www.eea.europa.eu/en/datahub/datahubitem-view/11db8d14-f167-4cd5-9205-95638dfd9618).
It extracts the names of biogeographical regions corresponding to each
reference grid cell, enabling the quantification of species presence
across these regions.

  
  

## CHELSA climate data

The
[`chelsa_process()`](https://biodt.github.io/IASDT.R/reference/CHELSA_data.md)
function manages the retrieval and processing of
[CHELSA](https://chelsa-climate.org/) (Climatologies at High Resolution
for the Earth’s Land Surface Areas) climate data for the study area
across multiple climate scenarios. CHELSA delivers high-resolution
global datasets encompassing a range of environmental variables for
current conditions and future projections. These processed data are
integrated into the species distribution models as predictor variables.
For each environmental variable, the dataset encompasses 45 future
scenarios, derived from the combination of five CMIP6 climate models,
three Shared Socioeconomic Pathways (SSPs), and three future time
periods (refer to the CHELSA [technical
specifications](https://chelsa-climate.org/wp-admin/download-page/chelsa_tech_specification_V2.pdf)
for details). Additionally, future climate model outputs are aggregated
to generate ensemble predictions for each SSP and time period
combination.

| Climate model | Institution                                          |
|:--------------|:-----------------------------------------------------|
| mpi-esm1-2-hr | Max Planck Institute for Meteorology, Germany        |
| ipsl-cm6a-lr  | Institut Pierre Simon Laplace, France                |
| ukesm1-0-ll   | Met Office Hadley Centre, UK                         |
| gfdl-esm4     | National Oceanic and Atmospheric Administration, USA |
| mri-esm2-0    | Meteorological Research Institute, Japan             |

| Shared Socioeconomic Pathway | Description                                  |
|:-----------------------------|:---------------------------------------------|
| ssp126                       | SSP1-RCP2.6 climate as simulated by the GCMs |
| ssp370                       | SSP3-RCP7 climate as simulated by the GCMs   |
| ssp585                       | SSP5-RCP8.5 climate as simulated by the GCMs |

  
  

## Railways and roads intensity

The
[`railway_intensity()`](https://biodt.github.io/IASDT.R/reference/Railway_Intensity.md)
and
[`road_intensity()`](https://biodt.github.io/IASDT.R/reference/Road_Intensity.md)
functions retrieve and process railway data (sourced from
[OpenRailwayMap](https://www.openrailwaymap.org/)) and road data
(sourced from the Global Roads Inventory Project;
[GRIP](https://www.globio.info/download-grip-dataset)), respectively,
for the study area. Road intensity reflects site accessibility, habitat
disturbance levels, and IAS dispersal potential, while railway density
serves as a proxy for IAS dispersal routes. The summed lengths of
railways and roads per grid cell, transformed to a logarithmic scale
(log₁₀), are incorporated into the models as predictor variables.

  
  

## River length

The
[`river_length()`](https://biodt.github.io/IASDT.R/reference/River_Length.md)
function processes data from the [EU-Hydro River Network
Database](https://land.copernicus.eu/en/products/eu-hydro/eu-hydro-river-network-database)
to compute river lengths categorized by Strahler order for each grid
cell. The [Strahler
order](https://en.wikipedia.org/wiki/Strahler_number), a hierarchical
classification of river networks, assigns higher numbers to larger, more
significant river segments. For each grid cell, the function calculates
the cumulative length of rivers at or above a given Strahler order
(e.g., for Strahler 5, it includes rivers with Strahler values of 5 or
greater). The total river length per grid cell for Strahler order 5 and
above, transformed to a logarithmic scale (log₁₀), serves as a potential
predictor variable in the species distribution models.

  
  

## Sampling efforts

To address the opportunistic bias inherent in the presence-only data,
the total number of vascular plant observations per grid cell from the
Global Biodiversity Information Facility ([GBIF](https://www.gbif.org/))
is employed as a proxy for sampling effort. The
[`efforts_process()`](https://biodt.github.io/IASDT.R/reference/Efforts_data.md)
function manages the request, retrieval, and processing of this sampling
effort data, encompassing over 275 million occurrences as of September
2025. Beyond total observations, the function also determines the number
of vascular plant species per grid cell. The processed data support two
primary applications:

- **grid cell filtering**: The number of species per grid cell enables
  optional filtering to exclude areas with insufficient sampling effort
  (e.g., grid cells with fewer than 100 observed vascular plant species
  are excluded).
- **sampling bias correction**: The total number of observations per
  grid cell, on a logarithmic scale (log₁₀), is incorporated as a
  predictor variable in the models to account for sampling bias. To
  mitigate this bias during predictions, this predictor is fixed at a
  constant value representing optimal sampling effort across the study
  area, as detailed in Warton *et al.*
  ([2013](https://doi.org/10.1371/journal.pone.0079168)).

------------------------------------------------------------------------

**Previous articles:**  
↠[1.
Overview](https://biodt.github.io/IASDT.R/articles/workflow_1_overview.md)  
**Next articles:**  
↠[3. Processing biotic
data](https://biodt.github.io/IASDT.R/articles/workflow_3_biotic_data.md)  
↠[4. Model
fitting](https://biodt.github.io/IASDT.R/articles/workflow_4_model_fitting.md)  
↠[5. Model
post-processing](https://biodt.github.io/IASDT.R/articles/workflow_5_model_postprocess.md)  
