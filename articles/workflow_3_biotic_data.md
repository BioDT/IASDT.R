# IASDT modelling workflow — 3. biotic data

  
  

This article details the processing of biotic data within the `IASDT`
modelling workflow. The biotic data are sourced from `GBIF`, `EASIN`,
and `eLTER`; see below. All data preparation adheres to the
[FAIR](https://doi.org/10.1038/sdata.2016.18) principles (Findable,
Accessible, Interoperable, and Reusable) to ensure scientific integrity
and reproducibility.

  
  

## IAS species list

The most recent checklist of naturalised terrestrial alien plant species
of non-European origin was retrieved from
[FloraVeg.EU](https://floraveg.eu/). This checklist was standardized
against the GBIF taxonomic backbone to ensure taxonomic consistency,
resulting in a total of 1,319 species (September, 2025). The table below
summarises the species count by taxonomic class and order. Occurrence
data for these species were compiled from available sources detailed in
subsequent sections.

| class          | order          | \# species | class          | order        | \# species |
|:---------------|:---------------|:-----------|:---------------|:-------------|:-----------|
| Liliopsida     | Acorales       | 2          | Magnoliopsida  | Gentianales  | 9          |
| Liliopsida     | Alismatales    | 8          | Magnoliopsida  | Geraniales   | 9          |
| Liliopsida     | Arecales       | 5          | Magnoliopsida  | Gunnerales   | 1          |
| Liliopsida     | Asparagales    | 74         | Magnoliopsida  | Lamiales     | 63         |
| Liliopsida     | Commelinales   | 4          | Magnoliopsida  | Laurales     | 1          |
| Liliopsida     | Liliales       | 6          | Magnoliopsida  | Magnoliales  | 1          |
| Liliopsida     | Poales         | 163        | Magnoliopsida  | Malpighiales | 41         |
| Liliopsida     | Zingiberales   | 7          | Magnoliopsida  | Malvales     | 10         |
| Lycopodiopsida | Selaginellales | 1          | Magnoliopsida  | Myrtales     | 74         |
| Magnoliopsida  | Apiales        | 15         | Magnoliopsida  | Oxalidales   | 19         |
| Magnoliopsida  | Asterales      | 148        | Magnoliopsida  | Piperales    | 3          |
| Magnoliopsida  | Boraginales    | 10         | Magnoliopsida  | Proteales    | 4          |
| Magnoliopsida  | Brassicales    | 13         | Magnoliopsida  | Ranunculales | 27         |
| Magnoliopsida  | Caryophyllales | 141        | Magnoliopsida  | Rosales      | 125        |
| Magnoliopsida  | Celastrales    | 4          | Magnoliopsida  | Sapindales   | 23         |
| Magnoliopsida  | Cornales       | 10         | Magnoliopsida  | Saxifragales | 38         |
| Magnoliopsida  | Cucurbitales   | 6          | Magnoliopsida  | Solanales    | 51         |
| Magnoliopsida  | Dipsacales     | 19         | Magnoliopsida  | Vitales      | 9          |
| Magnoliopsida  | Ericales       | 40         | Pinopsida      | Pinales      | 46         |
| Magnoliopsida  | Escalloniales  | 1          | Polypodiopsida | Cyatheales   | 1          |
| Magnoliopsida  | Fabales        | 52         | Polypodiopsida | Equisetales  | 1          |
| Magnoliopsida  | Fagales        | 20         | Polypodiopsida | Polypodiales | 12         |
| Magnoliopsida  | Garryales      | 1          | Polypodiopsida | Schizaeales  | 1          |

  
  

## GBIF

The
[`gbif_process()`](https://biodt.github.io/IASDT.R/reference/GBIF_data.md)
function manages the retrieval, processing, and visualisation of the
most current occurrence data for invasive alien species (IAS) from the
Global Biodiversity Information Facility
([GBIF](https://www.gbif.org/)), totalling approximately 10.4 million
occurrences as of September 2025. Occurrences deemed doubtful or
exhibiting high spatial uncertainty are excluded to ensure data quality.
Additionally, the function calculates the number of observations and
occupied grid cells per species. Its primary output consists of gridded
distribution maps for each species, delineating presence-only data
across the reference grid.

  
  

## EASIN

The European Alien Species Information Network
([EASIN](https://easin.jrc.ec.europa.eu/)) serves as a web-based
platform providing spatial data on approximately 14,000 alien species.
The
[`easin_process()`](https://biodt.github.io/IASDT.R/reference/EASIN_data.md)
function oversees the processing and visualisation of the latest IAS
occurrence data sourced from EASIN. A checklist of vascular plants from
EASIN’s database was retrieved via its API and standardized against the
GBIF taxonomic backbone. Occurrence data for taxa corresponding to the
IAS checklist were subsequently downloaded through the API. Of the 34
partner organizations contributing to EASIN (including GBIF), only
non-GBIF data are incorporated into the models to avoid redundancy.
Following processing and spatial filtering, the dataset includes over
427K observations for 465 IAS as of September 2025. Consistent with the
GBIF output, the function generates gridded distribution maps for each
species.

  
  

## eLTER

The Integrated European Long-Term Ecosystem, Critical Zone, and
Socio-Ecological Research Infrastructure ([eLTER](https://elter-ri.eu/))
comprises a network of sites dedicated to collecting ecological data for
long-term research across the European Union. Cleaned standardised data
encompassed 4K observations for 109 IAS (September 2025). The
[`elter_process()`](https://biodt.github.io/IASDT.R/reference/eLTER_Process.md)
function processed standardised eLTER data in a format ready to be used
in the models.

  
  

## Species distribution

The
[`naps_process()`](https://biodt.github.io/IASDT.R/reference/naps_data.md)
function integrates and visualises species distribution data sourced
from `GBIF`, `EASIN`, and `eLTER`. It combines pre-processed datasets
and generates presence-absence rasters for use in the species
distribution models.

------------------------------------------------------------------------

**Previous articles:**  
↠[1.
Overview](https://biodt.github.io/IASDT.R/articles/workflow_1_overview.md)  
↠[2. Processing abiotic
data](https://biodt.github.io/IASDT.R/articles/workflow_2_abiotic_data.md)  
**Next articles:**  
↠[4. Model
fitting](https://biodt.github.io/IASDT.R/articles/workflow_4_model_fitting.md)  
↠[5. Model
post-processing](https://biodt.github.io/IASDT.R/articles/workflow_5_model_postprocess.md)  
