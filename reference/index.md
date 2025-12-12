# Package index

## Package info

- [`IASDT.R`](https://biodt.github.io/IASDT.R/reference/IASDT.R-package.md)
  [`IASDT.R-package`](https://biodt.github.io/IASDT.R/reference/IASDT.R-package.md)
  : IASDT.R: Modelling the Distribution of Invasive Alien Plant Species
  in Europe

## Prepare species distribution data

- [`gbif_process()`](https://biodt.github.io/IASDT.R/reference/GBIF_data.md)
  [`gbif_download()`](https://biodt.github.io/IASDT.R/reference/GBIF_data.md)
  [`gbif_read_chunk()`](https://biodt.github.io/IASDT.R/reference/GBIF_data.md)
  [`gbif_species_data()`](https://biodt.github.io/IASDT.R/reference/GBIF_data.md)
  :

  Process GBIF occurrence data for the `IASDT`

- [`easin_process()`](https://biodt.github.io/IASDT.R/reference/EASIN_data.md)
  [`easin_taxonomy()`](https://biodt.github.io/IASDT.R/reference/EASIN_data.md)
  [`easin_download()`](https://biodt.github.io/IASDT.R/reference/EASIN_data.md)
  [`easin_plot()`](https://biodt.github.io/IASDT.R/reference/EASIN_data.md)
  :

  Process EASIN data for the `IASDT`

- [`elter_process()`](https://biodt.github.io/IASDT.R/reference/eLTER_Process.md)
  :

  Process eLTER data for the `IASDT`

- [`naps_process()`](https://biodt.github.io/IASDT.R/reference/naps_data.md)
  [`naps_distribution()`](https://biodt.github.io/IASDT.R/reference/naps_data.md)
  [`naps_plot()`](https://biodt.github.io/IASDT.R/reference/naps_data.md)
  [`naps_standardisation()`](https://biodt.github.io/IASDT.R/reference/naps_data.md)
  :

  Process and map Naturalized Alien Plant Species (NAPS) data for the
  `IASDT`

- [`get_species_name()`](https://biodt.github.io/IASDT.R/reference/get_species_name.md)
  :

  Get species name or information of an `IASDT` species ID

## Prepare abiotic data

- [`clc_process()`](https://biodt.github.io/IASDT.R/reference/CLC_Process.md)
  :

  Process Corine Land Cover (CLC) data for the `IASDT`

- [`chelsa_variables`](https://biodt.github.io/IASDT.R/reference/CHELSA_variables.md)
  : Detailed information on CHELSA climate variables

- [`chelsa_process()`](https://biodt.github.io/IASDT.R/reference/CHELSA_data.md)
  [`chelsa_prepare()`](https://biodt.github.io/IASDT.R/reference/CHELSA_data.md)
  [`chelsa_project()`](https://biodt.github.io/IASDT.R/reference/CHELSA_data.md)
  :

  Process CHELSA Climate Data for the `IASDT`

- [`efforts_process()`](https://biodt.github.io/IASDT.R/reference/Efforts_data.md)
  [`efforts_request()`](https://biodt.github.io/IASDT.R/reference/Efforts_data.md)
  [`efforts_download()`](https://biodt.github.io/IASDT.R/reference/Efforts_data.md)
  [`efforts_summarize()`](https://biodt.github.io/IASDT.R/reference/Efforts_data.md)
  [`efforts_split()`](https://biodt.github.io/IASDT.R/reference/Efforts_data.md)
  [`efforts_plot()`](https://biodt.github.io/IASDT.R/reference/Efforts_data.md)
  :

  Process GBIF sampling effort data for the `IASDT`

- [`railway_intensity()`](https://biodt.github.io/IASDT.R/reference/Railway_Intensity.md)
  :

  Calculate railway intensity based on `OpenStreetMap` data

- [`road_intensity()`](https://biodt.github.io/IASDT.R/reference/Road_Intensity.md)
  : Calculate road intensity per grid cell

- [`river_length()`](https://biodt.github.io/IASDT.R/reference/River_Length.md)
  : Calculate the length of rivers in each Strahler order per grid cell

- [`bioreg_process()`](https://biodt.github.io/IASDT.R/reference/BioReg_Process.md)
  : Process biogeographical regions dataset

- [`wetness_index_process()`](https://biodt.github.io/IASDT.R/reference/wetness_index_process.md)
  : Download and Process Topographic Wetness Index Data

- [`soil_density_process()`](https://biodt.github.io/IASDT.R/reference/soil_density_process.md)
  : Retrieve and project soil bulk density data

## Modelling functions

Functions for preparing data, running the models, and postprocessing of
model outputs

### Data preparation

Prepare input data and scripts for fitting Hmsc-HPC on GPU

- [`mod_cv_prepare()`](https://biodt.github.io/IASDT.R/reference/mod_CV_prepare.md)
  : Prepare spatial-block cross-validation folds for spatial analysis
- [`prepare_knots()`](https://biodt.github.io/IASDT.R/reference/prepare_knots.md)
  : Prepare knot locations for Hmsc GPP models
- [`mod_prepare_data()`](https://biodt.github.io/IASDT.R/reference/Mod_inputs.md)
  [`mod_prepare_hpc()`](https://biodt.github.io/IASDT.R/reference/Mod_inputs.md)
  : Prepare initial models for model fitting with Hmsc-HPC
- [`mod_slurm()`](https://biodt.github.io/IASDT.R/reference/Mod_SLURM.md)
  [`mod_slurm_refit()`](https://biodt.github.io/IASDT.R/reference/Mod_SLURM.md)
  : Prepare SLURM scripts for Hmsc-HPC model fitting
- [`mod_fit_windows()`](https://biodt.github.io/IASDT.R/reference/mod_fit_windows.md)
  : Fit Hmsc-HPC models on UFZ Windows Server
- [`mod_cv_fit()`](https://biodt.github.io/IASDT.R/reference/Mod_CV_Fit.md)
  : Prepare cross-validated Hmsc models for HPC fitting
- [`install_hmsc_windows()`](https://biodt.github.io/IASDT.R/reference/install_hmsc_windows.md)
  : Install Hmsc-HPC in a python virtual environment on Windows
- [`solve1()`](https://biodt.github.io/IASDT.R/reference/cpp_functions.md)
  [`solve2()`](https://biodt.github.io/IASDT.R/reference/cpp_functions.md)
  [`solve2vect()`](https://biodt.github.io/IASDT.R/reference/cpp_functions.md)
  [`fast_pnorm()`](https://biodt.github.io/IASDT.R/reference/cpp_functions.md)
  [`exp_neg_div()`](https://biodt.github.io/IASDT.R/reference/cpp_functions.md)
  : helper C++ functions for fast matrix computations
- [`fit_sdm_models()`](https://biodt.github.io/IASDT.R/reference/fit_sdm_models.md)
  : Species Distribution Modelling Workflow for Single-Species Models

### Model postprocessing

Postprocessing model outputs, including checking for convergence, making
spatial predictions, evaluation, and plotting.

- [`coda_to_tibble()`](https://biodt.github.io/IASDT.R/reference/Coda_to_tibble.md)
  : Convert a Coda object to a tibble with specified parameter
  transformations

- [`trim_hmsc()`](https://biodt.github.io/IASDT.R/reference/trim_hmsc.md)
  : Trim an Hmsc Model Object by Removing Specified Components

- [`mod_get_posteriors()`](https://biodt.github.io/IASDT.R/reference/mod_get_posteriors.md)
  :

  Combines posteriors exported by `Hmsc-HPC` into an Hmsc object

- [`mod_merge_chains()`](https://biodt.github.io/IASDT.R/reference/Mod_Merge_Chains.md)
  [`mod_merge_chains_cv()`](https://biodt.github.io/IASDT.R/reference/Mod_Merge_Chains.md)
  :

  Merge model chains into `Hmsc` and `coda` objects

- [`mod_summary()`](https://biodt.github.io/IASDT.R/reference/Mod_Summary.md)
  : Summary of Hmsc model parameters

- [`rc_prepare_data()`](https://biodt.github.io/IASDT.R/reference/Response_curves.md)
  [`rc_plot_species()`](https://biodt.github.io/IASDT.R/reference/Response_curves.md)
  [`rc_plot_species_all()`](https://biodt.github.io/IASDT.R/reference/Response_curves.md)
  [`rc_plot_sr()`](https://biodt.github.io/IASDT.R/reference/Response_curves.md)
  : Prepare and plot response curve data for Hmsc models

- [`predict_latent_factor()`](https://biodt.github.io/IASDT.R/reference/predict_latent_factor.md)
  : Draws samples from the conditional predictive distribution of latent
  factors

- [`plot_latent_factor()`](https://biodt.github.io/IASDT.R/reference/plot_latent_factor.md)
  : Plot spatial variation in site loadings of HMSC models

- [`predict_hmsc()`](https://biodt.github.io/IASDT.R/reference/Predict_Hmsc.md)
  : Calculates predicted values from a fitted Hmsc model

- [`predict_maps()`](https://biodt.github.io/IASDT.R/reference/Predict_Maps.md)
  [`predict_maps_cv()`](https://biodt.github.io/IASDT.R/reference/Predict_Maps.md)
  :

  Predict habitat suitability of `Hmsc` models

- [`plot_prediction()`](https://biodt.github.io/IASDT.R/reference/plot_prediction.md)
  :

  Plot species and level of invasion predictions as JPEG files using
  `ggplot2`

- [`mod_cv_evaluate()`](https://biodt.github.io/IASDT.R/reference/mod_cv_evaluate.md)
  : Cross-validation Model Evaluation and Plotting

- [`mod_cv_evaluate_plot()`](https://biodt.github.io/IASDT.R/reference/mod_cv_evaluate_plot.md)
  : Plot Evaluation Results for Cross-Validated Hmsc Models

- [`mod_cv_merge_predictions()`](https://biodt.github.io/IASDT.R/reference/mod_cv_merge_predictions.md)
  : Merge Cross-Validation Predictions and Generate Summary Maps

- [`mod_postprocess_1_cpu()`](https://biodt.github.io/IASDT.R/reference/Mod_postprocessing.md)
  [`mod_prepare_tf()`](https://biodt.github.io/IASDT.R/reference/Mod_postprocessing.md)
  [`mod_postprocess_2_cpu()`](https://biodt.github.io/IASDT.R/reference/Mod_postprocessing.md)
  [`mod_postprocess_cv_1_cpu()`](https://biodt.github.io/IASDT.R/reference/Mod_postprocessing.md)
  [`mod_postprocess_cv_2_cpu()`](https://biodt.github.io/IASDT.R/reference/Mod_postprocessing.md)
  : Model pipeline for post-processing fitted Hmsc models

- [`plot_evaluation()`](https://biodt.github.io/IASDT.R/reference/plot_evaluation.md)
  : Generate plots for the explanatory power of Hmsc models

- [`mod_heatmap_beta()`](https://biodt.github.io/IASDT.R/reference/Parameter_Heatmap.md)
  [`mod_heatmap_omega()`](https://biodt.github.io/IASDT.R/reference/Parameter_Heatmap.md)
  :

  Heatmaps for the `beta` and `omega` parameters of the Hmsc model

- [`convergence_plot()`](https://biodt.github.io/IASDT.R/reference/Convergence_plots.md)
  [`convergence_alpha()`](https://biodt.github.io/IASDT.R/reference/Convergence_plots.md)
  [`convergence_rho()`](https://biodt.github.io/IASDT.R/reference/Convergence_plots.md)
  [`convergence_beta_ranges()`](https://biodt.github.io/IASDT.R/reference/Convergence_plots.md)
  : Plot model convergence of a selected model

- [`convergence_plot_all()`](https://biodt.github.io/IASDT.R/reference/Convergence_Plot_All.md)
  : Plot model convergence of multiple modelling alternatives

- [`variance_partitioning_compute()`](https://biodt.github.io/IASDT.R/reference/Variance_partitioning.md)
  [`variance_partitioning_plot()`](https://biodt.github.io/IASDT.R/reference/Variance_partitioning.md)
  : Computes and visualise variance partitioning of Hmsc models

- [`plot_gelman()`](https://biodt.github.io/IASDT.R/reference/plot_gelman.md)
  [`plot_gelman_alpha()`](https://biodt.github.io/IASDT.R/reference/plot_gelman.md)
  [`plot_gelman_beta()`](https://biodt.github.io/IASDT.R/reference/plot_gelman.md)
  [`plot_gelman_omega()`](https://biodt.github.io/IASDT.R/reference/plot_gelman.md)
  [`plot_gelman_rho()`](https://biodt.github.io/IASDT.R/reference/plot_gelman.md)
  : Plot Gelman-Rubin-Brooks
