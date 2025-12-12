# IASDT modelling workflow — 4. Model fitting

  
  

This article outlines the preparation of input data for model fitting
and the subsequent process of fitting these models on GPUs within the
`IASDT` workflow.

  
  

## Model input data

The primary function for preparing model-fitting data and initialising
models is
[`mod_prepare_hpc()`](https://biodt.github.io/IASDT.R/reference/Mod_inputs.md).
It structures data for each habitat-specific model into distinct
directories (e.g., datasets/processed/model_fitting/HabX, where *X*
represents a [habitat
type](https://biodt.github.io/IASDT.R/articles/workflow_1_overview.html#models)).
This function orchestrates a suite of specialised sub-functions to
perform the following tasks:

  

- [`mod_prepare_data()`](https://biodt.github.io/IASDT.R/reference/Mod_inputs.md):
  prepare input data for modelling, with key arguments including:

|                       |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
|:----------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| hab_abb               | abbreviation of a single habitat type to be modelled                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| directory_name        | directory path for storing all model files                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| min_efforts_n_species | minimum number of vascular plant species per grid cell required for inclusion in model fitting. This reflects the total count of vascular plant species (including native species) recorded in GBIF across Europe, as computed during the [sampling effort preparation](https://biodt.github.io/IASDT.R/articles/workflow_2_abiotic_data.html#sampling-efforts) step ([`efforts_process()`](https://biodt.github.io/IASDT.R/reference/Efforts_data.md)). This argument filters out grid cells with insufficient sampling effort |
| exclude_cultivated    | whether to exclude countries with cultivated or casual observations for each species                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| exclude_0_habitat     | whether to exclude grid cells with zero [habitat coverage](https://biodt.github.io/IASDT.R/articles/workflow_2_abiotic_data.html#corine-land-cover-habitat-data) of the respective habitat type                                                                                                                                                                                                                                                                                                                                 |
| n_pres_per_species    | minimum number of presence grid cells required for a species to be included in the models, calculated after excluding grid cells with low sampling effort (min_efforts_n_species), zero habitat coverage (exclude_0_habitat), and countries with cultivated or casual observations (exclude_cultivated)                                                                                                                                                                                                                         |

  
  

- [`mod_cv_prepare()`](https://biodt.github.io/IASDT.R/reference/mod_CV_prepare.md):
  prepare and visualise options for spatial-block cross-validation. In
  the *cv_dist* strategy, block size is governed by the *cv_n_grids*
  argument, whereas in the *cv_large* strategy, the study area is
  partitioned into larger blocks based on the *cv_n_rows* and
  *cv_n_columns* arguments.

|                          |                                                                                                                                                                                                                                   |
|:-------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| cv_n_folds               | number of cross-validation folds                                                                                                                                                                                                  |
| cv_n_grids               | number of grid cells in each directions for the *cv_dist* cross-validation strategy (default: 20, yielding 20 × 20 grid cell blocks).                                                                                             |
| cv_n_rows / cv_n_columns | number of rows and columns defining in the *cv_large* cross-validation strategy, partitioning the study area into large blocks (default: *cv_n_rows = cv_n_columns = 2*, resulting in four blocks divided at median coordinates). |

  
  

- [`prepare_knots()`](https://biodt.github.io/IASDT.R/reference/prepare_knots.md):
  prepare and visualise knot locations for Gaussian Predictive Process
  (GPP) models, as described by Tikhonov *et al.*
  ([2019](https://doi.org/10.1002/ecy.2929)).

|                 |                                                                                                                                                                                                                                                                         |
|:----------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| gpp             | whether to incorporate spatial random effects using the Gaussian Predictive Process (GPP)                                                                                                                                                                               |
| gpp_dists       | distance (in kilometres; controlled by the `min_distance` argument of [`prepare_knots()`](https://biodt.github.io/IASDT.R/reference/prepare_knots.md)) specifying both the spacing between knots and the minimum distance between a knot and the nearest sampling point |
| gpp_plot        | whether to plot the coordinates of sampling units and knots                                                                                                                                                                                                             |
| min_lf / max_lf | minimum and maximum number of latent factors to be include                                                                                                                                                                                                              |
| alphapw         | prior specification for the alpha parameter                                                                                                                                                                                                                             |

  
  

- [`mod_slurm()`](https://biodt.github.io/IASDT.R/reference/Mod_SLURM.md):
  generate SLURM scripts to facilitate model fitting on GPUs using the
  `Hmsc-HPC` extension.

|                               |                                            |
|:------------------------------|:-------------------------------------------|
| job_name                      | name assigned to the SLURM job             |
| ntasks                        | number of tasks to execute                 |
| cpus_per_task / gpus_per_node | Number of CPUs and GPUs allocated per node |
| memory_per_cpu                | memory allocation per CPU                  |
| job_runtime                   | maximum duration for job execution         |
| hpc_partition                 | name of the HPC partition                  |
| n_array_jobs                  | number of jobs within each SLURM script    |

  
  

Other arguments:

- selection of predictors:

|                        |                                                                                                                                                                                                                                    |
|:-----------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| bio_variables          | names of [CHELSA](https://biodt.github.io/IASDT.R/articles/workflow_2_abiotic_data.html#chelsa-climate-data) variables to include in the model                                                                                     |
| quadratic_variables    | names of variables for which quadratic terms are incorporated                                                                                                                                                                      |
| efforts_as_predictor   | whether to include the (log₁₀-transformed) [sampling effort](https://biodt.github.io/IASDT.R/articles/workflow_2_abiotic_data.html#sampling-efforts) as a predictor                                                                |
| road_rail_as_predictor | whether to include the (log₁₀-transformed) [summed road and railway intensity](https://biodt.github.io/IASDT.R/articles/workflow_2_abiotic_data.html#railways-and-roads-intensity) as a predictor                                  |
| habitat_as_predictor   | whether to include the (log₁₀-transformed) [percentage coverage of the respective habitat type](https://biodt.github.io/IASDT.R/articles/workflow_2_abiotic_data.html#corine-land-cover-habitat-data) per grid cell as a predictor |
| river_as_predictor     | whether to include the (log₁₀-transformed) total [river length](https://biodt.github.io/IASDT.R/articles/workflow_2_abiotic_data.html#river-length) per grid cell as a predictor                                                   |

  
  

- model fitting options

|                       |                                                                                                                   |
|:----------------------|:------------------------------------------------------------------------------------------------------------------|
| mcmc_n_chains         | number of MCMC chains                                                                                             |
| mcmc_thin             | thinning value(s) in MCMC sampling                                                                                |
| mcmc_samples          | number of MCMC samples per chain                                                                                  |
| mcmc_transient_factor | transient multiplication factor. The value of transient will equal the multiplication of transientFactor and thin |
| mcmc_verbose          | interval at which MCMC sampling progress is reported                                                              |
| precision             | floating-point precision mode for Hmsc-HPC sampling                                                               |

  
  

- n_species_per_grid: minimum number of IAS per grid cell for a grid
  cell to be included in the analysis
- model_country: fit the model for a specific country or countries
- whether or not to use phylogenetic trees: use_phylo_tree and
  no_phylo_tree
- path_hmsc: directory path to `Hmsc-HPC` extension installation
- slurm_prepare: whether to prepare SLURM script for model fitting on
  GPU via
  [`mod_slurm()`](https://biodt.github.io/IASDT.R/reference/Mod_SLURM.md)

------------------------------------------------------------------------

## Model fitting on GPUs

Following the preparation of model input data and initialisation of
models, the subsequent phase involves fitting these models on GPUs. For
each habitat type, the
[`mod_prepare_hpc()`](https://biodt.github.io/IASDT.R/reference/Mod_inputs.md)
function produces:

- python commands (commands_to_fit.txt) for fitting model chains across
  all model variants on GPUs, with each line corresponding to a single
  chain.

»» Example model fitting commands

``` bash
python3 -m hmsc.run_gibbs_sampler --input datasets/processed/model_fitting/Mod_Riv_Hab1/InitMod4HPC/InitMod_GPP120_Tree_samp1000_th200.rds --output datasets/processed/model_fitting/Mod_Riv_Hab1/model_fitting_hpc/GPP120_Tree_samp1000_th200_Chain1_post.rds --samples 1000 --transient 100000 --thin 200 --verbose 200 --chain 0 --fp 64 >& datasets/processed/model_fitting/Mod_Riv_Hab1/model_fitting_hpc/GPP120_Tree_samp1000_th200_Chain1_Progress.txt
python3 -m hmsc.run_gibbs_sampler --input datasets/processed/model_fitting/Mod_Riv_Hab1/InitMod4HPC/InitMod_GPP120_Tree_samp1000_th200.rds --output datasets/processed/model_fitting/Mod_Riv_Hab1/model_fitting_hpc/GPP120_Tree_samp1000_th200_Chain2_post.rds --samples 1000 --transient 100000 --thin 200 --verbose 200 --chain 1 --fp 64 >& datasets/processed/model_fitting/Mod_Riv_Hab1/model_fitting_hpc/GPP120_Tree_samp1000_th200_Chain2_Progress.txt
```

  

- one or more SLURM script files (bash_fit.slurm) designed to submit all
  model-fitting commands (commands_to_fit.txt) as batch jobs on a
  high-performance computing (HPC) system.

  

»» Example SLURM script

``` bash
#!/bin/bash

# -----------------------------------------------------------
# Job array configuration
# -----------------------------------------------------------
#SBATCH --job-name=Hab1
#SBATCH --ntasks=1
#SBATCH --output=datasets/processed/model_fitting/Mod_Riv_Hab1/model_fitting_hpc/JobsLog/%x-%A-%a.out
#SBATCH --error=datasets/processed/model_fitting/Mod_Riv_Hab1/model_fitting_hpc/JobsLog/%x-%A-%a.out
#SBATCH --account=project_465001588
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=1
#SBATCH --time=3-00:00:00
#SBATCH --partition=small-g
#SBATCH --array=1-30

# -----------------------------------------------------------
# Job info
# -----------------------------------------------------------
echo "Start time = $(date)"
echo "Submitting directory = "$SLURM_SUBMIT_DIR
echo "working directory = "$PWD
echo "Project name = "$SLURM_JOB_ACCOUNT
echo "Job id = "$SLURM_JOB_ID
echo "Job name = "$SLURM_JOB_NAME
echo "memory per CPU = "$SLURM_MEM_PER_CPU
echo "The GPU IDs of GPUs in the job allocation (if any) = "$SLURM_JOB_GPUS
echo "Node running the job script = "$SLURMD_NODENAME
echo "Process ID of the process started for the task = "$SLURM_TASK_PID
echo "Dependency = "$SLURM_JOB_DEPENDENCY
echo "Number of nodes assigned to a job = "$SLURM_NNODES
echo "Number of tasks requested by the job = "$SLURM_NTASKS
echo "Number of cpus per task = "$SLURM_CPUS_PER_TASK
echo "Number of tasks in the array = "$SLURM_ARRAY_TASK_COUNT
echo "Array's maximum ID (index) number = "$SLURM_ARRAY_TASK_MIN
echo "Array's minimum ID (index) number = "$SLURM_ARRAY_TASK_MAX

# -----------------------------------------------------------
# File contains bash commands for model fitting
# -----------------------------------------------------------
File=datasets/processed/model_fitting/Mod_Riv_Hab1/Commands2Fit.txt

# -----------------------------------------------------------
# Loading Hmsc-HPC
# -----------------------------------------------------------
source /pfs/lustrep4/scratch/project_465001857/elgabbas/Hmsc_simplify_io/setup-env.sh

# -----------------------------------------------------------
# Check GPU
# -----------------------------------------------------------
export TF_CPP_MIN_LOG_LEVEL=3
PythonCheckGPU=references/LUMI_Check_GPU.py

# -----------------------------------------------------------
# Some checking
# -----------------------------------------------------------
path_python=$(which python3)
echo -e "Some Checking:\n  >>  Working directory: $PWD\n  >>  Python path:       $path_python\n  >>  Checking GPU:      $(python3 $PythonCheckGPU)\n"

# -----------------------------------------------------------
# Run array job
# -----------------------------------------------------------
head -n $SLURM_ARRAY_TASK_ID $File | tail -n 1 | bash

echo "End of program at `date`"


# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# This script was created on: 2025-02-13 19:43 CET
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
```

  
  

Batch jobs for model fitting can be submitted using the `sbatch`
command, for example:

``` bash
sbatch datasets/processed/model_fitting/Hab1/bash_fit.slurm
```

------------------------------------------------------------------------

**Previous articles:**  
↠[1.
Overview](https://biodt.github.io/IASDT.R/articles/workflow_1_overview.md)  
↠[2. Processing abiotic
data](https://biodt.github.io/IASDT.R/articles/workflow_2_abiotic_data.md)  
↠[3. Processing biotic
data](https://biodt.github.io/IASDT.R/articles/workflow_3_biotic_data.md)  
**Next articles:**  
↠[5. Model
post-processing](https://biodt.github.io/IASDT.R/articles/workflow_5_model_postprocess.md)  
