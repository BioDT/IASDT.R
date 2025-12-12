# IASDT modelling workflow — 5. Model post-processing

  

Post-processing of fitted models in the `IASDT` workflow involves
multiple steps, utilising both CPU and GPU computations to optimize
performance and manage memory constraints effectively.

  

## Step 1: CPU

The
**[`mod_postprocess_1_cpu()`](https://biodt.github.io/IASDT.R/reference/Mod_postprocessing.md)**
function begins the post-processing phase for each habitat type by
automating the following tasks:

[TABLE]

  

------------------------------------------------------------------------

  

### Computationally intensive tasks offloaded to GPU

Previous attempts to prepare response curve data, predict at new sites,
and compute variance partitioning using R on CPUs (such as the UFZ
Windows server and LUMI HPC) were limited by memory constraints. As a
result, these tasks are now offloaded to GPU-based computations using
`Python` and `TensorFlow`. The
[`mod_postprocess_1_cpu()`](https://biodt.github.io/IASDT.R/reference/Mod_postprocessing.md)
function calls the following sub-functions to generate the necessary
commands for GPU execution:

|                                                                                                         |                                                                   |
|:--------------------------------------------------------------------------------------------------------|:------------------------------------------------------------------|
| [`rc_prepare_data()`](https://biodt.github.io/IASDT.R/reference/Response_curves.md)                     | prepares data for predicting latent factors for response curves   |
| [`predict_maps()`](https://biodt.github.io/IASDT.R/reference/Predict_Maps.md)                           | prepares data for predicting latent factors at new sampling units |
| [`variance_partitioning_compute()`](https://biodt.github.io/IASDT.R/reference/Variance_partitioning.md) | prepares data for computing variance partitioning                 |

  

------------------------------------------------------------------------

  

### Preparing commands for GPU computations

> **Predicting latent factors:**

- Predictions of latent factors for response curves and new sampling
  units are performed using a `TensorFlow` script located at
  [inst/crossprod_solve.py](https://github.com/BioDT/IASDT.R/blob/main/inst/crossprod_solve.py).
- For these tasks, the corresponding R functions export multiple `.qs2`
  and `.feather` data files to the temp_pred subdirectory, which are
  essential for GPU computations. Additionally, they generate execution
  commands saved as lf_rc_commands\_*.txt* (for response curves) and
  lf_new_sites_commands\_.txt (for new sites).

»» Example lf_rc_commands.txt file

``` bash
Error: File not found.
```

  

»» Example lf_new_sites_Commands.txt file

``` bash
Error: File not found.
```

  
  

- Response curves
  - [`rc_prepare_data()`](https://biodt.github.io/IASDT.R/reference/Response_curves.md)
    extends the functionality of
    [`Hmsc::constructGradient()`](https://rdrr.io/pkg/Hmsc/man/constructGradient.html)
    and
    [`Hmsc::plotGradient()`](https://rdrr.io/pkg/Hmsc/man/plotGradient.html)
    by enabling the preparation of response curve data on GPUs when
    `lf_commands_only = TRUE`.
  - For predictions at mean coordinates (as specified by the
    `coordinates` argument in
    [`Hmsc::constructGradient()`](https://rdrr.io/pkg/Hmsc/man/constructGradient.html)),
    latent factor predictions — which are typically memory-intensive
    when using
    [`Hmsc::predictLatentFactor()`](https://rdrr.io/pkg/Hmsc/man/predictLatentFactor.html)
    — are computed on GPUs.
- Predicting at new sites
  - [`predict_maps()`](https://biodt.github.io/IASDT.R/reference/Predict_Maps.md)
    sets up GPU computations for predictions at new sites when both
    `lf_only = TRUE` and `lf_commands_only = TRUE`.

  
  

> **Computing variance partitioning:**

- Variance partitioning computations are performed on GPUs using
  `TensorFlow` scripts located at
  [inst/VP_geta.py](https://github.com/BioDT/IASDT.R/blob/main/inst/VP_geta.py),
  [inst/VP_getf.py](https://github.com/BioDT/IASDT.R/blob/main/inst/VP_getf.py),
  and
  [inst/VP_gemu.py](https://github.com/BioDT/IASDT.R/blob/main/inst/VP_gemu.py).
  These scripts implement functionality derived from
  [`Hmsc::computeVariancePartitioning()`](https://rdrr.io/pkg/Hmsc/man/computeVariancePartitioning.html),
  specifically the internal functions `geta`, `getf`, and `gemu`.
- The
  [`variance_partitioning_compute()`](https://biodt.github.io/IASDT.R/reference/Variance_partitioning.md)
  function exports the necessary files to the temp_vp subdirectory,
  including numerous `.qs2` and `.feather` files. It also generates
  execution commands saved as VP_A_Command.txt, VP_F_Command.txt, and
  VP_mu_Command.txt.

  
  

### Combining commands for GPU computations

Once
[`mod_postprocess_1_cpu()`](https://biodt.github.io/IASDT.R/reference/Mod_postprocessing.md)
has been executed for all habitat types, the
**[`mod_prepare_tf()`](https://biodt.github.io/IASDT.R/reference/Mod_postprocessing.md)**
function consolidates the batch scripts for GPU computations across all
habitat types:

- It aggregates the script files that contain commands for response
  curves and latent factor predictions, splitting them into multiple
  scripts (tf_chunk\_\*.txt) for batch processing. Additionally, it
  generates a SLURM script (lf_SLURM.slurm) for executing the latent
  factor predictions.

»» Example tf_chunk\_\*.txt file

``` bash
#!/bin/bash

# Load TensorFlow module and configure environment
ml use /appl/local/csc/modulefiles
ml tensorflow
export TF_CPP_MIN_LOG_LEVEL=3
export TF_ENABLE_ONEDNN_OPTS=0

# Verify GPU availability
python3 -c "import tensorflow as tf; print(\"Num GPUs Available:\", len(tf.config.list_physical_devices(\"GPU\")))"

# 20 commands to be executed:
python3 crossprod_solve.py --s1 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s1.feather' --s2 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s2.feather' --post_eta 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_postEta_ch001.feather' --path_out 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch001.feather' --denom 50000 --chunk_size 1000 --threshold_mb 2000 --solve_chunk_size 50 --verbose  >> datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch001.log 2>&1
python3 crossprod_solve.py --s1 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s1.feather' --s2 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s2.feather' --post_eta 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_postEta_ch002.feather' --path_out 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch002.feather' --denom 1470707 --chunk_size 1000 --threshold_mb 2000 --solve_chunk_size 50 --verbose  >> datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch002.log 2>&1
python3 crossprod_solve.py --s1 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s1.feather' --s2 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s2.feather' --post_eta 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_postEta_ch003.feather' --path_out 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch003.feather' --denom 1485354 --chunk_size 1000 --threshold_mb 2000 --solve_chunk_size 50 --verbose  >> datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch003.log 2>&1
python3 crossprod_solve.py --s1 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s1.feather' --s2 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s2.feather' --post_eta 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_postEta_ch004.feather' --path_out 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch004.feather' --denom 1456061 --chunk_size 1000 --threshold_mb 2000 --solve_chunk_size 50 --verbose  >> datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch004.log 2>&1
python3 crossprod_solve.py --s1 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s1.feather' --s2 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s2.feather' --post_eta 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_postEta_ch005.feather' --path_out 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch005.feather' --denom 1500000 --chunk_size 1000 --threshold_mb 2000 --solve_chunk_size 50 --verbose  >> datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch005.log 2>&1
python3 crossprod_solve.py --s1 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s1.feather' --s2 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s2.feather' --post_eta 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_postEta_ch006.feather' --path_out 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch006.feather' --denom 1500000 --chunk_size 1000 --threshold_mb 2000 --solve_chunk_size 50 --verbose  >> datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch006.log 2>&1
python3 crossprod_solve.py --s1 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s1.feather' --s2 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s2.feather' --post_eta 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_postEta_ch007.feather' --path_out 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch007.feather' --denom 1426768 --chunk_size 1000 --threshold_mb 2000 --solve_chunk_size 50 --verbose  >> datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch007.log 2>&1
python3 crossprod_solve.py --s1 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s1.feather' --s2 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s2.feather' --post_eta 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_postEta_ch008.feather' --path_out 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch008.feather' --denom 1470707 --chunk_size 1000 --threshold_mb 2000 --solve_chunk_size 50 --verbose  >> datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch008.log 2>&1
python3 crossprod_solve.py --s1 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s1.feather' --s2 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s2.feather' --post_eta 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_postEta_ch009.feather' --path_out 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch009.feather' --denom 1485354 --chunk_size 1000 --threshold_mb 2000 --solve_chunk_size 50 --verbose  >> datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch009.log 2>&1
python3 crossprod_solve.py --s1 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s1.feather' --s2 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s2.feather' --post_eta 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_postEta_ch010.feather' --path_out 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch010.feather' --denom 1456061 --chunk_size 1000 --threshold_mb 2000 --solve_chunk_size 50 --verbose  >> datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch010.log 2>&1
python3 crossprod_solve.py --s1 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s1.feather' --s2 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s2.feather' --post_eta 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_postEta_ch011.feather' --path_out 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch011.feather' --denom 1470707 --chunk_size 1000 --threshold_mb 2000 --solve_chunk_size 50 --verbose  >> datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch011.log 2>&1
python3 crossprod_solve.py --s1 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s1.feather' --s2 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s2.feather' --post_eta 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_postEta_ch012.feather' --path_out 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch012.feather' --denom 1412121 --chunk_size 1000 --threshold_mb 2000 --solve_chunk_size 50 --verbose  >> datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch012.log 2>&1
python3 crossprod_solve.py --s1 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s1.feather' --s2 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s2.feather' --post_eta 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_postEta_ch013.feather' --path_out 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch013.feather' --denom 1426768 --chunk_size 1000 --threshold_mb 2000 --solve_chunk_size 50 --verbose  >> datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch013.log 2>&1
python3 crossprod_solve.py --s1 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s1.feather' --s2 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s2.feather' --post_eta 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_postEta_ch014.feather' --path_out 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch014.feather' --denom 1426768 --chunk_size 1000 --threshold_mb 2000 --solve_chunk_size 50 --verbose  >> datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch014.log 2>&1
python3 crossprod_solve.py --s1 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s1.feather' --s2 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s2.feather' --post_eta 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_postEta_ch015.feather' --path_out 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch015.feather' --denom 1441414 --chunk_size 1000 --threshold_mb 2000 --solve_chunk_size 50 --verbose  >> datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch015.log 2>&1
python3 crossprod_solve.py --s1 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s1.feather' --s2 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s2.feather' --post_eta 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_postEta_ch016.feather' --path_out 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch016.feather' --denom 1456061 --chunk_size 1000 --threshold_mb 2000 --solve_chunk_size 50 --verbose  >> datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch016.log 2>&1
python3 crossprod_solve.py --s1 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s1.feather' --s2 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s2.feather' --post_eta 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_postEta_ch017.feather' --path_out 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch017.feather' --denom 1441414 --chunk_size 1000 --threshold_mb 2000 --solve_chunk_size 50 --verbose  >> datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch017.log 2>&1
python3 crossprod_solve.py --s1 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s1.feather' --s2 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s2.feather' --post_eta 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_postEta_ch018.feather' --path_out 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch018.feather' --denom 1382828 --chunk_size 1000 --threshold_mb 2000 --solve_chunk_size 50 --verbose  >> datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch018.log 2>&1
python3 crossprod_solve.py --s1 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s1.feather' --s2 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s2.feather' --post_eta 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_postEta_ch019.feather' --path_out 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch019.feather' --denom 1485354 --chunk_size 1000 --threshold_mb 2000 --solve_chunk_size 50 --verbose  >> datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch019.log 2>&1
python3 crossprod_solve.py --s1 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s1.feather' --s2 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_s2.feather' --post_eta 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_postEta_ch020.feather' --path_out 'datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch020.feather' --denom 1368182 --chunk_size 1000 --threshold_mb 2000 --solve_chunk_size 50 --verbose  >> datasets/processed/model_fitting/Mod_Riv_Hab1/temp_pred/lf_1_Test_etaPred_ch020.log 2>&1
```

  

»» Example lf_SLURM.slurm file

``` bash
Error: File not found.
```

  

- It combines the variance partitioning command files into a single
  VP_Commands.txt file and prepares a SLURM script (VP_SLURM.slurm) for
  the variance partitioning computations.

»» Example VP_Commands.txt file

``` bash
python3 VP_gemu.py --tr datasets/processed/model_fitting/Mod_Q_Hab1/temp_vp/VP_Tr.feather --gamma datasets/processed/model_fitting/Mod_Q_Hab1/temp_vp/VP_Gamma.feather --output datasets/processed/model_fitting/Mod_Q_Hab1/temp_vp/VP_Mu.feather --ncores 3 --chunk_size 50 >> datasets/processed/model_fitting/Mod_Q_Hab1/temp_vp/VP_Mu.log 2>&1
python3 VP_gemu.py --tr datasets/processed/model_fitting/Mod_Q_Hab2/temp_vp/VP_Tr.feather --gamma datasets/processed/model_fitting/Mod_Q_Hab2/temp_vp/VP_Gamma.feather --output datasets/processed/model_fitting/Mod_Q_Hab2/temp_vp/VP_Mu.feather --ncores 3 --chunk_size 50 >> datasets/processed/model_fitting/Mod_Q_Hab2/temp_vp/VP_Mu.log 2>&1
python3 VP_gemu.py --tr datasets/processed/model_fitting/Mod_Q_Hab3/temp_vp/VP_Tr.feather --gamma datasets/processed/model_fitting/Mod_Q_Hab3/temp_vp/VP_Gamma.feather --output datasets/processed/model_fitting/Mod_Q_Hab3/temp_vp/VP_Mu.feather --ncores 3 --chunk_size 50 >> datasets/processed/model_fitting/Mod_Q_Hab3/temp_vp/VP_Mu.log 2>&1
python3 VP_gemu.py --tr datasets/processed/model_fitting/Mod_Q_Hab4a/temp_vp/VP_Tr.feather --gamma datasets/processed/model_fitting/Mod_Q_Hab4a/temp_vp/VP_Gamma.feather --output datasets/processed/model_fitting/Mod_Q_Hab4a/temp_vp/VP_Mu.feather --ncores 3 --chunk_size 50 >> datasets/processed/model_fitting/Mod_Q_Hab4a/temp_vp/VP_Mu.log 2>&1
python3 VP_gemu.py --tr datasets/processed/model_fitting/Mod_Q_Hab4b/temp_vp/VP_Tr.feather --gamma datasets/processed/model_fitting/Mod_Q_Hab4b/temp_vp/VP_Gamma.feather --output datasets/processed/model_fitting/Mod_Q_Hab4b/temp_vp/VP_Mu.feather --ncores 3 --chunk_size 50 >> datasets/processed/model_fitting/Mod_Q_Hab4b/temp_vp/VP_Mu.log 2>&1
python3 VP_gemu.py --tr datasets/processed/model_fitting/Mod_Q_Hab10/temp_vp/VP_Tr.feather --gamma datasets/processed/model_fitting/Mod_Q_Hab10/temp_vp/VP_Gamma.feather --output datasets/processed/model_fitting/Mod_Q_Hab10/temp_vp/VP_Mu.feather --ncores 3 --chunk_size 50 >> datasets/processed/model_fitting/Mod_Q_Hab10/temp_vp/VP_Mu.log 2>&1
python3 VP_gemu.py --tr datasets/processed/model_fitting/Mod_Q_Hab12a/temp_vp/VP_Tr.feather --gamma datasets/processed/model_fitting/Mod_Q_Hab12a/temp_vp/VP_Gamma.feather --output datasets/processed/model_fitting/Mod_Q_Hab12a/temp_vp/VP_Mu.feather --ncores 3 --chunk_size 50 >> datasets/processed/model_fitting/Mod_Q_Hab12a/temp_vp/VP_Mu.log 2>&1
python3 VP_gemu.py --tr datasets/processed/model_fitting/Mod_Q_Hab12b/temp_vp/VP_Tr.feather --gamma datasets/processed/model_fitting/Mod_Q_Hab12b/temp_vp/VP_Gamma.feather --output datasets/processed/model_fitting/Mod_Q_Hab12b/temp_vp/VP_Mu.feather --ncores 3 --chunk_size 50 >> datasets/processed/model_fitting/Mod_Q_Hab12b/temp_vp/VP_Mu.log 2>&1
python3 VP_geta.py --tr datasets/processed/model_fitting/Mod_Q_Hab1/temp_vp/VP_Tr.feather --x datasets/processed/model_fitting/Mod_Q_Hab1/temp_vp/VP_X.feather --gamma datasets/processed/model_fitting/Mod_Q_Hab1/temp_vp/VP_Gamma.feather --output datasets/processed/model_fitting/Mod_Q_Hab1/temp_vp/VP_A.feather --ncores 3 --chunk_size 50 >> datasets/processed/model_fitting/Mod_Q_Hab1/temp_vp/VP_A.log 2>&1
python3 VP_geta.py --tr datasets/processed/model_fitting/Mod_Q_Hab2/temp_vp/VP_Tr.feather --x datasets/processed/model_fitting/Mod_Q_Hab2/temp_vp/VP_X.feather --gamma datasets/processed/model_fitting/Mod_Q_Hab2/temp_vp/VP_Gamma.feather --output datasets/processed/model_fitting/Mod_Q_Hab2/temp_vp/VP_A.feather --ncores 3 --chunk_size 50 >> datasets/processed/model_fitting/Mod_Q_Hab2/temp_vp/VP_A.log 2>&1
python3 VP_geta.py --tr datasets/processed/model_fitting/Mod_Q_Hab3/temp_vp/VP_Tr.feather --x datasets/processed/model_fitting/Mod_Q_Hab3/temp_vp/VP_X.feather --gamma datasets/processed/model_fitting/Mod_Q_Hab3/temp_vp/VP_Gamma.feather --output datasets/processed/model_fitting/Mod_Q_Hab3/temp_vp/VP_A.feather --ncores 3 --chunk_size 50 >> datasets/processed/model_fitting/Mod_Q_Hab3/temp_vp/VP_A.log 2>&1
python3 VP_geta.py --tr datasets/processed/model_fitting/Mod_Q_Hab4a/temp_vp/VP_Tr.feather --x datasets/processed/model_fitting/Mod_Q_Hab4a/temp_vp/VP_X.feather --gamma datasets/processed/model_fitting/Mod_Q_Hab4a/temp_vp/VP_Gamma.feather --output datasets/processed/model_fitting/Mod_Q_Hab4a/temp_vp/VP_A.feather --ncores 3 --chunk_size 50 >> datasets/processed/model_fitting/Mod_Q_Hab4a/temp_vp/VP_A.log 2>&1
python3 VP_geta.py --tr datasets/processed/model_fitting/Mod_Q_Hab4b/temp_vp/VP_Tr.feather --x datasets/processed/model_fitting/Mod_Q_Hab4b/temp_vp/VP_X.feather --gamma datasets/processed/model_fitting/Mod_Q_Hab4b/temp_vp/VP_Gamma.feather --output datasets/processed/model_fitting/Mod_Q_Hab4b/temp_vp/VP_A.feather --ncores 3 --chunk_size 50 >> datasets/processed/model_fitting/Mod_Q_Hab4b/temp_vp/VP_A.log 2>&1
python3 VP_geta.py --tr datasets/processed/model_fitting/Mod_Q_Hab10/temp_vp/VP_Tr.feather --x datasets/processed/model_fitting/Mod_Q_Hab10/temp_vp/VP_X.feather --gamma datasets/processed/model_fitting/Mod_Q_Hab10/temp_vp/VP_Gamma.feather --output datasets/processed/model_fitting/Mod_Q_Hab10/temp_vp/VP_A.feather --ncores 3 --chunk_size 50 >> datasets/processed/model_fitting/Mod_Q_Hab10/temp_vp/VP_A.log 2>&1
python3 VP_geta.py --tr datasets/processed/model_fitting/Mod_Q_Hab12a/temp_vp/VP_Tr.feather --x datasets/processed/model_fitting/Mod_Q_Hab12a/temp_vp/VP_X.feather --gamma datasets/processed/model_fitting/Mod_Q_Hab12a/temp_vp/VP_Gamma.feather --output datasets/processed/model_fitting/Mod_Q_Hab12a/temp_vp/VP_A.feather --ncores 3 --chunk_size 50 >> datasets/processed/model_fitting/Mod_Q_Hab12a/temp_vp/VP_A.log 2>&1
python3 VP_geta.py --tr datasets/processed/model_fitting/Mod_Q_Hab12b/temp_vp/VP_Tr.feather --x datasets/processed/model_fitting/Mod_Q_Hab12b/temp_vp/VP_X.feather --gamma datasets/processed/model_fitting/Mod_Q_Hab12b/temp_vp/VP_Gamma.feather --output datasets/processed/model_fitting/Mod_Q_Hab12b/temp_vp/VP_A.feather --ncores 3 --chunk_size 50 >> datasets/processed/model_fitting/Mod_Q_Hab12b/temp_vp/VP_A.log 2>&1
python3 VP_getf.py --x datasets/processed/model_fitting/Mod_Q_Hab1/temp_vp/VP_X.feather --beta_dir datasets/processed/model_fitting/Mod_Q_Hab1/temp_vp --output datasets/processed/model_fitting/Mod_Q_Hab1/temp_vp/VP_F.feather --ncores 3 >> datasets/processed/model_fitting/Mod_Q_Hab1/temp_vp/VP_F.log 2>&1
python3 VP_getf.py --x datasets/processed/model_fitting/Mod_Q_Hab2/temp_vp/VP_X.feather --beta_dir datasets/processed/model_fitting/Mod_Q_Hab2/temp_vp --output datasets/processed/model_fitting/Mod_Q_Hab2/temp_vp/VP_F.feather --ncores 3 >> datasets/processed/model_fitting/Mod_Q_Hab2/temp_vp/VP_F.log 2>&1
python3 VP_getf.py --x datasets/processed/model_fitting/Mod_Q_Hab3/temp_vp/VP_X.feather --beta_dir datasets/processed/model_fitting/Mod_Q_Hab3/temp_vp --output datasets/processed/model_fitting/Mod_Q_Hab3/temp_vp/VP_F.feather --ncores 3 >> datasets/processed/model_fitting/Mod_Q_Hab3/temp_vp/VP_F.log 2>&1
python3 VP_getf.py --x datasets/processed/model_fitting/Mod_Q_Hab4a/temp_vp/VP_X.feather --beta_dir datasets/processed/model_fitting/Mod_Q_Hab4a/temp_vp --output datasets/processed/model_fitting/Mod_Q_Hab4a/temp_vp/VP_F.feather --ncores 3 >> datasets/processed/model_fitting/Mod_Q_Hab4a/temp_vp/VP_F.log 2>&1
python3 VP_getf.py --x datasets/processed/model_fitting/Mod_Q_Hab4b/temp_vp/VP_X.feather --beta_dir datasets/processed/model_fitting/Mod_Q_Hab4b/temp_vp --output datasets/processed/model_fitting/Mod_Q_Hab4b/temp_vp/VP_F.feather --ncores 3 >> datasets/processed/model_fitting/Mod_Q_Hab4b/temp_vp/VP_F.log 2>&1
python3 VP_getf.py --x datasets/processed/model_fitting/Mod_Q_Hab10/temp_vp/VP_X.feather --beta_dir datasets/processed/model_fitting/Mod_Q_Hab10/temp_vp --output datasets/processed/model_fitting/Mod_Q_Hab10/temp_vp/VP_F.feather --ncores 3 >> datasets/processed/model_fitting/Mod_Q_Hab10/temp_vp/VP_F.log 2>&1
python3 VP_getf.py --x datasets/processed/model_fitting/Mod_Q_Hab12a/temp_vp/VP_X.feather --beta_dir datasets/processed/model_fitting/Mod_Q_Hab12a/temp_vp --output datasets/processed/model_fitting/Mod_Q_Hab12a/temp_vp/VP_F.feather --ncores 3 >> datasets/processed/model_fitting/Mod_Q_Hab12a/temp_vp/VP_F.log 2>&1
python3 VP_getf.py --x datasets/processed/model_fitting/Mod_Q_Hab12b/temp_vp/VP_X.feather --beta_dir datasets/processed/model_fitting/Mod_Q_Hab12b/temp_vp --output datasets/processed/model_fitting/Mod_Q_Hab12b/temp_vp/VP_F.feather --ncores 3 >> datasets/processed/model_fitting/Mod_Q_Hab12b/temp_vp/VP_F.log 2>&1
```

  

»» Example VP_SLURM.slurm file

``` bash
#!/bin/bash
#SBATCH --job-name=VP_TF
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --account=project_465001588
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=1
#SBATCH --time=01:30:00
#SBATCH --partition=small-g
#SBATCH --output=datasets/processed/model_fitting/Mod_Q_Hab_TF/log/%x-%A-%a.out
#SBATCH --error=datasets/processed/model_fitting/Mod_Q_Hab_TF/log/%x-%A-%a.out
#SBATCH --array=1-24

# File containing commands to be executed
File=datasets/processed/model_fitting/Mod_Q_Hab_TF/VP_Commands.txt

# Load TensorFlow module and configure environment
ml use /appl/local/csc/modulefiles
ml tensorflow
export TF_CPP_MIN_LOG_LEVEL=3
export TF_ENABLE_ONEDNN_OPTS=0

# Verify GPU availability
python3 -c "import tensorflow as tf; print(\"Num GPUs Available:\", len(tf.config.list_physical_devices(\"GPU\")))"

# Run array job
head -n $SLURM_ARRAY_TASK_ID $File | tail -n 1 | bash

echo End of program at `date`
```

------------------------------------------------------------------------

## Step 2: GPU

In this step, latent factor predictions and variance partitioning are
computed on GPUs. The batch jobs for these computations can be submitted
using the `sbatch` command:

``` bash
# Submit SLURM jobs for variance partitioning and latent factor predictions
sbatch datasets/processed/model_fitting/Mod_Q_Hab_TF/VP_SLURM.slurm
sbatch datasets/processed/model_fitting/Mod_Q_Hab_TF/lf_SLURM.slurm
```

Additionally, cross-validated models are fitted by submitting the
corresponding SLURM scripts for each cross-validation strategy:

``` bash
# Submit SLURM jobs for cross-validated model fitting
#
# cross-validation method "cv_dist"
sbatch datasets/processed/model_fitting/HabX/model_fitting_cv/cv_bash_fit_dist.slurm
# cross-validation method "cv_large"
sbatch datasets/processed/model_fitting/HabX/model_fitting_cv/cv_bash_fit_large.slurm
```

------------------------------------------------------------------------

## Step 3: CPU

To continue the post-processing of the fitted models on CPUs, two
functions need to be executed:

### 1. **`mod_postprocess_2_cpu()`**

The
[`mod_postprocess_2_cpu()`](https://biodt.github.io/IASDT.R/reference/Mod_postprocessing.md)
function progresses the post-processing pipeline for HMSC models on the
CPU by automating the following tasks:

[TABLE]

  

### 2. **`mod_postprocess_cv_1_cpu()`**

The
[`mod_postprocess_cv_1_cpu()`](https://biodt.github.io/IASDT.R/reference/Mod_postprocessing.md)
function begins the post-processing of cross-validated models on the CPU
by automating the following tasks:

|                                                                                          |                                                                                                                                                                                                                                          |
|:-----------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`mod_merge_chains_cv()`](https://biodt.github.io/IASDT.R/reference/Mod_Merge_Chains.md) | merges fitted cross-validated model chains into `Hmsc` model objects and saves them to disk.                                                                                                                                             |
| [`predict_maps_cv()`](https://biodt.github.io/IASDT.R/reference/Predict_Maps.md)         | prepares scripts for predicting latent factors for each cross-validation strategy at new sampling units (evaluation folds). The arguments `lf_only` and `lf_commands_only` are set to `TRUE` to prepare only the necessary script files. |

Once
[`predict_maps_cv()`](https://biodt.github.io/IASDT.R/reference/Predict_Maps.md)
has been completed, the function combines the computation commands into
multiple text script files (tf_chunk\_\*.txt) in each model’s
model_fitting_cv/lf_TF_commands subdirectory. These scripts need to be
executed on GPUs using a single batch job submitted via a SLURM script
in the same directory, lf_SLURM.slurm.

------------------------------------------------------------------------

## Step 4: GPU

In this step, the computation of latent factors for cross-validated
models is performed on GPUs using SLURM scripts.

``` bash
sbatch datasets/processed/model_fitting/HabX/model_fitting_cv/lf_TF_commands/lf_SLURM.slurm
```

------------------------------------------------------------------------

## Step 5: CPU

The final step of the post-processing pipeline is carried out on CPUs
using the
[`mod_postprocess_cv_2_cpu()`](https://biodt.github.io/IASDT.R/reference/Mod_postprocessing.md)
function. This function automates the following tasks:

- Predicting habitat suitability at the testing cross-validation folds
  using the
  [`predict_maps_cv()`](https://biodt.github.io/IASDT.R/reference/Predict_Maps.md)
  function.
- Computing the model’s predictive power (using spatially independent
  testing data) with the same function, based on four metrics: AUC (area
  under the ROC curve); RMSE (root mean square error); continuous Boyce
  index; and Tjur R².
- Plotting the model’s evaluation, including:
  - Predictive power values for each evaluation metric versus the mean
    number of testing presences
  - Explanatory versus predictive power for each evaluation metric

------------------------------------------------------------------------

**Previous articles:**  
↠[1.
Overview](https://biodt.github.io/IASDT.R/articles/workflow_1_overview.md)  
↠[2. Processing abiotic
data](https://biodt.github.io/IASDT.R/articles/workflow_2_abiotic_data.md)  
↠[3. Processing biotic
data](https://biodt.github.io/IASDT.R/articles/workflow_3_biotic_data.md)  
↠[4. Model
fitting](https://biodt.github.io/IASDT.R/articles/workflow_4_model_fitting.md)  
