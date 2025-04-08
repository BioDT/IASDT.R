#!/bin/bash
#SBATCH --job-name=VP_TF
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --account=project_465001588
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=1
#SBATCH --time=01:30:00
#SBATCH --partition=small-g
#SBATCH --output=datasets/processed/model_fitting/TF_postprocess/log/%x-%A-%a.out
#SBATCH --error=datasets/processed/model_fitting/TF_postprocess/log/%x-%A-%a.out
#SBATCH --array=1-24

# File containing commands to be executed
File=datasets/processed/model_fitting/VP_Commands.txt

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
