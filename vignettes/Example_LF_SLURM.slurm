#!/bin/bash
#SBATCH --job-name=PP_lf
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --account=project_465001588
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=1
#SBATCH --time=01:00:00
#SBATCH --partition=small-g
#SBATCH --output=datasets/processed/model_fitting/TMP/%x-%A-%a.out
#SBATCH --error=datasets/processed/model_fitting/TMP/%x-%A-%a.out
#SBATCH --array=1-186

# Define directories
OutputDir="datasets/processed/model_fitting/TF_BatchFiles"

# Find all the split files and sort them explicitly
SplitFiles=($(find "$OutputDir" -type f -name "tf_chunk_*.txt" | sort -V))

# Check if files were found
if [ ${#SplitFiles[@]} -eq 0 ]; then
    echo "Error: No files matching tf_chunk_*.txt found in $OutputDir"
    exit 1
fi

# Ensure no more than `, NumFiles, ` files are processed
MaxFiles=186
if [ ${#SplitFiles[@]} -gt $MaxFiles ]; then
    SplitFiles=("${SplitFiles[@]:0:$MaxFiles}")
    echo "More than $MaxFiles files found, limiting to the first $MaxFiles files."
fi

# Get the index of the current task based on SLURM_ARRAY_TASK_ID
TaskIndex=$((SLURM_ARRAY_TASK_ID - 1))

# Validate TaskIndex
if [ $TaskIndex -ge ${#SplitFiles[@]} ] || [ $TaskIndex -lt 0 ]; then
    echo "Error: TaskIndex $TaskIndex is out of range. Valid range: 0 to $((${#SplitFiles[@]} - 1))"
    exit 1
fi

# Get the specific split file to process based on the job array task ID
SplitFile="${SplitFiles[$TaskIndex]}"

# Verify the selected split file
if [ -z "$SplitFile" ] || [ ! -f "$SplitFile" ]; then
    echo "Error: File $SplitFile does not exist or is invalid."
    exit 1
fi

# Processing file
echo "Processing file: $SplitFile"

# Run the selected split file
bash "$SplitFile"

echo End of program at `date`
