#!/bin/bash
#SBATCH --job-name=pipedream_array
#SBATCH --partition=cs05r
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=5120
#SBATCH --array=0-<N>
#SBATCH --output=pipedream_array_%A_%a.out

# Load modules
source /etc/profile
module load buster
module load graphviz
export LD_LIBRARY_PATH=/dls_sw/apps/graphviz/rhel8/12.2/lib:$LD_LIBRARY_PATH

# Path to your datasets CSV and YAML
DATASETS_CSV="/path/to/datasets.csv"
YAML_FILE="/path/to/pipedream_parameters.yaml"

# Get the dataset info for this array task
TASK_ID=$SLURM_ARRAY_TASK_ID
# Extract the dataset row (skip header)
DATASET_LINE=$(awk -v n=$((TASK_ID+2)) 'NR==n' "$DATASETS_CSV")

# Parse CSV columns (adjust if your CSV has different columns/order)
IFS=',' read -r CrystalName CompoundCode RefinementMTZfree DimplePathToPDB DimplePathToMTZ RefinementCIF <<< "$DATASET_LINE"

# Set up input/output directories as in your Python logic
PROCESSING_DIR="/path/to/processing"
DATASET_DIR="$PROCESSING_DIR/$CrystalName"
INPUT_DIR="$DATASET_DIR/input_files"
OUTPUT_DIR="$DATASET_DIR/Pipedream_array_$TASK_ID"
mkdir -p "$INPUT_DIR" "$OUTPUT_DIR"

# Copy input files (adjust as needed)
cp "$PROCESSING_DIR/analysis/model_building/$CrystalName/$RefinementMTZfree" "$INPUT_DIR/"
cp "$DimplePathToPDB" "$INPUT_DIR/"
cp "$DimplePathToMTZ" "$INPUT_DIR/"
if [ -n "$RefinementCIF" ]; then
  cp "$RefinementCIF" "$INPUT_DIR/"
fi

# Build the Pipedream command (add other options as needed)
CMD="/dls_sw/apps/GPhL/BUSTER/20240123/scripts/pipedream -nolmr \
  -hklin $INPUT_DIR/$(basename $RefinementMTZfree) \
  -xyzin $INPUT_DIR/$(basename $DimplePathToPDB) \
  -hklref $INPUT_DIR/$(basename $DimplePathToMTZ) \
  -d $OUTPUT_DIR"

# Run Pipedream
srun $CMD
