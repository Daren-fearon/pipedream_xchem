import os
import sys
import paramiko  # type: ignore
import time
import yaml  # type: ignore
import pandas as pd  # type: ignore
import sqlite3
import logging
from datetime import datetime
import argparse
import shutil
import json

# Constants
CLUSTER_BASTION = "wilson.diamond.ac.uk"
CLUSTER_USER = os.environ.get("CLUSTER_USER", os.getlogin())

def parse_args():
    parser = argparse.ArgumentParser(description="Submit SLURM jobs for Pipedream analysis.")
    parser.add_argument(
        "-y", "--yaml",
        type=str,
        required=True,
        help="Path to the pipedream_parameters.yaml file"
    )
    return parser.parse_args()

def read_yaml(yaml_file):
    try:
        with open(yaml_file, 'r') as file:
            params = yaml.safe_load(file)
            logging.info("YAML file read successfully.")
            return params
    except Exception as e:
        logging.error(f"Error reading YAML file: {e}")
        raise

def validate_params(params):
    required_keys = ['Database_path', 'Mode', 'Processing_directory', 'Refinement_parameters']
    refinement_keys = ['keepwater', 'mrefine', 'remediate', 'sidechainrebuild', 'runpepflip', 'rhocommands']

    for key in required_keys:
        if key not in params:
            raise ValueError(f"Missing required key in YAML parameters: {key}")

    refinement_params = params.get('Refinement_parameters', {})
    for key in refinement_keys:
        if key not in refinement_params:
            raise ValueError(f"Missing required refinement parameter: {key}")

    if params['Mode'] == 'specific_datasets' and 'Dataset_csv_path' not in params:
        raise ValueError("Missing 'Dataset_csv_path' for 'specific_datasets' mode")

def get_datasets(params):
    try:
        conn = sqlite3.connect(params['Database_path'])

        if params['Mode'] == 'specific_datasets':
            # Only require CrystalName in the CSV
            csv_path = params.get('Dataset_csv_path') or params.get('Filtered_dataset_csv_path')
            if not csv_path:
                raise ValueError("Missing 'Dataset_csv_path' for 'specific_datasets' mode")
            csv_df = pd.read_csv(csv_path)
            if 'CrystalName' not in csv_df.columns:
                raise ValueError("CSV must contain a 'CrystalName' column.")
            crystal_names = csv_df['CrystalName'].dropna().unique().tolist()
            if not crystal_names:
                raise ValueError("No CrystalName values found in the CSV.")

            query = f"""
                SELECT CrystalName, CompoundCode, RefinementMTZfree, 
                       DimplePathToPDB, DimplePathToMTZ, RefinementCIF 
                FROM mainTable 
                WHERE CrystalName IN ({','.join(['?']*len(crystal_names))})
            """
            datasets = pd.read_sql(query, conn, params=crystal_names)
            # Do not overwrite the original CSV
        else:
            query = """
                SELECT CrystalName, CompoundCode, RefinementMTZfree, 
                       DimplePathToPDB, DimplePathToMTZ, RefinementCIF 
                FROM mainTable 
                WHERE RefinementOutcome = '1 - Analysis Pending'
            """
            datasets = pd.read_sql(query, conn)
            # Save to CSV for array job compatibility
            csv_path = os.path.join(params['Processing_directory'], 'datasets_metadata.csv')
            datasets.to_csv(csv_path, index=False)
            params['Filtered_dataset_csv_path'] = csv_path

        conn.close()

        # Check for missing or empty required fields
        required_fields = ['CompoundCode', 'RefinementMTZfree', 'DimplePathToPDB', 'DimplePathToMTZ', 'RefinementCIF']
        missing_info = []
        for idx, row in datasets.iterrows():
            missing = [field for field in required_fields if pd.isna(row[field]) or str(row[field]).strip() == '']
            if missing:
                missing_info.append({
                    'CrystalName': row.get('CrystalName', f'row {idx}'),
                    'MissingFields': missing
                })
                logging.warning(f"CrystalName: {row.get('CrystalName', f'row {idx}')} | Missing: {', '.join(missing)}")
        # Do not raise error, just log
        if missing_info:
            msg = ("One or more datasets are missing required fields. "
                   "These will be processed, but may fail downstream.\n" +
                   "Details:\n" +
                   '\n'.join([f"CrystalName: {item['CrystalName']} | Missing: {', '.join(item['MissingFields'])}" for item in missing_info]))
            logging.warning(msg)

        logging.info(f"{len(datasets)} datasets retrieved successfully.")
        return datasets

    except Exception as e:
        logging.error(f"Error retrieving datasets: {e}")
        raise

def fetch_password():
    import getpass
    if not sys.stdin.isatty():
        logging.error("Standard input is not a TTY. Password input may not be secure.")
        raise RuntimeError("Insecure environment for password input.")
    return getpass.getpass("Enter your password: ")

def submit_jobs(datasets, params):
    """
    Prepares input files and directories for each dataset. All job submission is handled via a single SLURM array job script.
    Assigns a unique Pipedream_N output directory for each dataset at preparation time.
    """
    try:
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        output_base = params.get('Output_directory') or params['Processing_directory']
        parent_dir = params.get('Output_directory') or f"{params['Processing_directory']}/analysis/Pipedream/Pipedream_{timestamp}"
        os.makedirs(parent_dir, exist_ok=True)

        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)

        logging.basicConfig(
            filename=f'{parent_dir}/pipedream.log',
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )

        output_yaml = {}
        filtered_datasets = []

        for index, row in datasets.iterrows():
            crystal_name = row['CrystalName']
            compound_code = row['CompoundCode']
            refinement_mtzfree = row['RefinementMTZfree']
            dimple_path_pdb = row['DimplePathToPDB']
            dimple_path_mtz = row['DimplePathToMTZ']
            refinement_cif = row.get('RefinementCIF')
            use_rhofit = bool(refinement_cif and str(refinement_cif).strip())

            required_paths = {
                "RefinementMTZfree": refinement_mtzfree,
                "DimplePathToPDB": dimple_path_pdb,
                "DimplePathToMTZ": dimple_path_mtz
            }
            missing_keys = [k for k, v in required_paths.items() if v is None or str(v).strip() == ""]
            if missing_keys:
                logging.error(f"Skipping dataset {crystal_name} due to missing required file paths: {missing_keys} | Values: {required_paths}")
                continue

            # Construct the full path for RefinementMTZfree
            refinement_mtzfree_path = os.path.realpath(os.path.join(params['Processing_directory'], "analysis", "model_building", crystal_name, refinement_mtzfree))
            dimple_path_pdb = os.path.realpath(dimple_path_pdb) if dimple_path_pdb else None
            dimple_path_mtz = os.path.realpath(dimple_path_mtz) if dimple_path_mtz else None

            if use_rhofit:
                if refinement_cif is None or str(refinement_cif).strip() == "":
                    logging.error(f"Refinement CIF is missing or empty for {crystal_name}. Skipping.")
                    logging.info(f"Skipping dataset {crystal_name} due to missing refinement_cif")
                    continue
                refinement_cif = resolve_dataset_path(refinement_cif, crystal_name, params['Processing_directory'])

            dataset_dir = f"{parent_dir}/{crystal_name}"
            input_dir = f"{dataset_dir}/input_files"
            filtered_datasets.append(row)

        # Prepare input files and output_yaml for filtered datasets only
        output_dir_list = []
        for row in filtered_datasets:
            crystal_name = row['CrystalName']
            compound_code = row['CompoundCode']
            refinement_mtzfree = row['RefinementMTZfree']
            dimple_path_pdb = row['DimplePathToPDB']
            dimple_path_mtz = row['DimplePathToMTZ']
            refinement_cif = row.get('RefinementCIF')
            use_rhofit = bool(refinement_cif and str(refinement_cif).strip())

            refinement_mtzfree_path = os.path.realpath(os.path.join(params['Processing_directory'], "analysis", "model_building", crystal_name, refinement_mtzfree))
            dimple_path_pdb = os.path.realpath(dimple_path_pdb) if dimple_path_pdb else None
            dimple_path_mtz = os.path.realpath(dimple_path_mtz) if dimple_path_mtz else None

            if use_rhofit:
                refinement_cif = resolve_dataset_path(refinement_cif, crystal_name, params['Processing_directory'])

            dataset_dir = f"{parent_dir}/{crystal_name}"
            input_dir = f"{dataset_dir}/input_files"
            os.makedirs(input_dir, exist_ok=True)

            # Find next available Pipedream_N output dir for this dataset
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_dir = f"{dataset_dir}/Pipedream_{timestamp}"
            output_dir_list.append(output_dir)

            try:
                if params.get("Copy_files_to_working_dir", True):
                    shutil.copy(refinement_mtzfree_path, os.path.join(input_dir, os.path.basename(refinement_mtzfree_path)))
                    shutil.copy(dimple_path_pdb, os.path.join(input_dir, os.path.basename(dimple_path_pdb)))
                    shutil.copy(dimple_path_mtz, os.path.join(input_dir, os.path.basename(dimple_path_mtz)))
                    if use_rhofit:
                        shutil.copy(refinement_cif, os.path.join(input_dir, os.path.basename(refinement_cif)))
                        # Check for <CompoundCode>.png in the same directory as RefinementCIF
                        if refinement_cif:
                            png_path = os.path.join(os.path.dirname(refinement_cif), f"{compound_code}.png")
                            if os.path.exists(png_path):
                                shutil.copy(png_path, os.path.join(input_dir, os.path.basename(png_path)))
                        else:
                            logging.warning(f"refinement_cif is None for {crystal_name} when checking for PNG.")

                # Remove crystallisation components if requested
                if params.get("Remove_crystallisation_components", False):
                    pdb_path = os.path.join(input_dir, os.path.basename(dimple_path_pdb))
                    if pdb_path and os.path.exists(pdb_path):
                        with open(pdb_path, 'r') as f:
                            lines = f.readlines()
                        with open(pdb_path, 'w') as f:
                            for line in lines:
                                if not any(res in line for res in ["DMS", "EDO", "GOL", "SO4", "PO4", "PEG"]):
                                    f.write(line)
                    else:
                        logging.warning(f"PDB path missing or does not exist for {crystal_name}: {pdb_path}")

            except FileNotFoundError as e:
                logging.error(f"Missing file for {crystal_name}: {e}")
                logging.info(f"Skipping dataset {crystal_name} due to missing file: {e}")
                continue
            except TypeError as e:
                logging.error(f"TypeError for {crystal_name}: {e}. Variables: refinement_mtzfree_path={refinement_mtzfree_path}, dimple_path_pdb={dimple_path_pdb}, dimple_path_mtz={dimple_path_mtz}, refinement_cif={refinement_cif}")
                logging.info(f"Skipping dataset {crystal_name} due to TypeError: {e}")
                continue

            # --- Build correct Pipedream command for SLURM script ---
            hklin_file = os.path.join(input_dir, os.path.basename(refinement_mtzfree_path))
            xyzin_file = os.path.join(input_dir, os.path.basename(dimple_path_pdb))
            hklref_file = os.path.join(input_dir, os.path.basename(dimple_path_mtz))
            rhofit_file = os.path.join(input_dir, os.path.basename(refinement_cif)) if use_rhofit else None
            rhocommands = params['Refinement_parameters'].get('rhocommands', None)
            pipedream_cmd = f"/dls_sw/apps/GPhL/BUSTER/20240123/scripts/pipedream -nolmr -hklin {hklin_file} -xyzin {xyzin_file} -hklref {hklref_file} -d {output_dir} "
            refinement_args = refinement_params_to_args(params.get('Refinement_parameters', {}), rhofit_file=rhofit_file, rhocommands=rhocommands if use_rhofit else None)
            if refinement_args:
                pipedream_cmd += refinement_args + " "
            # Store the command for later use if needed (not used in this function, but could be logged)
            output_yaml[crystal_name] = {
                'CrystalName': crystal_name,
                'CompoundCode': compound_code,
                'PipedreamDirectory': output_dir,
                'ReportHTML': f"{output_dir}/report-{compound_code}/index.html",
                'LigandReportHTML': f"{output_dir}/report-{compound_code}/ligand/index.html",
                'ExpectedSummary': f"{output_dir}/pipedream_summary.json",
                'PipedreamCommand': pipedream_cmd.strip()
            }

            if use_rhofit and os.path.exists(os.path.join(input_dir, f"{compound_code}.png")):
                output_yaml[crystal_name]['LigandPNG'] = os.path.join(input_dir, f"{compound_code}.png")

        # Save output_yaml to file as JSON
        output_json_name = f"Pipedream_{timestamp}_output.json"
        with open(f"{parent_dir}/{output_json_name}", 'w') as json_file:
            json.dump(output_yaml, json_file, indent=2)

        # Save filtered datasets to CSV for array job compatibility, including output_dir and full command
        filtered_df = pd.DataFrame(filtered_datasets)
        filtered_df['PipedreamDirectory'] = output_dir_list
        # Add the full command for each dataset
        filtered_df['PipedreamCommand'] = [output_yaml[row['CrystalName']]['PipedreamCommand'] for row in filtered_datasets]
        # Only keep the required columns for the SLURM script
        required_cols = ['CrystalName', 'CompoundCode', 'RefinementMTZfree', 'DimplePathToPDB', 'DimplePathToMTZ', 'RefinementCIF', 'PipedreamDirectory', 'PipedreamCommand']
        filtered_df = filtered_df[required_cols]
        filtered_csv_path = os.path.join(parent_dir, 'datasets_metadata.csv')
        filtered_df.to_csv(filtered_csv_path, index=False)
        params['Filtered_dataset_csv_path'] = filtered_csv_path

        # Print the number of valid datasets for SLURM
        print(f"Number of valid datasets for SLURM: {filtered_df.shape[0]}")
        logging.info(f"Number of valid datasets for SLURM: {filtered_df.shape[0]}")

        # Write SLURM array index to CrystalName mapping for debugging in array_logs dir
        array_logs_dir = os.path.join(parent_dir, 'array_logs')
        os.makedirs(array_logs_dir, exist_ok=True)
        # Use timestamp for unique map file name
        map_timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        slurm_map_path = os.path.join(array_logs_dir, f'slurm_array_index_map_{map_timestamp}.csv')
        with open(slurm_map_path, 'w') as f:
            f.write('SLURM_ARRAY_TASK_ID,CrystalName\n')
            for idx, row in enumerate(filtered_datasets):
                f.write(f'{idx},{row["CrystalName"]}\n')

        return parent_dir, filtered_df.shape[0]

    except Exception as e:
        logging.error(f"Error preparing input files: {e}")
        raise

def refinement_params_to_args(refinement_params, rhofit_file=None, rhocommands=None):
    """
    Convert the Refinement_parameters dict from YAML to a string of Pipedream command-line arguments.
    Handles booleans, strings, and lists (e.g. rhocommands).
    If rhofit_file and rhocommands are provided, include them as well.
    Automatically prepends a dash to each rhocommand if not present, and ensures correct command line for -rhocommands.
    Handles both single and multiple rhocommands robustly.
    """
    args = []
    for key, value in refinement_params.items():
        if key == 'rhocommands':
            continue  # handled separately
        if isinstance(value, bool):
            if value:
                args.append(f"-{key}")
        elif isinstance(value, str):
            if value.strip():
                args.append(f"-{key} {value}")
        elif isinstance(value, list):
            # For lists (other than rhocommands), join as space-separated
            joined = ' '.join(str(v) for v in value if str(v).strip())
            if joined:
                args.append(f"-{key} {joined}")
    if rhofit_file:
        args.append(f"-rhofit {rhofit_file}")
    if rhocommands:
        # Accept both list and string
        def ensure_dash(cmd):
            cmd = str(cmd).strip()
            return cmd if cmd.startswith('-') else f'-{cmd}' if cmd else ''
        if isinstance(rhocommands, list):
            rhocmds = [ensure_dash(cmd) for cmd in rhocommands if str(cmd).strip()]
            rhocmd_str = ' '.join(rhocmds)
        else:
            rhocmd_str = ensure_dash(rhocommands)
        if rhocmd_str:
            # Remove any leading/trailing quotes to avoid double quoting
            rhocmd_str = rhocmd_str.strip('"\'')
            # Collapse multiple spaces
            rhocmd_str = ' '.join(rhocmd_str.split())
            args.append(f'-rhocommands {rhocmd_str}')
    return " ".join(args)

def write_array_job_script(datasets_csv, output_dir, num_datasets, script_path, refinement_params=None):
    # Only write array line if there is at least one dataset
    if num_datasets < 1:
        raise ValueError("No datasets to process: cannot write SLURM array job script.")
    # Cap at 100 concurrent jobs
    max_concurrent = 100
    script = f"""#!/bin/bash
#SBATCH --job-name=pipedream_array
#SBATCH --partition=cs05r
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4096
#SBATCH --array=0-{num_datasets-1}%{max_concurrent}
#SBATCH --output=array_logs/pipedream_array_%A_%a.out

# Create array_logs directory for SLURM output
mkdir -p array_logs

# Load modules (silence module system errors)
source /etc/profile 2>/dev/null
module load buster >/dev/null 2>&1
module load graphviz >/dev/null 2>&1
export LD_LIBRARY_PATH=/dls_sw/apps/graphviz/rhel8/12.2/lib:$LD_LIBRARY_PATH

# Path to your datasets CSV
DATASETS_CSV=\"{datasets_csv}\"
# Path to processing directory for input files
PROCESSING_DIR=\"{params['Processing_directory']}\"

# Get the dataset info for this array task
TASK_ID=$SLURM_ARRAY_TASK_ID
# Extract the dataset row (skip header)
DATASET_LINE=$(awk -v n=$((TASK_ID+2)) 'NR==n' "$DATASETS_CSV")

# Check for empty or invalid dataset line
if [ -z "$DATASET_LINE" ] || [[ "$DATASET_LINE" =~ ^[0-9]+$ ]]; then
  echo "[ERROR] No valid dataset line found for SLURM_ARRAY_TASK_ID $TASK_ID. Exiting."
  exit 3
fi

# Parse CSV columns (adjust if your CSV has different columns/order)
IFS=',' read -r CrystalName CompoundCode RefinementMTZfree DimplePathToPDB DimplePathToMTZ RefinementCIF PipedreamDirectory PipedreamCommand <<< "$DATASET_LINE"

# Check for invalid or missing CrystalName
if [ -z "$CrystalName" ] || [[ "$CrystalName" =~ ^[0-9]+$ ]]; then
  echo "[ERROR] Invalid or missing CrystalName for SLURM_ARRAY_TASK_ID $TASK_ID. Exiting."
  exit 4
fi

# Print CrystalName at the top of the log for easier debugging
echo "[SLURM] CrystalName: $CrystalName (SLURM_ARRAY_TASK_ID: $TASK_ID)"

# Set up input/output directories as in your Python logic
OUTPUT_DIR_BASE=\"{output_dir}\"
DATASET_DIR="$OUTPUT_DIR_BASE/$CrystalName"
INPUT_DIR="$DATASET_DIR/input_files"
OUTPUT_DIR="$PipedreamDirectory"

# Redirect all output to a unique per-dataset log file in the dataset directory (using output dir name)
mkdir -p "$DATASET_DIR"
LOG_BASENAME="${{CrystalName}}_slurm_$(basename $OUTPUT_DIR).out"
exec > "$DATASET_DIR/$LOG_BASENAME" 2>&1

# Check for empty output dir
if [ -z "$OUTPUT_DIR" ]; then
  echo "[ERROR] OUTPUT_DIR is empty for $CrystalName. Check your datasets.csv and Python logic."
  exit 2
fi

# Only create input dir, not output dir
mkdir -p "$INPUT_DIR"

# If output dir exists (shouldn't, but double check), skip
if [ -d "$OUTPUT_DIR" ]; then
  echo "[INFO] Skipping $CrystalName: output directory $OUTPUT_DIR already exists."
  exit 0
fi

# Copy input files from processing directory (not output directory)
cp "$PROCESSING_DIR/analysis/model_building/$CrystalName/$RefinementMTZfree" "$INPUT_DIR/"
cp "$DimplePathToPDB" "$INPUT_DIR/"
cp "$DimplePathToMTZ" "$INPUT_DIR/"
if [ -n "$RefinementCIF" ]; then
  cp "$RefinementCIF" "$INPUT_DIR/"
fi

# Use the full command from the CSV (already quoted as needed)
CMD="$PipedreamCommand"

# Collapse multiple spaces to a single space
CMD=$(echo "$CMD" | tr -s ' ')

# Run Pipedream (this will create the output dir)
echo "$CMD"
eval "$CMD"

# Move any __*.setvar.lis files to array_logs dir for tidiness
if compgen -G "__*.setvar.lis" > /dev/null; then
  mv __*.setvar.lis array_logs/
fi
"""
    with open(script_path, "w") as f:
        f.write(script)

def submit_sbatch_on_wilson(script_path, processing_dir):
    """
    SSH to wilson and submit the job script using sbatch, using password authentication.
    Uses $CLUSTER_USER if set, otherwise current user.
    Uses paramiko for SSH and SFTP.
    """
    import getpass
    user = os.environ.get("CLUSTER_USER", os.getlogin())
    wilson_host = "wilson.diamond.ac.uk"
    remote_script_name = os.path.basename(script_path)
    remote_script_path = os.path.join(processing_dir, remote_script_name)
    password = getpass.getpass(f"Enter password for {user}@wilson: ")
    print(f"[INFO] Connecting to wilson and copying job script...")
    try:
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(wilson_host, username=user, password=password)
        sftp = ssh.open_sftp()
        # Ensure remote directory exists
        try:
            sftp.chdir(processing_dir)
        except IOError:
            # Directory does not exist, try to create it
            parts = processing_dir.strip('/').split('/')
            path = ''
            for part in parts:
                path += '/' + part
                try:
                    sftp.chdir(path)
                except IOError:
                    sftp.mkdir(path)
                    sftp.chdir(path)
        # Copy the script
        sftp.put(script_path, remote_script_path)
        sftp.close()
        print(f"[INFO] Script copied to wilson:{remote_script_path}")
        # Submit the job
        stdin, stdout, stderr = ssh.exec_command(f"cd {processing_dir} && sbatch {remote_script_name}")
        out = stdout.read().decode()
        err = stderr.read().decode()
        print(f"[INFO] sbatch output:\n{out}")
        # Only print sbatch error if it is not the known harmless module error
        harmless_module_error = "Module ERROR: no such variable"
        if err and harmless_module_error not in err:
            print(f"[ERROR] sbatch error:\n{err}")
        else:
            if err:
                # Log the harmless error for reference
                logging.info(f"[sbatch harmless stderr suppressed from terminal]: {err.strip()}")
            print("[INFO] Job script submitted on wilson via sbatch.")
        ssh.close()
    except Exception as e:
        print(f"[ERROR] Failed to submit job on wilson: {e}")
        print(f"You can submit manually with:\n  scp {script_path} {user}@wilson:{remote_script_path}\n  ssh {user}@wilson 'cd {processing_dir} && sbatch {remote_script_name}'")

def resolve_dataset_path(path, crystal_name, processing_dir):
    """
    Resolves a dataset path, handling both absolute and relative paths.
    If relative, constructs the full path using the processing directory and crystal name.
    """
    if os.path.isabs(path):
        return os.path.realpath(path)
    else:
        return os.path.realpath(os.path.join(
            processing_dir, "analysis", "model_building", crystal_name, path
        ))

if __name__ == "__main__":
    try:
        args = parse_args()
        yaml_file = args.yaml
        params = read_yaml(yaml_file)
        validate_params(params)
        datasets = get_datasets(params)
        if datasets is not None and not datasets.empty:
            parent_dir, num_valid_datasets = submit_jobs(datasets, params)
            # Write and submit array job script
            datasets_csv = params['Filtered_dataset_csv_path']
            output_dir = params.get('Output_directory') or params['Processing_directory']
            script_path = os.path.join(output_dir, 'pipedream_array_job.sh')
            write_array_job_script(datasets_csv, output_dir, num_valid_datasets, script_path, params.get('Refinement_parameters'))
            print(f"Submitting SLURM array job for {num_valid_datasets} datasets...")
            submit_sbatch_on_wilson(script_path, output_dir)
        else:
            logging.warning("No datasets to process. Exiting.")
    except Exception as e:
        logging.error(f"Script execution failed: {e}")

