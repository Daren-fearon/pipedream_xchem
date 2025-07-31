"""
Pipedream XChem Integration Script

This script submits SLURM jobs for Pipedream analysis with ligand restraints 
generated from SMILES using grade2. It processes crystallographic datasets 
from an SQLite database and prepares them for automated refinement.

Author: DFearon
Date: July 2025
"""

import os
import sys
import getpass
import paramiko  # type: ignore
import yaml  # type: ignore
import pandas as pd  # type: ignore
import sqlite3
import logging
from datetime import datetime
import argparse
import shutil
import json
from typing import Dict, List, Optional, Tuple, Any

# Constants
CLUSTER_BASTION = "wilson.diamond.ac.uk"
CLUSTER_USER = os.environ.get("CLUSTER_USER", os.getlogin())
MAX_CONCURRENT_JOBS = 100
VERSION = "1.0.0"


def setup_logging(log_dir: str = None, log_level: str = "INFO", verbose: bool = False) -> None:
    """
    Configure logging for the application.
    
    Args:
        log_dir: Directory to save log files. If None, uses current directory.
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR)
        verbose: Enable verbose console output
    """
    if log_dir:
        os.makedirs(log_dir, exist_ok=True)
        log_file = os.path.join(log_dir, 'pipedream.log')
    else:
        log_file = 'pipedream.log'
    
    # Convert string level to logging constant
    numeric_level = getattr(logging, log_level.upper(), logging.INFO)
    
    # Clear any existing handlers
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    
    handlers = [logging.FileHandler(log_file)]
    
    # Add console handler based on verbose setting
    if verbose:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(numeric_level)
        handlers.append(console_handler)
    else:
        # Minimal console output for non-verbose mode
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.WARNING)  # Only warnings and errors to console
        handlers.append(console_handler)
    
    logging.basicConfig(
        level=numeric_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=handlers,
        force=True  # Force reconfiguration
    )


def safe_file_copy(src: str, dst: str, crystal_name: str = "Unknown") -> bool:
    """
    Safely copy a file with proper error handling and logging.
    
    Args:
        src: Source file path
        dst: Destination file path
        crystal_name: Crystal name for logging context
        
    Returns:
        True if successful, False otherwise
    """
    try:
        shutil.copy2(src, dst)
        logging.debug(f"Copied {src} -> {dst} for crystal {crystal_name}")
        return True
    except Exception as e:
        logging.error(f"Failed to copy {src} -> {dst} for crystal {crystal_name}: {e}")
        return False


def validate_file_paths(paths: Dict[str, str], crystal_name: str = "Unknown") -> List[str]:
    """
    Validate that required file paths exist.
    
    Args:
        paths: Dictionary of path descriptions to file paths
        crystal_name: Crystal name for logging context
        
    Returns:
        List of missing file descriptions
    """
    missing_files = []
    for description, path in paths.items():
        if not path or not os.path.exists(path):
            missing_files.append(f"{description}: {path}")
            logging.warning(f"Missing file for crystal {crystal_name} - {description}: {path}")
    return missing_files


def extract_dataset_fields(row: pd.Series) -> Tuple[str, str, str, str, str, str]:
    """Extract required dataset fields from a pandas Series row."""
    return (
        row['CrystalName'],
        row['CompoundCode'],
        row['RefinementMTZfree'],
        row['DimplePathToPDB'],
        row['DimplePathToMTZ'],
        row['CompoundSMILES']
    )

def resolve_path(base_dir: str, *parts: str) -> str:
    """Resolve a file path from base directory and path components."""
    return os.path.realpath(os.path.join(base_dir, *parts))

def get_smiles_for_compound(compound_code: str, smiles_from_db: Optional[str], input_dir: str) -> Optional[str]:
    """
    Get SMILES string for a compound from database and save it as a file.
    
    Args:
        compound_code (str): Compound code
        smiles_from_db (str): SMILES from database (may be None/empty)
        input_dir (str): Directory to save SMILES file (required)
    
    Returns:
        str: SMILES string or None if not found
    """
    # Get SMILES from database
    if smiles_from_db and str(smiles_from_db).strip():
        smiles_string = str(smiles_from_db).strip()
        logging.info(f"Using SMILES from database for {compound_code}: {smiles_string}")
        
        # Always save SMILES to file in input directory
        os.makedirs(input_dir, exist_ok=True)
        smiles_file_path = os.path.join(input_dir, f"{compound_code}.smiles")
        try:
            with open(smiles_file_path, 'w', encoding='utf-8') as f:
                f.write(f"{smiles_string}\n")
            logging.info(f"Saved SMILES to file: {smiles_file_path}")
        except Exception as e:
            logging.warning(f"Failed to save SMILES file for {compound_code}: {e}")
        
        return smiles_string
    
    logging.warning(f"No SMILES found in database for {compound_code}")
    return None

def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Submit SLURM jobs for Pipedream analysis with ligand restraints generated from SMILES using grade2.",
        epilog="Example: python pipedream_xchem.py --parameters parameters.yaml --output-json /path/to/results.json --dry-run",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Version
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {VERSION}"
    )
    
    # Required arguments
    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument(
        "--parameters", "-p",
        type=str,
        required=True,
        help="Path to the parameters file (pipedream_parameters.yaml)"
    )
    
    # Output options
    output_group = parser.add_argument_group('output options')
    output_group.add_argument(
        "--output-json",
        type=str,
        help="Custom output path for the JSON metadata file (default: auto-generated with timestamp in output directory)"
    )
    output_group.add_argument(
        "--output-csv", 
        type=str,
        help="Custom output path for the CSV datasets file (default: datasets_metadata.csv in output directory)"
    )
    
    # Execution options
    execution_group = parser.add_argument_group('execution options')
    execution_group.add_argument(
        "--dry-run",
        action="store_true",
        help="Prepare files and generate commands but don't submit SLURM job (useful for testing)"
    )
    
    # Logging options
    logging_group = parser.add_argument_group('logging options')
    logging_group.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose console output"
    )
    logging_group.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Set logging level (default: INFO)"
    )
    
    args = parser.parse_args()
    
    # Validate custom output paths
    if args.output_json and not args.output_json.endswith('.json'):
        parser.error("--output-json must have a .json extension")
    if args.output_csv and not args.output_csv.endswith('.csv'):
        parser.error("--output-csv must have a .csv extension")
    
    return args

def read_yaml(yaml_file: str) -> Dict[str, Any]:
    """Read and parse a YAML configuration file."""
    try:
        with open(yaml_file, 'r', encoding='utf-8') as file:
            params = yaml.safe_load(file)
            logging.info("YAML file read successfully.")
            return params
    except Exception as e:
        logging.error(f"Error reading YAML file: {e}")
        raise

def validate_params(params: Dict[str, Any]) -> None:
    """Validate required parameters in the configuration."""
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

def get_datasets(params: Dict[str, Any]) -> pd.DataFrame:
    """Retrieve datasets from the database based on configuration parameters."""
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

            # Try to include CompoundSMILES, but don't fail if column doesn't exist
            try:
                query = f"""
                    SELECT CrystalName, CompoundCode, RefinementMTZfree, 
                           DimplePathToPDB, DimplePathToMTZ, CompoundSMILES
                    FROM mainTable 
                    WHERE CrystalName IN ({','.join(['?']*len(crystal_names))})
                """
                datasets = pd.read_sql(query, conn, params=crystal_names)
            except Exception as e:
                logging.warning(f"CompoundSMILES column not available in database: {e}")
                query = f"""
                    SELECT CrystalName, CompoundCode, RefinementMTZfree, 
                           DimplePathToPDB, DimplePathToMTZ, NULL as CompoundSMILES
                    FROM mainTable 
                    WHERE CrystalName IN ({','.join(['?']*len(crystal_names))})
                """
                datasets = pd.read_sql(query, conn, params=crystal_names)
            # Do not overwrite the original CSV
        else:
            # Try to include CompoundSMILES, but don't fail if column doesn't exist
            try:
                query = """
                    SELECT CrystalName, CompoundCode, RefinementMTZfree, 
                           DimplePathToPDB, DimplePathToMTZ, CompoundSMILES
                    FROM mainTable 
                    WHERE RefinementOutcome = '1 - Analysis Pending'
                """
                datasets = pd.read_sql(query, conn)
            except Exception as e:
                logging.warning(f"CompoundSMILES column not available in database: {e}")
                query = """
                    SELECT CrystalName, CompoundCode, RefinementMTZfree, 
                           DimplePathToPDB, DimplePathToMTZ, NULL as CompoundSMILES
                    FROM mainTable 
                    WHERE RefinementOutcome = '1 - Analysis Pending'
                """
                datasets = pd.read_sql(query, conn)
            # Save to CSV for array job compatibility
            csv_path = os.path.join(params['Processing_directory'], 'datasets_metadata.csv')
            datasets.to_csv(csv_path, index=False)
            params['Filtered_dataset_csv_path'] = csv_path

        conn.close()

        # Check for missing or empty required fields (RefinementCIF is no longer required for SMILES-based workflow)
        required_fields = ['CompoundCode', 'RefinementMTZfree', 'DimplePathToPDB', 'DimplePathToMTZ']
        missing_info = []
        for idx, row in datasets.iterrows():
            missing = [field for field in required_fields if pd.isna(row[field]) or str(row[field]).strip() == '']
            if missing:
                missing_info.append({
                    'CrystalName': row.get('CrystalName', f'row {idx}'),
                    'MissingFields': missing
                })
        # Do not raise error, just log summary
        if missing_info:
            msg = ("One or more datasets are missing required fields. "
                   "These will be processed, but may fail downstream.\n" +
                   "Details:\n" +
                   '\n'.join([f"CrystalName: {item['CrystalName']} | Missing: {', '.join(item['MissingFields'])}" for item in missing_info]))
            logging.warning(msg)

        logging.info(f"{len(datasets)} datasets retrieved successfully.")
        print(f"Retrieved {len(datasets)} datasets from database.")
        if missing_info:
            print(f"Warning: {len(missing_info)} datasets have missing required fields and may fail processing.")
        return datasets

    except Exception as e:
        logging.error(f"Error retrieving datasets: {e}")
        raise

def fetch_password() -> str:
    """Securely fetch password from user input."""
    if not sys.stdin.isatty():
        logging.error("Standard input is not a TTY. Password input may not be secure.")
        raise RuntimeError("Insecure environment for password input.")
    return getpass.getpass("Enter your password: ")

def setup_dataset_directory(crystal_name: str, parent_dir: str) -> Tuple[str, str]:
    """Create dataset directory and input subdirectory."""
    dataset_dir = f"{parent_dir}/{crystal_name}"
    input_dir = f"{dataset_dir}/input_files"
    os.makedirs(input_dir, exist_ok=True)
    return dataset_dir, input_dir


def copy_input_files_and_prepare_for_restraints(
    input_dir: str, 
    refinement_mtzfree_path: str, 
    dimple_path_pdb: str, 
    dimple_path_mtz: str,
    compound_code: str, 
    smiles_string: Optional[str]
) -> Tuple[Optional[str], Optional[str]]:
    """
    Copy required input files to the dataset input directory. 
    Restraint generation will happen on cluster.
    """
    # Copy the main input files (MTZ, PDB, MTZ reference)
    files_to_copy = [
        (refinement_mtzfree_path, os.path.join(input_dir, os.path.basename(refinement_mtzfree_path))),
        (dimple_path_pdb, os.path.join(input_dir, os.path.basename(dimple_path_pdb))),
        (dimple_path_mtz, os.path.join(input_dir, os.path.basename(dimple_path_mtz)))
    ]
    
    for src, dst in files_to_copy:
        if not safe_file_copy(src, dst, compound_code):
            logging.error(f"Failed to copy required file for {compound_code}: {src}")
            return None, None
    
    # Return expected restraint file paths (will be generated on cluster)
    if smiles_string:
        logging.info(f"SMILES available for {compound_code} - restraints will be generated on cluster")
        expected_cif = os.path.join(input_dir, f"{compound_code}.cif")
        expected_pdb = os.path.join(input_dir, f"{compound_code}.pdb")
        return expected_cif, expected_pdb
    else:
        logging.warning(f"No SMILES available for {compound_code}, cannot generate restraints")
        return None, None


def process_pdb_file(input_dir: str, dimple_path_pdb: str, crystal_name: str, params: Dict[str, Any]) -> None:
    """Remove crystallisation components from PDB file if requested."""
    if params.get("Remove_crystallisation_components", False):
        pdb_path = os.path.join(input_dir, os.path.basename(dimple_path_pdb))
        if pdb_path and os.path.exists(pdb_path):
            with open(pdb_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()
            with open(pdb_path, 'w', encoding='utf-8') as f:
                for line in lines:
                    if not any(res in line for res in ["DMS", "EDO", "GOL", "SO4", "PO4", "PEG"]):
                        f.write(line)
        else:
            logging.warning(f"PDB path missing or does not exist for {crystal_name}: {pdb_path}")


def submit_jobs(
    datasets: pd.DataFrame, 
    params: Dict[str, Any], 
    output_json_path: Optional[str] = None, 
    output_csv_path: Optional[str] = None, 
    dry_run: bool = False
) -> Tuple[str, int, str]:
    """
    Prepares input files and directories for each dataset. 
    All job submission is handled via a single SLURM array job script.
    Assigns a unique Pipedream_N output directory for each dataset at preparation time.
    
    Returns:
        Tuple of (parent_dir, num_datasets, json_output_path)
    """
    try:
        print("Preparing input files and directories for datasets...")
        # Use single timestamp for all operations to ensure consistency
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

        print(f"Processing {len(datasets)} datasets...")
        for index, row in datasets.iterrows():
            crystal_name, compound_code, refinement_mtzfree, dimple_path_pdb, dimple_path_mtz, compound_smiles = extract_dataset_fields(row)

            required_paths = {
                "RefinementMTZfree": refinement_mtzfree,
                "DimplePathToPDB": dimple_path_pdb,
                "DimplePathToMTZ": dimple_path_mtz
            }
            missing_keys = [k for k, v in required_paths.items() if v is None or str(v).strip() == ""]
            if missing_keys:
                logging.error(f"Skipping dataset {crystal_name} due to missing required file paths: {missing_keys} | Values: {required_paths}")
                continue

            # Check if SMILES is available (we need it for grade2)
            if not compound_smiles or str(compound_smiles).strip() == '' or str(compound_smiles).strip().lower() in ['none', 'null']:
                logging.error(f"Skipping dataset {crystal_name} - no SMILES available for ligand restraint generation")
                continue

            filtered_datasets.append(row)

        # Prepare input files and output_yaml for filtered datasets only
        print(f"Generating restraints and copying input files for {len(filtered_datasets)} valid datasets...")
        successful_count = 0
        for row in filtered_datasets:
            crystal_name, compound_code, refinement_mtzfree, dimple_path_pdb, dimple_path_mtz, compound_smiles = extract_dataset_fields(row)
            
            logging.info(f"Processing dataset {crystal_name} (compound: {compound_code})")

            refinement_mtzfree_path = resolve_path(params['Processing_directory'], "analysis", "model_building", crystal_name, refinement_mtzfree)
            dimple_path_pdb = resolve_dataset_path(dimple_path_pdb, crystal_name, params['Processing_directory']) if dimple_path_pdb else None
            dimple_path_mtz = resolve_dataset_path(dimple_path_mtz, crystal_name, params['Processing_directory']) if dimple_path_mtz else None

            dataset_dir, input_dir = setup_dataset_directory(crystal_name, parent_dir)

            # Get SMILES for the compound (only call once in second loop)
            smiles_string = get_smiles_for_compound(
                compound_code, 
                compound_smiles,
                input_dir
            )

            # Find next available Pipedream_N output dir for this dataset
            # Use the same timestamp for consistency across all paths
            output_dir = f"{dataset_dir}/Pipedream_{timestamp}"

            try:
                # Check if all required files exist before proceeding
                file_paths = {
                    "RefinementMTZfree": refinement_mtzfree_path,
                    "DimplePathToPDB": dimple_path_pdb,
                    "DimplePathToMTZ": dimple_path_mtz
                }
                missing_files = validate_file_paths(file_paths, crystal_name)
                
                if missing_files:
                    logging.error(f"Skipping dataset {crystal_name} due to missing files: {'; '.join(missing_files)}")
                    continue
                
                # Prepare input files and get expected restraint paths (generation happens on cluster)
                generated_cif, generated_pdb = copy_input_files_and_prepare_for_restraints(
                    input_dir, refinement_mtzfree_path, dimple_path_pdb, dimple_path_mtz, 
                    compound_code, smiles_string
                )
                
                if not (generated_cif and generated_pdb):
                    logging.error(f"Cannot prepare restraints for {crystal_name}. Skipping.")
                    continue
                
                process_pdb_file(input_dir, dimple_path_pdb, crystal_name, params)

            except FileNotFoundError as e:
                logging.error(f"Missing file for {crystal_name}: {e}")
                logging.info(f"Skipping dataset {crystal_name} due to missing file: {e}")
                continue
            except TypeError as e:
                logging.error(f"TypeError for {crystal_name}: {e}. Variables: refinement_mtzfree_path={refinement_mtzfree_path}, dimple_path_pdb={dimple_path_pdb}, dimple_path_mtz={dimple_path_mtz}")
                logging.info(f"Skipping dataset {crystal_name} due to TypeError: {e}")
                continue

            # --- Build correct Pipedream command for SLURM script ---
            hklin_file = os.path.join(input_dir, os.path.basename(refinement_mtzfree_path))
            xyzin_file = os.path.join(input_dir, os.path.basename(dimple_path_pdb))
            hklref_file = os.path.join(input_dir, os.path.basename(dimple_path_mtz))
            # Use the generated restraint files
            rhofit_file = generated_cif
            rhocommands = params['Refinement_parameters'].get('rhocommands', None)
            pipedream_cmd = f"/dls_sw/apps/GPhL/BUSTER/20240123/scripts/pipedream -nolmr -hklin {hklin_file} -xyzin {xyzin_file} -hklref {hklref_file} -d {output_dir} "
            refinement_args = refinement_params_to_args(params.get('Refinement_parameters', {}), rhofit_file=rhofit_file, rhocommands=rhocommands)
            if refinement_args:
                pipedream_cmd += refinement_args + " "
            # Build output metadata for JSON (only fields needed by collate_pipedream_results.py)
            output_yaml[crystal_name] = {
                'Input_dir': input_dir,
                'CompoundCode': compound_code,
                'PipedreamDirectory': output_dir,
                'ReportHTML': f"{output_dir}/report-{compound_code}/index.html",
                'LigandReportHTML': f"{output_dir}/report-{compound_code}/ligand/index.html",
                'ExpectedSummary': f"{output_dir}/pipedream_summary.json",
                'PipedreamCommand': pipedream_cmd.strip(),
                'GeneratedCIF': generated_cif,
                'GeneratedPDB': generated_pdb,
                'InputSMILES': smiles_string
            }
            successful_count += 1
            logging.info(f"Successfully prepared dataset {crystal_name} ({successful_count} total)")

        # Log summary of processing results
        logging.info(f"Dataset processing complete: {successful_count}/{len(filtered_datasets)} datasets successfully prepared")
        print(f"Successfully prepared {successful_count}/{len(filtered_datasets)} datasets")

        # Save output_yaml to file as JSON
        if output_json_path:
            json_output_path = output_json_path
        else:
            output_json_name = f"Pipedream_{timestamp}_output.json"
            json_output_path = f"{parent_dir}/{output_json_name}"
        
        with open(json_output_path, 'w', encoding='utf-8') as json_file:
            json.dump(output_yaml, json_file, indent=2)
        logging.info(f"JSON metadata saved to: {json_output_path}")
        print(f"Metadata saved to: {json_output_path}")

        # Save filtered datasets to CSV for array job compatibility, including output_dir and full command
        # Only include datasets that were successfully processed (exist in output_yaml)
        successfully_processed_datasets = [row for row in filtered_datasets if row['CrystalName'] in output_yaml]
        
        if not successfully_processed_datasets:
            logging.warning("No datasets were successfully processed. No CSV will be created.")
            if dry_run:
                logging.info("DRY RUN: No datasets processed successfully.")
                print("DRY RUN: No datasets processed successfully.")
                return parent_dir, 0, json_output_path
            else:
                raise ValueError("No datasets were successfully processed.")
        
        filtered_df = pd.DataFrame(successfully_processed_datasets)
        # Extract output directories for successfully processed datasets only
        successful_output_dirs = [output_yaml[row['CrystalName']]['PipedreamDirectory'] for row in successfully_processed_datasets]
        filtered_df['PipedreamDirectory'] = successful_output_dirs
        # Add the full command for each dataset
        filtered_df['PipedreamCommand'] = [output_yaml[row['CrystalName']]['PipedreamCommand'] for row in successfully_processed_datasets]
        # Only keep the required columns for the SLURM script
        required_cols = ['CrystalName', 'CompoundCode', 'RefinementMTZfree', 'DimplePathToPDB', 'DimplePathToMTZ', 'CompoundSMILES', 'PipedreamDirectory', 'PipedreamCommand']
        filtered_df = filtered_df[required_cols]
        
        if output_csv_path:
            filtered_csv_path = output_csv_path
        else:
            filtered_csv_path = os.path.join(parent_dir, 'datasets_metadata.csv')
        
        filtered_df.to_csv(filtered_csv_path, index=False)
        params['Filtered_dataset_csv_path'] = filtered_csv_path
        logging.info(f"CSV datasets saved to: {filtered_csv_path}")

        # Print the number of valid datasets for SLURM
        logging.info(f"Number of valid datasets for SLURM: {filtered_df.shape[0]}")
        print(f"Prepared {filtered_df.shape[0]} datasets for processing.")
        
        if dry_run:
            logging.info("DRY RUN: Files prepared but no SLURM job will be submitted.")
            print("DRY RUN: Files prepared but no SLURM job will be submitted.")
            return parent_dir, filtered_df.shape[0], json_output_path

        # Write SLURM array index to CrystalName mapping for debugging in array_logs dir
        array_logs_dir = os.path.join(parent_dir, 'array_logs')
        os.makedirs(array_logs_dir, exist_ok=True)
        # Use the same timestamp for consistency across all paths
        slurm_map_path = os.path.join(array_logs_dir, f'slurm_array_index_map_{timestamp}.csv')
        with open(slurm_map_path, 'w', encoding='utf-8') as f:
            f.write('SLURM_ARRAY_TASK_ID,CrystalName\n')
            for idx, row in enumerate(successfully_processed_datasets):
                f.write(f'{idx},{row["CrystalName"]}\n')

        return parent_dir, filtered_df.shape[0], json_output_path

    except Exception as e:
        logging.error(f"Error preparing input files: {e}")
        raise

def refinement_params_to_args(
    refinement_params: Dict[str, Any], 
    rhofit_file: Optional[str] = None, 
    rhocommands: Optional[Any] = None
) -> str:
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

def write_array_job_script(
    datasets_csv: str, 
    output_dir: str, 
    num_datasets: int, 
    script_path: str, 
    processing_dir: str, 
    refinement_params: Optional[Dict[str, Any]] = None
) -> None:
    """Write SLURM array job script for batch processing."""
    # Only write array line if there is at least one dataset
    if num_datasets < 1:
        raise ValueError("No datasets to process: cannot write SLURM array job script.")
    
    # Cap at maximum concurrent jobs
    script_content = """#!/bin/bash
#SBATCH --job-name=pipedream_array
#SBATCH --partition=cs05r
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4096
#SBATCH --array=0-{num_datasets}%{max_jobs}
#SBATCH --output=array_logs/pipedream_array_%A_%a.out

# Create array_logs directory for SLURM output
mkdir -p array_logs

# Load modules (silence module system errors)
source /etc/profile 2>/dev/null
module load buster >/dev/null 2>&1
module load graphviz >/dev/null 2>&1
export LD_LIBRARY_PATH=/dls_sw/apps/graphviz/rhel8/12.2/lib:$LD_LIBRARY_PATH

# Set up CSD environment for grade2
export CSDHOME=/dls_sw/apps/CSDS/2024.1.0/
export BDG_TOOL_MOGUL=/dls_sw/apps/CSDS/2024.1.0/ccdc-software/mogul/bin/mogul

# Path to your datasets CSV
DATASETS_CSV="{datasets_csv}"
# Path to processing directory for input files
PROCESSING_DIR="{processing_dir}"

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
IFS=',' read -r CrystalName CompoundCode RefinementMTZfree DimplePathToPDB DimplePathToMTZ CompoundSMILES PipedreamDirectory PipedreamCommand <<< "$DATASET_LINE"

# Check for invalid or missing CrystalName
if [ -z "$CrystalName" ] || [[ "$CrystalName" =~ ^[0-9]+$ ]]; then
  echo "[ERROR] Invalid or missing CrystalName for SLURM_ARRAY_TASK_ID $TASK_ID. Exiting."
  exit 4
fi

# Print CrystalName at the top of the log for easier debugging
echo "[SLURM] CrystalName: $CrystalName (SLURM_ARRAY_TASK_ID: $TASK_ID)"

# Set up input/output directories as in your Python logic
OUTPUT_DIR_BASE="{output_dir}"
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

# If output dir exists, remove it to allow reprocessing
if [ -d "$OUTPUT_DIR" ]; then
  echo "[INFO] Output directory $OUTPUT_DIR already exists for $CrystalName. Removing for reprocessing..."
  rm -rf "$OUTPUT_DIR"
fi

# Copy input files from processing directory (not output directory)
cp "$PROCESSING_DIR/analysis/model_building/$CrystalName/$RefinementMTZfree" "$INPUT_DIR/"
cp "$DimplePathToPDB" "$INPUT_DIR/"
cp "$DimplePathToMTZ" "$INPUT_DIR/"

# Generate ligand restraints using grade2 if SMILES is available
if [ -n "$CompoundSMILES" ] && [ "$CompoundSMILES" != "None" ] && [ "$CompoundSMILES" != "NULL" ]; then
  echo "Generating ligand restraints for $CompoundCode using grade2..."
  echo "Using SMILES: $CompoundSMILES"
  
  # Check CSD availability for grade2
  echo "Checking CSD environment for grade2..."
  echo "BDG_TOOL_MOGUL: $BDG_TOOL_MOGUL"
  echo "BDG_TOOL_CSD_PYTHON_API: $BDG_TOOL_CSD_PYTHON_API"
  
  # Create SMILES file from database data
  SMILES_FILE="$INPUT_DIR/$CompoundCode.smiles"
  echo "$CompoundSMILES" > "$SMILES_FILE"
  
  # Create temporary SMILES file for grade2
  TEMP_SMILES="$INPUT_DIR/${{CompoundCode}}_temp.smi"
  echo "$CompoundSMILES" > "$TEMP_SMILES"
  
  # Run grade2 to generate restraints
  cd "$INPUT_DIR"
  echo "Running: grade2 -i $TEMP_SMILES -r LIG -o $CompoundCode -f"
  grade2 -i "$TEMP_SMILES" -r LIG -o "$CompoundCode" -f
  GRADE2_EXIT=$?
  
  # Clean up temporary file
  rm -f "$TEMP_SMILES"
  
  # Check if grade2 succeeded and rename files to expected names
  if [ $GRADE2_EXIT -eq 0 ] && [ -f "${{CompoundCode}}.restraints.cif" ] && [ -f "${{CompoundCode}}.xyz.pdb" ]; then
    echo "Successfully generated restraints for $CompoundCode"
    echo "Generated files: ${{CompoundCode}}.restraints.cif, ${{CompoundCode}}.xyz.pdb"
    
    # Rename files to match expected names for Pipedream
    mv "${{CompoundCode}}.restraints.cif" "${{CompoundCode}}.cif"
    mv "${{CompoundCode}}.xyz.pdb" "${{CompoundCode}}.pdb"
    
    echo "Renamed files to: ${{CompoundCode}}.cif, ${{CompoundCode}}.pdb"
  else
    echo "[ERROR] Failed to generate restraints for $CompoundCode. grade2 exit code: $GRADE2_EXIT"
    echo "Expected files from grade2: ${{CompoundCode}}.restraints.cif, ${{CompoundCode}}.xyz.pdb"
    echo "Note: If CSD is not available, consider using the Grade Webserver at http://grade.globalphasing.org"
    ls -la "$INPUT_DIR" | grep "$CompoundCode"
    exit 5
  fi
else
  echo "[ERROR] No SMILES available for $CompoundCode. CompoundSMILES value: '$CompoundSMILES'"
  exit 6
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
""".format(
        num_datasets=num_datasets-1,
        max_jobs=MAX_CONCURRENT_JOBS,
        datasets_csv=datasets_csv,
        processing_dir=processing_dir,
        output_dir=output_dir
    )
    
    # Write the script to file
    try:
        with open(script_path, "w", encoding='utf-8') as f:
            f.write(script_content)
        logging.info(f"SLURM script written to: {script_path} ({len(script_content)} characters)")
        
        # Verify the file was written correctly
        with open(script_path, "r", encoding='utf-8') as f:
            written_content = f.read()
        if len(written_content) == 0:
            raise ValueError(f"Script file is empty after writing: {script_path}")
        logging.info(f"SLURM script verification: file contains {len(written_content)} characters")
            
    except Exception as e:
        logging.error(f"Error writing SLURM script to {script_path}: {e}")
        raise

def submit_sbatch_on_wilson(script_path: str, processing_dir: str) -> None:
    """
    SSH to wilson and submit the job script using sbatch, using password authentication.
    Uses $CLUSTER_USER if set, otherwise current user.
    Uses paramiko for SSH and SFTP.
    """
    user = os.environ.get("CLUSTER_USER", os.getlogin())
    wilson_host = "wilson.diamond.ac.uk"
    remote_script_name = os.path.basename(script_path)
    remote_script_path = os.path.join(processing_dir, remote_script_name)
    print(f"Connecting to {wilson_host} as {user}...")
    password = getpass.getpass(f"Enter password for {user}@wilson: ")
    logging.info(f"[INFO] Connecting to wilson and copying job script...")
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
        logging.info(f"[INFO] Script copied to wilson:{remote_script_path}")
        print(f"Job script copied to wilson. Submitting job...")
        # Submit the job
        stdin, stdout, stderr = ssh.exec_command(f"cd {processing_dir} && sbatch {remote_script_name}")
        out = stdout.read().decode()
        err = stderr.read().decode()
        logging.info(f"[INFO] sbatch output:\n{out}")
        print(f"SLURM job submitted successfully!")
        if out.strip():
            print(f"Job submission output: {out.strip()}")
        # Only print sbatch error if it is not the known harmless module error
        harmless_module_error = "Module ERROR: no such variable"
        if err and harmless_module_error not in err:
            logging.info(f"[ERROR] sbatch error:\n{err}")
        else:
            if err:
                # Log the harmless error for reference
                logging.info(f"[sbatch harmless stderr suppressed from terminal]: {err.strip()}")
            logging.info("[INFO] Job script submitted on wilson via sbatch.")
        ssh.close()
    except Exception as e:
        logging.info(f"[ERROR] Failed to submit job on wilson: {e}")
        print(f"ERROR: Failed to submit job on wilson: {e}")
        logging.info(f"You can submit manually with:\n  scp {script_path} {user}@wilson:{remote_script_path}\n  ssh {user}@wilson 'cd {processing_dir} && sbatch {remote_script_name}'")

def resolve_dataset_path(path: str, crystal_name: str, processing_dir: str) -> str:
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
        parameters_file = args.parameters
        
        # Set up logging early with user preferences
        setup_logging(log_level=args.log_level, verbose=args.verbose)
        
        if args.verbose:
            print(f"Starting Pipedream XChem pipeline with parameters:")
            print(f"  Parameters file: {parameters_file}")
            print(f"  Log level: {args.log_level}")
            print(f"  Dry run: {args.dry_run}")
            if args.output_json:
                print(f"  Custom JSON output: {args.output_json}")
            if args.output_csv:
                print(f"  Custom CSV output: {args.output_csv}")
        
        print("Reading parameters and validating configuration...")
        params = read_yaml(parameters_file)
        validate_params(params)
        print("Getting datasets from database...")
        datasets = get_datasets(params)
        if datasets is not None and not datasets.empty:
            parent_dir, num_valid_datasets, json_output_path = submit_jobs(
                datasets, 
                params, 
                output_json_path=args.output_json,
                output_csv_path=args.output_csv,
                dry_run=args.dry_run
            )
            
            if not args.dry_run:
                # Write and submit array job script
                datasets_csv = params['Filtered_dataset_csv_path']
                output_dir = params.get('Output_directory') or params['Processing_directory']
                script_path = os.path.join(parent_dir, 'pipedream_array_job.sh')
                print("Writing SLURM array job script...")
                write_array_job_script(datasets_csv, output_dir, num_valid_datasets, script_path, params['Processing_directory'], params.get('Refinement_parameters'))
                logging.info(f"Submitting SLURM array job for {num_valid_datasets} datasets...")
                print(f"Submitting SLURM array job for {num_valid_datasets} datasets...")
                submit_sbatch_on_wilson(script_path, parent_dir)
            else:
                logging.info(f"DRY RUN complete. JSON metadata: {json_output_path}")
                logging.info(f"DRY RUN complete. CSV datasets: {params['Filtered_dataset_csv_path']}")
                print(f"DRY RUN complete.")
                print(f"JSON metadata: {json_output_path}")
                print(f"CSV datasets: {params['Filtered_dataset_csv_path']}")
        else:
            logging.warning("No datasets to process. Exiting.")
            print("No datasets to process. Exiting.")
            
    except Exception as e:
        import traceback
        error_details = traceback.format_exc()
        logging.error(f"Script execution failed: {e}")
        logging.error(f"Full traceback: {error_details}")
        print(f"ERROR: Script execution failed: {e}")
        print(f"Check the log file for full error details.")
        sys.exit(1)

