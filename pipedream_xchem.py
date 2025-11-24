"""
Pipedream XChem Integration Script

This script submits SLURM jobs for Pipedream analysis with ligand restraints 
generated from SMILES using grade2. It processes crystallographic datasets 
from an SQLite database and prepares them for automated refinement.

Author: DFearon
Date: November 2025
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
from pathlib import Path
from contextlib import contextmanager
from concurrent.futures import ThreadPoolExecutor, as_completed

# Add tqdm with graceful fallback
try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False
    # Fallback: tqdm does nothing if not installed
    def tqdm(iterable, **kwargs):
        return iterable

# Configuration Constants
VERSION = "1.0.2"


# Exit codes documentation
EXIT_CODES = {
    0: "Success",
    1: "General error or user cancelled",
    2: "Empty OUTPUT_DIR",
    3: "Invalid dataset line in CSV",
    4: "Invalid CrystalName",
    5: "Failed to generate restraints with grade2",
    6: "No SMILES available"
}

# Cluster Configuration
CLUSTER_BASTION = "wilson.diamond.ac.uk"
CLUSTER_USER = os.environ.get("CLUSTER_USER", os.getlogin())
MAX_CONCURRENT_JOBS = 100
MAX_ARRAY_SIZE = 1000
DEFAULT_SBATCH_TIMEOUT = 30
DEFAULT_CPUS_PER_TASK = 4
DEFAULT_MEM_PER_CPU_MB = 4096
DEFAULT_CLUSTER_PARTITION = "cs05r"
DEFAULT_JOB_PRIORITY = "normal"
VALID_PARTITIONS = ["cs05r", "cs04r"]
VALID_PRIORITIES = ["normal", "low", "high"]

# Software Paths (fallback if module load buster fails)
BUSTER_PIPEDREAM_PATH_FALLBACK = "/dls_sw/apps/GPhL/BUSTER/20250717/scripts/pipedream"
BUSTER_REPORT_PATH_FALLBACK = "/dls_sw/apps/GPhL/BUSTER/20240123/scripts/buster-report"
CSDHOME_PATH = "/dls_sw/apps/CSDS/2024.1.0/"
MOGUL_PATH = "/dls_sw/apps/CSDS/2024.1.0/ccdc-software/mogul/bin/mogul"
GRAPHVIZ_LIB_PATH = "/dls_sw/apps/graphviz/rhel8/12.2/lib"


def setup_logging(log_dir: str = None, log_level: str = "INFO", verbose: bool = False, log_file: str = 'pipedream.log') -> None:
    """Configure logging with optional custom log file name."""
    if log_dir:
        os.makedirs(log_dir, exist_ok=True)
        log_path = os.path.join(log_dir, log_file)
    else:
        log_path = log_file
    
    # Convert string level to logging constant
    numeric_level = getattr(logging, log_level.upper(), logging.INFO)
    
    # Clear any existing handlers
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    
    handlers = [logging.FileHandler(log_path)]
    
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
    """
    try:
        file_size = os.path.getsize(src)
        if file_size > 100 * 1024 * 1024:  # Log for files > 100MB
            logging.info(f"Copying large file ({file_size/1024/1024:.1f} MB): {src}")
        
        shutil.copy2(src, dst)
        logging.debug(f"Copied {src} -> {dst} for crystal {crystal_name}")
        return True
    except FileNotFoundError:
        logging.error(f"Source file not found: {src} for crystal {crystal_name}")
        return False
    except PermissionError:
        logging.error(f"Permission denied copying {src} for crystal {crystal_name}")
        return False
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

def resolve_path(base_dir: str, *parts: str, must_exist: bool = False) -> Optional[str]:
    """
    Unified path resolution function.
    
    Handles both:
    - Absolute paths (uses as-is)
    - Relative paths (joins with base_dir)
    - Multiple path components (joins all parts)
    
    Args:
        base_dir: Base directory for relative paths
        *parts: Path components to join (can be relative or absolute)
        must_exist: If True, return None if path doesn't exist
        
    Returns:
        Absolute path or None if validation fails
        
    Examples:
        resolve_path("/base", "subdir", "file.txt")
        resolve_path("/base", "/absolute/path/file.txt")  # Returns absolute path
        resolve_path("/base", "relative/path", must_exist=True)
    """
    if not parts:
        logging.debug("resolve_path called with no path parts")
        return None
    
    # Join all parts
    path_str = os.path.join(*parts) if len(parts) > 1 else parts[0]
    
    # If already absolute, use as-is; otherwise join with base_dir
    if os.path.isabs(path_str):
        path = os.path.realpath(path_str)
    else:
        path = os.path.realpath(os.path.join(base_dir, path_str))
    
    # Validate existence if required
    if must_exist and not os.path.exists(path):
        logging.debug(f"Path does not exist: {path}")
        return None
    
    return path

def get_smiles_for_compound(
    compound_code: str, 
    smiles_from_db: Optional[str], 
    input_dir: str
) -> Optional[str]:
    """Get SMILES string for a compound and save to file."""
    if not validate_smiles(smiles_from_db, compound_code):
        return None
    
    smiles_string = str(smiles_from_db).strip()
    logging.debug(f"Using SMILES from database for {compound_code}: {smiles_string}")
    
    # Always save SMILES to file (required for downstream processing)
    os.makedirs(input_dir, exist_ok=True)
    smiles_file_path = os.path.join(input_dir, f"{compound_code}.smiles")
    try:
        with open(smiles_file_path, 'w', encoding='utf-8') as f:
            f.write(f"{smiles_string}\n")
        logging.debug(f"Saved SMILES to file: {smiles_file_path}")
    except Exception as e:
        logging.error(f"Failed to save SMILES file for {compound_code}: {e}")
        return None  # Return None if we can't save (critical for processing)
    
    return smiles_string

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

def validate_yaml_paths(params: Dict[str, Any]) -> bool:
    """
    Check if YAML contains placeholder paths that need to be replaced.
    
    Args:
        params: Configuration parameters dictionary
        
    Returns:
        True if no placeholders found, False otherwise
        
    Common placeholders:
        - <proposal>
        - <visit>
        - <Pipedream_run>
    """
    placeholders = ['<proposal>', '<visit>', '<Pipedream_run>']
    found_placeholders = []
    
    def check_value(key: str, value: Any, path: str = "") -> None:
        """Recursively check for placeholders in nested structures."""
        full_key = f"{path}.{key}" if path else key
        
        if isinstance(value, str):
            for placeholder in placeholders:
                if placeholder in value:
                    found_placeholders.append({
                        'key': full_key,
                        'value': value,
                        'placeholder': placeholder
                    })
        elif isinstance(value, dict):
            for k, v in value.items():
                check_value(k, v, full_key)
        elif isinstance(value, list):
            for i, item in enumerate(value):
                check_value(f"{key}[{i}]", item, path)
    
    # Check all parameters
    for key, value in params.items():
        check_value(key, value)
    
    # Report findings
    if found_placeholders:
        logging.error("Configuration contains placeholder paths that must be replaced:")
        print("\n" + "="*70)
        print("ERROR: YAML Configuration Contains Placeholder Paths")
        print("="*70)
        for item in found_placeholders:
            logging.error(f"  {item['key']}: {item['value']}")
            print(f"  Parameter: {item['key']}")
            print(f"  Value: {item['value']}")
            print(f"  Found placeholder: {item['placeholder']}")
            print()
        print("Please replace these placeholders with actual paths before running.")
        print("="*70 + "\n")
        return False
    
    return True


def validate_params(params: Dict[str, Any]) -> None:
    """Validate required parameters in the configuration."""
    required_keys = ['Database_path', 'Mode', 'Processing_directory', 'Refinement_parameters']
    refinement_keys = ['keepwater', 'TLS', 'remediate', 'sidechainrebuild', 'runpepflip', 'rhocommands']

    for key in required_keys:
        if key not in params:
            raise ValueError(f"Missing required key in YAML parameters: {key}")

    refinement_params = params.get('Refinement_parameters', {})
    for key in refinement_keys:
        if key not in refinement_params:
            raise ValueError(f"Missing required refinement parameter: {key}")

    if params['Mode'] == 'specific_datasets' and 'Dataset_csv_path' not in params:
        raise ValueError("Missing 'Dataset_csv_path' for 'specific_datasets' mode")
    
    # Validate optional cluster configuration
    if 'Cluster_partition' in params:
        partition = params['Cluster_partition']
        if partition not in VALID_PARTITIONS:
            raise ValueError(f"Invalid Cluster_partition '{partition}'. Must be one of: {', '.join(VALID_PARTITIONS)}")
    
    if 'Job_priority' in params:
        priority = params['Job_priority']
        if priority not in VALID_PRIORITIES:
            raise ValueError(f"Invalid Job_priority '{priority}'. Must be one of: {', '.join(VALID_PRIORITIES)}")
    
    # Check for placeholder paths
    if not validate_yaml_paths(params):
        raise ValueError("Configuration contains unresolved placeholder paths. Please update your YAML file with actual paths.")

def get_datasets(params: Dict[str, Any]) -> pd.DataFrame:
    """Retrieve datasets from the database based on configuration parameters."""
    try:
        conn = sqlite3.connect(params['Database_path'])

        # Define base query components
        base_columns = """
            SELECT CrystalName, CompoundCode, RefinementMTZfree, 
                   DimplePathToPDB, DimplePathToMTZ, CompoundSMILES
            FROM mainTable
        """
        
        # Build WHERE clause based on mode
        if params['Mode'] == 'specific_datasets':
            csv_path = params.get('Dataset_csv_path') or params.get('Filtered_dataset_csv_path')
            if not csv_path:
                raise ValueError("Missing 'Dataset_csv_path' for 'specific_datasets' mode")
            csv_df = pd.read_csv(csv_path)
            if 'CrystalName' not in csv_df.columns:
                raise ValueError("CSV must contain a 'CrystalName' column.")
            crystal_names = csv_df['CrystalName'].dropna().unique().tolist()
            if not crystal_names:
                raise ValueError("No CrystalName values found in the CSV.")
            
            where_clause = f"WHERE CrystalName IN ({','.join(['?']*len(crystal_names))})"
            query_params = crystal_names
        else:
            where_clause = "WHERE RefinementOutcome = '1 - Analysis Pending'"
            query_params = None
        
        # Try to get data with SMILES, fallback without
        try:
            query = base_columns + where_clause
            datasets = pd.read_sql(query, conn, params=query_params) if query_params else pd.read_sql(query, conn)
        except Exception as e:
            logging.warning(f"CompoundSMILES column not available: {e}")
            query_no_smiles = base_columns.replace(", CompoundSMILES", ", NULL as CompoundSMILES") + where_clause
            datasets = pd.read_sql(query_no_smiles, conn, params=query_params) if query_params else pd.read_sql(query_no_smiles, conn)
        
        conn.close()
        
        # Save CSV if not in specific_datasets mode
        if params['Mode'] != 'specific_datasets':
            csv_path = os.path.join(params['Processing_directory'], 'datasets_metadata.csv')
            datasets.to_csv(csv_path, index=False)
            params['Filtered_dataset_csv_path'] = csv_path
        
        # Check for missing or empty required fields
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
) -> bool:
    """
    Copy required input files to the dataset input directory.
    Restraint generation will happen on cluster.
    
    Returns:
        True if files copied successfully and SMILES is available, False otherwise
    """
    logging.debug(f"Copying input files for {compound_code}:")
    logging.debug(f"RefinementMTZfree: {refinement_mtzfree_path} -> {input_dir}/")
    logging.debug(f"DimplePathToPDB: {dimple_path_pdb} -> {input_dir}/")
    logging.debug(f"DimplePathToMTZ: {dimple_path_mtz} -> {input_dir}/")
    
    # Copy the main input files
    files_to_copy = [
        (refinement_mtzfree_path, os.path.join(input_dir, os.path.basename(refinement_mtzfree_path))),
        (dimple_path_pdb, os.path.join(input_dir, os.path.basename(dimple_path_pdb))),
        (dimple_path_mtz, os.path.join(input_dir, os.path.basename(dimple_path_mtz)))
    ]
    
    for src, dst in files_to_copy:
        if not safe_file_copy(src, dst, compound_code):
            logging.error(f"Failed to copy required file for {compound_code}: {src}")
            return False
    
    logging.debug(f"Successfully copied all input files for {compound_code}")
    
    # Check if SMILES is available for restraint generation
    if smiles_string:
        logging.debug(f"SMILES available for {compound_code} - restraints will be generated on cluster")
        return True
    else:
        logging.warning(f"No SMILES available for {compound_code}, cannot generate restraints")
        return False


def process_pdb_file(input_dir: str, dimple_path_pdb: str, crystal_name: str, params: Dict[str, Any]) -> None:
    """Remove crystallisation components from PDB file if requested."""
    if params.get("Remove_crystallisation_components", False):
        pdb_path = os.path.join(input_dir, os.path.basename(dimple_path_pdb))
        logging.debug(f"Remove_crystallisation_components enabled for {crystal_name}")
        logging.debug(f"Processing PDB file: {pdb_path}")
        
        if pdb_path and os.path.exists(pdb_path):
            with open(pdb_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()
            
            # Count removals by component type
            original_count = len(lines)
            components_to_remove = ["DMS", "EDO", "GOL", "SO4", "PO4", "PEG"]
            removed_counts = {comp: 0 for comp in components_to_remove}
            
            kept_lines = []
            for line in lines:
                if any(res in line for res in components_to_remove):
                    # Count which component was found
                    for comp in components_to_remove:
                        if comp in line:
                            removed_counts[comp] += 1
                            break
                else:
                    kept_lines.append(line)
            
            # Write cleaned file
            with open(pdb_path, 'w', encoding='utf-8') as f:
                f.writelines(kept_lines)
            
            removed_total = original_count - len(kept_lines)
            if removed_total > 0:
                component_summary = ", ".join([f"{comp}: {count}" for comp, count in removed_counts.items() if count > 0])
                logging.debug(f"Removed {removed_total} lines from {crystal_name} ({component_summary})")
            else:
                logging.debug(f"No crystallisation components found in {crystal_name}")
        else:
            logging.warning(f"PDB path missing or does not exist for {crystal_name}: {pdb_path}")


def submit_jobs(
    datasets: pd.DataFrame, 
    params: Dict[str, Any], 
    output_json_path: Optional[str] = None, 
    output_csv_path: Optional[str] = None, 
    dry_run: bool = False,
    log_level: str = "INFO",
    verbose: bool = False
) -> Tuple[str, int, str]:
    """Prepares input files and directories for each dataset."""
    try:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        parent_dir = params.get('Output_directory') or f"{params['Processing_directory']}/analysis/Pipedream/Pipedream_{timestamp}"
        os.makedirs(parent_dir, exist_ok=True)

        # setup_logging(log_dir=parent_dir, log_level=log_level, verbose=verbose)

        logging.info(f"Processing {len(datasets)} datasets...")
        logging.info(f"Output directory: {parent_dir}")
        
        output_yaml = {}
        filtered_datasets = []

        # Progress bar setup
        if HAS_TQDM and not verbose:
            # Quiet mode: show progress bar
            print(f"Processing {len(datasets)} datasets...")
            iterator = tqdm(
                datasets.iterrows(), 
                total=len(datasets),
                desc="Processing datasets",
                unit="dataset",
                ncols=80,
                bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]'
            )
        else:
            # Verbose mode or no tqdm: show detailed output
            print(f"Processing {len(datasets)} datasets...")
            iterator = datasets.iterrows()
        
        # OPTION 1: Keep sequential (current, safe)
        for index, row in iterator:
            crystal_name, compound_code, refinement_mtzfree, dimple_path_pdb, dimple_path_mtz, compound_smiles = extract_dataset_fields(row)

            # Validate required paths
            required_paths = {
                "RefinementMTZfree": refinement_mtzfree,
                "DimplePathToPDB": dimple_path_pdb,
                "DimplePathToMTZ": dimple_path_mtz
            }
            missing_keys = [k for k, v in required_paths.items() if v is None or str(v).strip() == ""]
            if missing_keys:
                logging.error(f"Skipping {crystal_name} - missing paths: {missing_keys}")
                continue

            # Check SMILES availability
            if not validate_smiles(compound_smiles, crystal_name):
                logging.error(f"Skipping {crystal_name} - no valid SMILES")
                continue

            # Process dataset
            try:
                metadata = process_single_dataset(row, params, parent_dir, timestamp)
                if metadata:
                    output_yaml[crystal_name] = metadata
                    filtered_datasets.append(row)
                    logging.debug(f"Successfully prepared {crystal_name} ({len(filtered_datasets)} total)")
                    
                    # Update progress bar description with success count
                    if HAS_TQDM and not verbose and hasattr(iterator, 'set_postfix'):
                        iterator.set_postfix({'✓': len(filtered_datasets), '✗': index + 1 - len(filtered_datasets)})
            except Exception as e:
                logging.error(f"Error processing {crystal_name}: {e}")
                continue

        # Summary
        logging.debug(f"Processing complete: {len(filtered_datasets)}/{len(datasets)} datasets prepared")
        print(f"\n✓ Successfully prepared {len(filtered_datasets)}/{len(datasets)} datasets.")

        # Save JSON metadata
        json_output_path = output_json_path or f"{parent_dir}/Pipedream_{timestamp}_output.json"
        with open(json_output_path, 'w', encoding='utf-8') as json_file:
            json.dump(output_yaml, json_file, indent=2)
        logging.debug(f"JSON metadata saved to: {json_output_path}")
        print(f"Metadata saved to: {json_output_path}")

        # Handle empty results
        if not filtered_datasets:
            logging.warning("No datasets successfully processed.")
            if dry_run:
                print("Dry run: No datasets processed successfully.")
                return parent_dir, 0, json_output_path
            else:
                raise ValueError("No datasets were successfully processed.")
        
        # Save CSV with metadata
        filtered_df = pd.DataFrame(filtered_datasets)
        filtered_df['PipedreamDirectory'] = [output_yaml[row['CrystalName']]['PipedreamDirectory'] for row in filtered_datasets]
        filtered_df['PipedreamCommand'] = [output_yaml[row['CrystalName']]['PipedreamCommand'] for row in filtered_datasets]
        
        required_cols = ['CrystalName', 'CompoundCode', 'RefinementMTZfree', 'DimplePathToPDB', 
                        'DimplePathToMTZ', 'CompoundSMILES', 'PipedreamDirectory', 'PipedreamCommand']
        filtered_df = filtered_df[required_cols]
        
        filtered_csv_path = output_csv_path or os.path.join(parent_dir, 'datasets_metadata.csv')
        filtered_df.to_csv(filtered_csv_path, index=False)
        params['Filtered_dataset_csv_path'] = filtered_csv_path
        logging.debug(f"CSV datasets saved to: {filtered_csv_path}")
        
        if dry_run:
            logging.debug("Dry run: Files prepared but no SLURM jobs will be submitted.")
            print("Dry run: Files prepared but no SLURM jobs will be submitted.")
            return parent_dir, len(filtered_df), json_output_path

        # Create SLURM array index mapping
        array_logs_dir = os.path.join(parent_dir, 'array_logs')
        os.makedirs(array_logs_dir, exist_ok=True)
        slurm_map_path = os.path.join(array_logs_dir, f'slurm_array_index_map_{timestamp}.csv')
        with open(slurm_map_path, 'w', encoding='utf-8') as f:
            f.write('SLURM_ARRAY_TASK_ID,CrystalName\n')
            for idx, row in enumerate(filtered_datasets):
                f.write(f'{idx},{row["CrystalName"]}\n')

        return parent_dir, len(filtered_df), json_output_path

    except Exception as e:
        logging.error(f"Error preparing input files: {e}")
        raise

def validate_smiles(smiles: Optional[str], compound_code: str = "Unknown") -> bool:
    """
    Validate SMILES string (consolidated validation logic).
    
    Args:
        smiles: SMILES string to validate
        compound_code: Compound code for logging context (optional)
        
    Returns:
        True if valid SMILES string, False otherwise
        
    Examples:
        >>> validate_smiles("CCO")
        True
        >>> validate_smiles(None)
        False
        >>> validate_smiles("none")
        False
    """
    if not smiles or str(smiles).strip().lower() in ['none', 'null', 'nan', '']:
        return False
    
    smiles_str = str(smiles).strip().lower()
    
    # Check for null placeholders
    invalid_values = ['none', 'null', 'nan', '']
    is_valid = smiles_str not in invalid_values
    
    if not is_valid and compound_code != "Unknown":
        logging.debug(f"Invalid SMILES placeholder for {compound_code}: {smiles_str}")
    
    return is_valid

def refinement_params_to_args(
    refinement_params: Dict[str, Any], 
    rhofit_file: Optional[str] = None, 
    rhocommands: Optional[Any] = None
) -> str:
    """
    Convert Refinement_parameters dict to Pipedream command-line arguments.
    
    Special handling:
    - TLS and WaterUpdatePkmaps are combined into single -mrefine flag
    - rhocommands are properly formatted with dashes
    
    Returns:
        Space-separated command line arguments
    """
    args = []
    
    # Build -mrefine flag (combines TLS and WaterUpdatePkmaps)
    tls_value = refinement_params.get('TLS', '')
    water_update_pkmaps = refinement_params.get('WaterUpdatePkmaps', False)
    
    if tls_value or water_update_pkmaps:
        mrefine_parts = []
        if tls_value:
            mrefine_parts.append(str(tls_value).strip())
        if water_update_pkmaps:
            mrefine_parts.append('WaterUpdatePkmaps')
        if mrefine_parts:
            args.append(f"-mrefine {','.join(mrefine_parts)}")
    
    # Process other parameters (skip special cases)
    skip_keys = {'rhocommands', 'WaterUpdatePkmaps', 'TLS'}
    
    for key, value in refinement_params.items():
        if key in skip_keys:
            continue
            
        if isinstance(value, bool):
            if value:
                args.append(f"-{key}")
        elif isinstance(value, str) and value.strip():
            args.append(f"-{key} {value}")
        elif isinstance(value, list):
            joined = ' '.join(str(v) for v in value if str(v).strip())
            if joined:
                args.append(f"-{key} {joined}")
    
    # Add rhofit file
    if rhofit_file:
        args.append(f"-rhofit {rhofit_file}")
    
    # Add rhocommands with proper formatting
    if rhocommands:
        def ensure_dash(cmd):
            cmd = str(cmd).strip()
            return cmd if cmd.startswith('-') else f'-{cmd}' if cmd else ''
        
        if isinstance(rhocommands, list):
            rhocmds = [ensure_dash(cmd) for cmd in rhocommands if str(cmd).strip()]
            rhocmd_str = ' '.join(rhocmds)
        else:
            rhocmd_str = ensure_dash(rhocommands)
        
        if rhocmd_str:
            # Clean up formatting
            rhocmd_str = rhocmd_str.strip('"\'')
            rhocmd_str = ' '.join(rhocmd_str.split())
            args.append(f'-rhocommands {rhocmd_str}')
    
    return " ".join(args)

def generate_slurm_header(num_datasets: int, array_offset: int = 0, partition: str = DEFAULT_CLUSTER_PARTITION, priority: str = DEFAULT_JOB_PRIORITY) -> str:
    """
    Generate SLURM job script header with directives.
    
    Args:
        num_datasets: Number of datasets to process
        array_offset: Offset for chunked arrays (default 0)
        partition: Cluster partition to use (default: cs05r)
        priority: Job priority level (default: normal)
        
    Returns:
        String with SLURM header directives
    """
    array_max = num_datasets - 1
    max_concurrent = min(MAX_CONCURRENT_JOBS, num_datasets)
    
    # Build priority directive based on priority level
    # Use --nice for priority (positive = lower priority, negative = higher priority)
    priority_directive = ""
    if priority == "low":
        priority_directive = "#SBATCH --nice=1000\n"  # Low priority (runs after others)
    elif priority == "high":
        priority_directive = "#SBATCH --nice=-100\n"  # High priority
    # normal priority doesn't need a directive
    
    return f"""#!/bin/bash
#SBATCH --job-name=pipedream_array
#SBATCH --partition={partition}
{priority_directive}#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4096
#SBATCH --array=0-{array_max}%{max_concurrent}
#SBATCH --output=array_logs/pipedream_array_%A_%a.out

# Exit codes: {', '.join(f'{k}={v}' for k, v in EXIT_CODES.items())}
"""


def generate_environment_setup() -> str:
    """
    Generate environment setup section for SLURM script.
    
    Returns:
        String with module loads and environment variables
    """
    return f"""
mkdir -p array_logs
source /etc/profile 2>/dev/null
module load buster >/dev/null 2>&1
module load graphviz >/dev/null 2>&1
export LD_LIBRARY_PATH={GRAPHVIZ_LIB_PATH}:$LD_LIBRARY_PATH
export CSDHOME={CSDHOME_PATH}
export BDG_TOOL_MOGUL={MOGUL_PATH}
"""


def generate_csv_parsing_section(datasets_csv: str, processing_dir: str, array_offset: int) -> str:
    """
    Generate CSV parsing and validation section.
    
    Args:
        datasets_csv: Path to CSV file with dataset metadata
        processing_dir: Base processing directory
        array_offset: Offset for chunked arrays
        
    Returns:
        String with bash code for CSV parsing
    """
    csv_line_offset = array_offset + 2  # +2 for header and 0-based indexing
    
    return f"""
DATASETS_CSV="{datasets_csv}"
PROCESSING_DIR="{processing_dir}"
TASK_ID=$SLURM_ARRAY_TASK_ID

# Calculate CSV line number
CSV_LINE_NUM=$(($TASK_ID + {csv_line_offset}))
DATASET_LINE=$(awk -v n=$CSV_LINE_NUM 'NR==n' "$DATASETS_CSV")

# Validate dataset line
if [ -z "$DATASET_LINE" ] || [[ "$DATASET_LINE" =~ ^[0-9]+$ ]]; then
  echo "[ERROR] No valid dataset at CSV line $CSV_LINE_NUM (array task $TASK_ID, offset {array_offset})"
  exit 3
fi

# Parse CSV line
IFS=',' read -r CrystalName CompoundCode RefinementMTZfree DimplePathToPDB DimplePathToMTZ CompoundSMILES PipedreamDirectory PipedreamCommand <<< "$DATASET_LINE"

# Validate crystal name
if [ -z "$CrystalName" ] || [[ "$CrystalName" =~ ^[0-9]+$ ]]; then
  echo "[ERROR] Invalid CrystalName at CSV line $CSV_LINE_NUM"
  exit 4
fi

echo "[SLURM] Processing: $CrystalName (array task $TASK_ID, CSV line $CSV_LINE_NUM)"
"""


def generate_directory_setup_section(output_dir: str) -> str:
    """
    Generate directory setup and logging redirection.
    
    Args:
        output_dir: Base output directory
        
    Returns:
        String with bash code for directory setup
    """
    return f"""
# Setup directories
OUTPUT_DIR_BASE="{output_dir}"
DATASET_DIR="$OUTPUT_DIR_BASE/$CrystalName"
INPUT_DIR="$DATASET_DIR/input_files"
OUTPUT_DIR="$PipedreamDirectory"

mkdir -p "$DATASET_DIR"
LOG_BASENAME="${{CrystalName}}_slurm_$(basename $OUTPUT_DIR).out"
exec > "$DATASET_DIR/$LOG_BASENAME" 2>&1

# Validate output directory
if [ -z "$OUTPUT_DIR" ]; then
  echo "[ERROR] OUTPUT_DIR is empty for $CrystalName"
  exit 2
fi

mkdir -p "$INPUT_DIR"

# Clean existing output
if [ -d "$OUTPUT_DIR" ]; then
  echo "[INFO] Removing existing output directory..."
  rm -rf "$OUTPUT_DIR"
fi
"""


def generate_file_copy_section(processing_dir: str) -> str:
    """
    Generate file copying section.
    
    Args:
        processing_dir: Base processing directory
        
    Returns:
        String with bash code for copying input files
    """
    return f"""
# Copy input files
cp "$PROCESSING_DIR/analysis/model_building/$CrystalName/$RefinementMTZfree" "$INPUT_DIR/"
cp "$DimplePathToPDB" "$INPUT_DIR/"
cp "$DimplePathToMTZ" "$INPUT_DIR/"
"""


def generate_restraint_generation_section() -> str:
    """
    Generate ligand restraint generation section using grade2.
    
    Note: SMILES validation happens in Python before job submission.
    This script assumes valid SMILES is available.
    
    Returns:
        String with bash code for grade2 restraint generation
    """
    return """
# Generate ligand restraints using SMILES
# Note: SMILES validation already performed before job submission
echo "Generating ligand restraints for $CompoundCode using grade2..."
echo "Using SMILES: $CompoundSMILES"

TEMP_SMILES="$INPUT_DIR/${CompoundCode}_temp.smi"
echo "$CompoundSMILES" > "$TEMP_SMILES"

cd "$INPUT_DIR"
echo "Running: grade2 -i $TEMP_SMILES -r LIG -o $CompoundCode -f"
grade2 -i "$TEMP_SMILES" -r LIG -o "$CompoundCode" -f
GRADE2_EXIT=$?

rm -f "$TEMP_SMILES"

# Validate grade2 output
if [ $GRADE2_EXIT -eq 0 ] && [ -f "${CompoundCode}.restraints.cif" ] && [ -f "${CompoundCode}.xyz.pdb" ]; then
  echo "Successfully generated restraints for $CompoundCode"
  mv "${CompoundCode}.restraints.cif" "${CompoundCode}.cif"
  mv "${CompoundCode}.xyz.pdb" "${CompoundCode}.pdb"
  echo "Renamed files to: ${CompoundCode}.cif, ${CompoundCode}.pdb"
else
  echo "[ERROR] Failed to generate restraints for $CompoundCode. grade2 exit code: $GRADE2_EXIT"
  ls -la "$INPUT_DIR" | grep "$CompoundCode"
  exit 5
fi
"""


def generate_pipedream_execution_section() -> str:
    """
    Generate Pipedream execution, post-processing, and cleanup section.
    
    Post-processing includes:
    - Map generation (2fofc and fofc) using gemmi
    - Edstats analysis for ligand RSR calculation
    
    Returns:
        String with bash code for running Pipedream and post-processing
    """
    return """
# Run Pipedream
CMD=$(echo "$PipedreamCommand" | sed 's/^"//;s/"$//' | tr -s ' ')
echo "Running command: $CMD"
bash -c "$CMD"
PIPEDREAM_EXIT=$?

if [ $PIPEDREAM_EXIT -ne 0 ]; then
  echo "[ERROR] Pipedream failed with exit code: $PIPEDREAM_EXIT"
  exit $PIPEDREAM_EXIT
fi

echo "[INFO] Pipedream completed successfully"

# Post-processing: Generate maps and run edstats
echo "[INFO] Starting post-processing (map generation and edstats)..."

# Find the postrefine directory
POSTREFINE_DIR="$OUTPUT_DIR/postrefine-$CompoundCode"

if [ -d "$POSTREFINE_DIR" ]; then
  MTZ_FILE="$POSTREFINE_DIR/refine.mtz"
  PDB_FILE="$POSTREFINE_DIR/refine.pdb"
  
  if [ -f "$MTZ_FILE" ] && [ -f "$PDB_FILE" ]; then
    # Generate 2fofc map
    echo "[INFO] Generating 2fofc map..."
    MAP_2FOFC="$POSTREFINE_DIR/refine_2fofc.map"
    gemmi sf2map --sample 5 "$MTZ_FILE" "$MAP_2FOFC" 2>&1
    
    if [ -f "$MAP_2FOFC" ]; then
      echo "[INFO] Successfully generated 2fofc map: $MAP_2FOFC"
    else
      echo "[WARNING] Failed to generate 2fofc map"
    fi
    
    # Generate fofc map
    echo "[INFO] Generating fofc (difference) map..."
    MAP_FOFC="$POSTREFINE_DIR/refine_fofc.map"
    gemmi sf2map --sample 5 -d "$MTZ_FILE" "$MAP_FOFC" 2>&1
    
    if [ -f "$MAP_FOFC" ]; then
      echo "[INFO] Successfully generated fofc map: $MAP_FOFC"
    else
      echo "[WARNING] Failed to generate fofc map"
    fi
    
    # Run edstats if both maps exist
    if [ -f "$MAP_2FOFC" ] && [ -f "$MAP_FOFC" ]; then
      echo "[INFO] Running edstats..."
      EDSTATS_OUT="$POSTREFINE_DIR/edstats.out"
      
      # Check if CCP4 is loaded (edstats is from CCP4)
      if command -v edstats &> /dev/null; then
        # Extract resolution values from pipedream_summary.json
        SUMMARY_JSON="$OUTPUT_DIR/pipedream_summary.json"
        
        if [ -f "$SUMMARY_JSON" ]; then
          # Use Python to extract resolution values from JSON
          RESOLUTION_VALUES=$(python3 -c "
import json
import sys
try:
    with open('$SUMMARY_JSON', 'r') as f:
        data = json.load(f)
    reslo = data.get('dataprocessing', {}).get('inputdata', {}).get('reslo', None)
    reshi = data.get('dataprocessing', {}).get('inputdata', {}).get('reshigh', None)
    if reslo is not None and reshi is not None:
        print(f'{reslo} {reshi}')
    else:
        sys.exit(1)
except Exception as e:
    sys.stderr.write(f'Error extracting resolution: {e}\n')
    sys.exit(1)
" 2>&1)
          
          if [ $? -eq 0 ] && [ -n "$RESOLUTION_VALUES" ]; then
            read RESLO RESHI <<< "$RESOLUTION_VALUES"
            echo "[INFO] Using resolution range: RESLO=$RESLO RESHI=$RESHI"
            
            # Run edstats with proper CCP4 syntax
            edstats XYZIN "$PDB_FILE" MAPIN1 "$MAP_2FOFC" MAPIN2 "$MAP_FOFC" OUT "$EDSTATS_OUT" << EOF
RESLO=$RESLO
RESHI=$RESHI
END
EOF
            EDSTATS_EXIT=$?
            
            if [ $EDSTATS_EXIT -eq 0 ] && [ -f "$EDSTATS_OUT" ]; then
              echo "[INFO] Successfully ran edstats: $EDSTATS_OUT"
            else
              echo "[WARNING] edstats failed with exit code: $EDSTATS_EXIT"
            fi
          else
            echo "[WARNING] Could not extract resolution values from $SUMMARY_JSON"
            echo "[INFO] Skipping edstats - resolution values required"
          fi
        else
          echo "[WARNING] Summary JSON not found: $SUMMARY_JSON"
          echo "[INFO] Skipping edstats - cannot determine resolution range"
        fi
      else
        echo "[WARNING] edstats command not found - skipping RSR calculation"
        echo "[INFO] Load CCP4 module if edstats is needed: module load ccp4"
      fi
    else
      echo "[WARNING] Cannot run edstats - maps not generated"
    fi
  else
    echo "[WARNING] MTZ or PDB file not found in postrefine directory"
  fi
else
  echo "[WARNING] Postrefine directory not found: $POSTREFINE_DIR"
fi

echo "[INFO] Post-processing complete"

# Cleanup setvar files
if compgen -G "__*.setvar.lis" > /dev/null; then
  mv __*.setvar.lis array_logs/ 2>/dev/null
fi
"""


def write_array_job_script(
    datasets_csv: str, 
    output_dir: str, 
    num_datasets: int, 
    script_path: str, 
    processing_dir: str, 
    refinement_params: Optional[Dict[str, Any]] = None,
    array_offset: int = 0,
    cluster_partition: str = DEFAULT_CLUSTER_PARTITION,
    job_priority: str = DEFAULT_JOB_PRIORITY
) -> None:
    """
    Write SLURM array job script for batch processing.
    
    Args:
        datasets_csv: Path to CSV file with dataset metadata
        output_dir: Base output directory
        num_datasets: Number of datasets in this chunk
        script_path: Where to write the script
        processing_dir: Base processing directory
        refinement_params: Refinement parameters dict (currently unused)
        array_offset: Offset for chunked job arrays (default 0)
        cluster_partition: Cluster partition to use (default: cs05r)
        job_priority: Job priority level (default: normal)
        
    Raises:
        ValueError: If num_datasets exceeds MAX_ARRAY_SIZE or is < 1
        ValueError: If datasets_csv doesn't exist
    """
    # Validation
    if num_datasets < 1:
        raise ValueError("No datasets to process: cannot write SLURM array job script.")
    
    if num_datasets > MAX_ARRAY_SIZE:
        raise ValueError(
            f"num_datasets ({num_datasets}) exceeds MAX_ARRAY_SIZE ({MAX_ARRAY_SIZE}). "
            "Use chunking in main()."
        )
    
    if not datasets_csv or not os.path.exists(datasets_csv):
        raise ValueError(f"Datasets CSV file does not exist: {datasets_csv}")
    
    # Ensure script directory exists
    script_dir = os.path.dirname(script_path)
    if not os.path.exists(script_dir):
        logging.info(f"Creating script directory: {script_dir}")
        os.makedirs(script_dir, exist_ok=True)
    
    # Log configuration
    logging.info(f"Writing SLURM script: {script_path}")
    logging.info(f"  Cluster partition: {cluster_partition}")
    logging.info(f"  Job priority: {job_priority}")
    logging.info(f"  Array range: 0-{num_datasets-1} (offset: {array_offset})")
    logging.info(f"  CSV line offset: {array_offset + 2}")
    logging.info(f"  Max concurrent: {min(MAX_CONCURRENT_JOBS, num_datasets)}")
    
    # Build script by concatenating sections
    script_sections = [
        generate_slurm_header(num_datasets, array_offset, cluster_partition, job_priority),
        generate_environment_setup(),
        generate_csv_parsing_section(datasets_csv, processing_dir, array_offset),
        generate_directory_setup_section(output_dir),
        generate_file_copy_section(processing_dir),
        generate_restraint_generation_section(),
        generate_pipedream_execution_section()
    ]
    script_content = '\n'.join(script_sections)
    
    # Write script with atomic operation
    try:
        temp_script_path = script_path + '.tmp'
        with open(temp_script_path, "w", encoding='utf-8') as f:
            f.write(script_content)
        
        # Validate temporary file
        if not os.path.exists(temp_script_path):
            raise IOError(f"Failed to create temporary script file: {temp_script_path}")
        
        temp_size = os.path.getsize(temp_script_path)
        if temp_size == 0:
            raise ValueError(f"Temporary script file is empty: {temp_script_path}")
        
        # Atomic move
        shutil.move(temp_script_path, script_path)
        
        # Final validation
        final_size = os.path.getsize(script_path)
        logging.info(f"SLURM script written successfully: {script_path} ({final_size} bytes)")
        
    except Exception as e:
        logging.error(f"Error writing SLURM script to {script_path}: {e}")
        if os.path.exists(temp_script_path):
            try:
                os.remove(temp_script_path)
            except:
                pass
        raise

def ensure_remote_directory(sftp: paramiko.SFTPClient, remote_dir: str) -> None:
    """
    Ensure remote directory exists, creating it if necessary.
    
    Args:
        sftp: Active SFTP client
        remote_dir: Remote directory path to create/verify
        
    Raises:
        IOError: If directory cannot be created
    """
    try:
        sftp.chdir(remote_dir)
        logging.debug(f"Remote directory exists: {remote_dir}")
    except IOError:
        # Directory does not exist, create it recursively
        logging.info(f"Creating remote directory: {remote_dir}")
        parts = remote_dir.strip('/').split('/')
        path = ''
        for part in parts:
            path += '/' + part
            try:
                sftp.chdir(path)
            except IOError:
                sftp.mkdir(path)
                sftp.chdir(path)
        logging.info(f"Successfully created remote directory: {remote_dir}")


def copy_file_to_remote(
    sftp: paramiko.SFTPClient, 
    local_path: str, 
    remote_path: str
) -> None:
    """
    Copy a local file to remote location via SFTP.
    
    Args:
        sftp: Active SFTP client
        local_path: Local file path to copy
        remote_path: Remote destination path
        
    Raises:
        IOError: If file copy fails
    """
    remote_name = os.path.basename(remote_path)
    print(f"\nCopying {remote_name} to wilson...")
    
    try:
        sftp.put(local_path, remote_path)
        logging.info(f"File copied to wilson:{remote_path}")
    except Exception as e:
        logging.error(f"Failed to copy file to wilson: {e}")
        raise


def execute_sbatch_command(
    ssh_client: paramiko.SSHClient, 
    script_name: str, 
    working_dir: str
) -> Tuple[int, str, str]:
    """
    Execute sbatch command on remote host and return results.
    
    Args:
        ssh_client: Active SSH client
        script_name: Name of the script to submit
        working_dir: Working directory for sbatch execution
        
    Returns:
        Tuple of (exit_status, stdout, stderr)
    """
    command = f"cd {working_dir} && sbatch {script_name}"
    logging.debug(f"Executing remote command: {command}")
    
    print(f"Submitting job...")
    stdin, stdout, stderr = ssh_client.exec_command(command)
    
    # Wait for command to complete and get results
    exit_status = stdout.channel.recv_exit_status()
    out = stdout.read().decode().strip()
    err = stderr.read().decode().strip()
    
    logging.info(f"sbatch exit status: {exit_status}")
    if out:
        logging.info(f"sbatch stdout: {out}")
    if err:
        # Filter out harmless module errors
        harmless_error = "Module ERROR: no such variable"
        if harmless_error in err:
            logging.debug(f"sbatch harmless stderr (suppressed): {err}")
        else:
            logging.warning(f"sbatch stderr: {err}")
    
    return exit_status, out, err


def validate_sbatch_result(exit_status: int, stdout: str, stderr: str) -> None:
    """
    Validate sbatch execution results and raise errors if needed.
    
    Args:
        exit_status: Exit code from sbatch
        stdout: Standard output from sbatch
        stderr: Standard error from sbatch
        
    Raises:
        RuntimeError: If sbatch command failed
    """
    if exit_status != 0:
        logging.error(f"sbatch failed with exit status {exit_status}")
        if stderr:
            logging.error(f"sbatch error output: {stderr}")
        print(f"ERROR: sbatch failed with exit status {exit_status}")
        if stderr:
            print(f"Error output: {stderr}")
        raise RuntimeError(f"sbatch command failed with exit status {exit_status}")
    
    if not stdout:
        logging.error("sbatch produced no output")
        print("ERROR: sbatch produced no output")
        raise RuntimeError("sbatch command produced no output")
    
    # Success!
    print(f"✓ Job submitted successfully!")
    print(f"{stdout}")


def submit_sbatch_on_wilson(
    script_path: str, 
    processing_dir: str, 
    ssh_client: paramiko.SSHClient
) -> None:
    """
    Copy script to wilson and submit via sbatch.
    
    This orchestrates the complete submission workflow:
    1. Open SFTP connection
    2. Ensure remote directory exists
    3. Copy script file to remote
    4. Execute sbatch command
    5. Validate and report results
    
    Args:
        script_path: Local path to the job script
        processing_dir: Remote directory to copy script to
        ssh_client: Existing paramiko SSH client connection
        
    Raises:
        RuntimeError: If submission fails
        IOError: If file operations fail
    """
    remote_script_name = os.path.basename(script_path)
    remote_script_path = os.path.join(processing_dir, remote_script_name)
    
    try:
        # Step 1: Open SFTP connection
        sftp = ssh_client.open_sftp()
        
        try:
            # Step 2: Ensure remote directory exists
            ensure_remote_directory(sftp, processing_dir)
            
            # Step 3: Copy script to remote
            copy_file_to_remote(sftp, script_path, remote_script_path)
            
        finally:
            # Always close SFTP connection
            sftp.close()
        
        # Step 4: Execute sbatch command
        exit_status, stdout, stderr = execute_sbatch_command(
            ssh_client, 
            remote_script_name, 
            processing_dir
        )
        
        # Step 5: Validate and report results
        validate_sbatch_result(exit_status, stdout, stderr)
        
    except Exception as e:
        logging.error(f"Failed to submit job on wilson: {e}")
        raise


def establish_wilson_connection() -> paramiko.SSHClient:
    """Establish SSH connection with automatic key discovery or password fallback."""
    user = os.environ.get("CLUSTER_USER", os.getlogin())
    wilson_host = "wilson.diamond.ac.uk"
    
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    
    # Try key-based authentication first (let paramiko auto-discover keys)
    try:
        logging.info(f"Attempting SSH connection to {wilson_host} as {user}...")
        ssh.connect(
            wilson_host,
            username=user,
            timeout=30,
            banner_timeout=30,
            look_for_keys=True,  # Auto-discover SSH keys
            allow_agent=True     # Use SSH agent if available
        )
        logging.info("SSH connection established using key-based authentication")
        return ssh
        
    except paramiko.AuthenticationException:
        logging.info("Key-based authentication failed, trying password authentication...")
        
        # Fallback to password authentication
        try:
            password = getpass.getpass(f"Password for {user}@{wilson_host}: ")
            ssh.connect(
                wilson_host,
                username=user,
                password=password,
                timeout=30,
                banner_timeout=30
            )
            logging.info("SSH connection established using password authentication")
            return ssh
            
        except paramiko.AuthenticationException:
            raise RuntimeError(
                f"SSH authentication failed for {user}@{wilson_host}\n"
                "Please check your credentials or SSH key configuration."
            )
    
    except paramiko.SSHException as e:
        raise RuntimeError(f"SSH connection failed: {e}")
    except Exception as e:
        raise RuntimeError(f"Unexpected error connecting to wilson: {e}")

def resolve_dataset_paths(
    crystal_name: str,
    refinement_mtzfree: str,
    dimple_path_pdb: str,
    dimple_path_mtz: str,
    base_dir: str
) -> Optional[Dict[str, str]]:
    """
    Resolve and validate all file paths for a dataset.
    
    Handles inconsistent database paths - works with:
    - Full absolute paths
    - Relative paths (joins with base_dir)
    - Just filenames (assumes in model_building subdir)
    
    Args:
        crystal_name: Crystal identifier
        refinement_mtzfree: RefinementMTZfree path from database
        dimple_path_pdb: DimplePathToPDB from database
        dimple_path_mtz: DimplePathToMTZ from database
        base_dir: Base processing directory
        
    Returns:
        Dictionary with resolved paths or None if validation fails
        Keys: 'refinement_mtzfree', 'dimple_pdb', 'dimple_mtz'
    """
    # Log raw values from database
    logging.debug(f"Raw database values for {crystal_name}:")
    logging.debug(f"  RefinementMTZfree: {repr(refinement_mtzfree)}")
    logging.debug(f"  DimplePathToPDB: {repr(dimple_path_pdb)}")
    logging.debug(f"  DimplePathToMTZ: {repr(dimple_path_mtz)}")
    
    # Validate we have actual values
    if not refinement_mtzfree or pd.isna(refinement_mtzfree):
        logging.error(f"Missing RefinementMTZfree for {crystal_name}")
        return None
    
    if not dimple_path_pdb or pd.isna(dimple_path_pdb):
        logging.error(f"Missing DimplePathToPDB for {crystal_name}")
        return None
        
    if not dimple_path_mtz or pd.isna(dimple_path_mtz):
        logging.error(f"Missing DimplePathToMTZ for {crystal_name}")
        return None
    
    # Resolve RefinementMTZfree: handle both filename and full path
    refinement_mtzfree_clean = str(refinement_mtzfree).strip()
    if os.path.isabs(refinement_mtzfree_clean):
        # Database has full path - extract just the filename
        logging.debug(f"RefinementMTZfree is absolute path, extracting basename")
        refinement_mtzfree_filename = os.path.basename(refinement_mtzfree_clean)
    else:
        # Database has relative path or filename - use as-is
        refinement_mtzfree_filename = refinement_mtzfree_clean
    
    refinement_mtzfree_path = resolve_path(
        base_dir, 
        "analysis", "model_building", crystal_name, 
        refinement_mtzfree_filename,
        must_exist=True
    )
    
    # Resolve Dimple paths: handle both absolute and relative paths
    dimple_path_pdb_clean = str(dimple_path_pdb).strip()
    if os.path.isabs(dimple_path_pdb_clean):
        logging.debug(f"DimplePathToPDB is absolute path")
        dimple_path_pdb_resolved = dimple_path_pdb_clean if os.path.exists(dimple_path_pdb_clean) else None
    else:
        dimple_path_pdb_resolved = resolve_path(base_dir, dimple_path_pdb_clean, must_exist=True)
    
    dimple_path_mtz_clean = str(dimple_path_mtz).strip()
    if os.path.isabs(dimple_path_mtz_clean):
        logging.debug(f"DimplePathToMTZ is absolute path")
        dimple_path_mtz_resolved = dimple_path_mtz_clean if os.path.exists(dimple_path_mtz_clean) else None
    else:
        dimple_path_mtz_resolved = resolve_path(base_dir, dimple_path_mtz_clean, must_exist=True)
    
    # Validate all paths exist
    missing_paths = []
    if not refinement_mtzfree_path:
        expected = f"{base_dir}/analysis/model_building/{crystal_name}/{refinement_mtzfree_filename}"
        missing_paths.append(f"RefinementMTZfree (tried: {expected})")
    if not dimple_path_pdb_resolved:
        missing_paths.append(f"DimplePathToPDB (tried: {dimple_path_pdb_clean})")
    if not dimple_path_mtz_resolved:
        missing_paths.append(f"DimplePathToMTZ (tried: {dimple_path_mtz_clean})")
    
    if missing_paths:
        logging.error(f"Missing required file paths for {crystal_name}:")
        for path in missing_paths:
            logging.error(f"{path}")
        return None
    
    # Log resolved paths
    logging.debug(f"Resolved paths for {crystal_name}:")
    logging.debug(f"RefinementMTZfree: {refinement_mtzfree_path}")
    logging.debug(f"DimplePathToPDB: {dimple_path_pdb_resolved}")
    logging.debug(f"DimplePathToMTZ: {dimple_path_mtz_resolved}")
    
    return {
        'refinement_mtzfree': refinement_mtzfree_path,
        'dimple_pdb': dimple_path_pdb_resolved,
        'dimple_mtz': dimple_path_mtz_resolved
    }


def build_pipedream_command(
    paths: Dict[str, str],
    input_dir: str,
    output_dir: str,
    compound_code: str,
    refinement_params: Dict[str, Any]
) -> str:
    """
    Build Pipedream command from resolved paths and parameters.
    
    Args:
        paths: Dictionary with resolved file paths (from resolve_dataset_paths)
        input_dir: Input files directory
        output_dir: Pipedream output directory
        compound_code: Compound code
        refinement_params: Refinement parameters dictionary
        
    Returns:
        Complete Pipedream command string
    """
    pipedream_path, _ = get_buster_path('pipedream')
    
    # Build file paths relative to input_dir
    hklin_file = os.path.join(input_dir, os.path.basename(paths['refinement_mtzfree']))
    xyzin_file = os.path.join(input_dir, os.path.basename(paths['dimple_pdb']))
    hklref_file = os.path.join(input_dir, os.path.basename(paths['dimple_mtz']))
    expected_cif = os.path.join(input_dir, f"{compound_code}.cif")
    
    # Build base command
    pipedream_cmd = (
        f"{pipedream_path} -nolmr "
        f"-hklin {hklin_file} "
        f"-xyzin {xyzin_file} "
        f"-hklref {hklref_file} "
        f"-d {output_dir} "
    )
    
    # Add refinement arguments
    refinement_args = refinement_params_to_args(
        refinement_params,
        rhofit_file=expected_cif,
        rhocommands=refinement_params.get('rhocommands')
    )
    
    if refinement_args:
        pipedream_cmd += refinement_args
    
    return pipedream_cmd.strip()


def create_dataset_metadata(
    input_dir: str,
    output_dir: str,
    compound_code: str,
    smiles_string: str,
    pipedream_cmd: str
) -> Dict[str, Any]:
    """
    Create metadata dictionary for a processed dataset.
    
    Args:
        input_dir: Input files directory
        output_dir: Pipedream output directory
        compound_code: Compound code
        smiles_string: SMILES string
        pipedream_cmd: Full Pipedream command
        
    Returns:
        Dictionary with dataset metadata
    """
    return {
        'Input_dir': input_dir,
        'CompoundCode': compound_code,
        'PipedreamDirectory': output_dir,
        'ReportHTML': f"{output_dir}/report-{compound_code}/index.html",
        'LigandReportHTML': f"{output_dir}/report-{compound_code}/ligand/index.html",
        'ExpectedSummary': f"{output_dir}/pipedream_summary.json",
        'PipedreamCommand': pipedream_cmd,
        'ExpectedCIF': os.path.join(input_dir, f"{compound_code}.cif"),
        'ExpectedPDB': os.path.join(input_dir, f"{compound_code}.pdb"),
        'InputSMILES': smiles_string
    }


# Separate logic from I/O
def validate_dataset_fields(row: pd.Series) -> Tuple[bool, List[str]]:
    """Pure validation logic - easy to test."""
    missing = []
    if pd.isna(row.get('CrystalName')):
        missing.append('CrystalName')
    # ... more checks
    return len(missing) == 0, missing

def process_single_dataset(row, params, parent_dir, timestamp):
    """Coordinates I/O - harder to test but smaller."""
    is_valid, errors = validate_dataset_fields(row)
    if not is_valid:
        logging.error(f"Invalid dataset: {errors}")
        return None
    # Extract dataset fields
    crystal_name, compound_code, refinement_mtzfree, dimple_path_pdb, dimple_path_mtz, compound_smiles = extract_dataset_fields(row)
    
    logging.debug(f"Processing dataset {crystal_name} (compound: {compound_code})")
    
    # Validate crystal name
    if not crystal_name or pd.isna(crystal_name):
        logging.error("Skipping dataset with missing CrystalName")
        return None
    
    # Resolve and validate all file paths
    paths = resolve_dataset_paths(
        crystal_name,
        refinement_mtzfree,
        dimple_path_pdb,
        dimple_path_mtz,
        params['Processing_directory']
    )
    
    if not paths:
        logging.error(f"Failed to resolve paths for {crystal_name}")
        return None
    
    # Setup directories
    dataset_dir, input_dir = setup_dataset_directory(crystal_name, parent_dir)
    
    # Get and save SMILES
    smiles_string = get_smiles_for_compound(
        compound_code,
        compound_smiles,
        input_dir
    )
    
    if not smiles_string:
        logging.error(f"Skipping {crystal_name} - no valid SMILES available")
        return None
    
    # Generate output directory
    output_dir = f"{dataset_dir}/Pipedream_{timestamp}"
    
    # Copy files and prepare for restraint generation
    files_ready = copy_input_files_and_prepare_for_restraints(
        input_dir,
        paths['refinement_mtzfree'],
        paths['dimple_pdb'],
        paths['dimple_mtz'],
        compound_code,
        smiles_string
    )
    
    if not files_ready:
        logging.error(f"Failed to prepare input files for {crystal_name}")
        return None
    
    # Process PDB file if needed
    process_pdb_file(input_dir, paths['dimple_pdb'], crystal_name, params)
    
    # Build Pipedream command
    pipedream_cmd = build_pipedream_command(
        paths,
        input_dir,
        output_dir,
        compound_code,
        params.get('Refinement_parameters', {})
    )
    
    # Create and return metadata
    return create_dataset_metadata(
        input_dir,
        output_dir,
        compound_code,
        smiles_string,
        pipedream_cmd
    )

def get_buster_path(command_name: str) -> Tuple[str, str]:
    """
    Get the path to a BUSTER command, preferring module-loaded version.
    
    Args:
        command_name: Name of the BUSTER command (e.g., 'buster-report', 'pipedream')
        
    Returns:
        Tuple of (full path to command, source description: 'PATH', 'fallback', or 'not_found')
    """
    import subprocess
    
    try:
        # Try to find the command in PATH (assumes module load buster has been run)
        result = subprocess.run(
            ['which', command_name],
            capture_output=True,
            text=True,
            timeout=5
        )
        
        if result.returncode == 0 and result.stdout.strip():
            path = result.stdout.strip()
            logging.debug(f"Found {command_name} in PATH (from module): {path}")
            return path, 'PATH'
    except Exception as e:
        logging.debug(f"Could not search PATH for {command_name}: {e}")
    
    # Fallback to hardcoded paths
    fallback_paths = {
        'buster-report': BUSTER_REPORT_PATH_FALLBACK,
        'pipedream': BUSTER_PIPEDREAM_PATH_FALLBACK
    }
    
    fallback_path = fallback_paths.get(command_name)
    if fallback_path and os.path.exists(fallback_path):
        logging.debug(f"Using fallback path for {command_name}: {fallback_path}")
        return fallback_path, 'fallback'
    elif fallback_path:
        logging.warning(f"Fallback path for {command_name} does not exist: {fallback_path}")
    
    # Last resort: return command name and hope it's in PATH
    logging.warning(f"Could not find {command_name}, will try to use from PATH")
    return command_name, 'not_found'


def check_buster_dependencies(verbose: bool = False) -> bool:
    """Check if all BUSTER dependencies are available after loading module."""
    try:
        import subprocess
        print("\nChecking BUSTER dependencies...")
        
        # Step 1: Try to load buster module
        print("Loading BUSTER module...")
        logging.info("Attempting to load buster module")
        
        try:
            # Source /etc/profile and load buster module
            load_module_cmd = "source /etc/profile 2>/dev/null && module load buster 2>&1"
            module_result = subprocess.run(
                load_module_cmd,
                shell=True,
                executable='/bin/bash',
                capture_output=True,
                text=True,
                timeout=10
            )
            
            if module_result.returncode == 0:
                print("✓ BUSTER module loaded successfully")
                logging.info("BUSTER module loaded successfully")
            else:
                print(f"⚠️  Warning: module load buster returned non-zero exit code")
                logging.warning(f"module load buster stderr: {module_result.stderr}")
                if verbose:
                    print(f"   Output: {module_result.stderr}")
                
        except subprocess.TimeoutExpired:
            print("⚠️  Warning: module load buster timed out")
            logging.warning("module load buster timed out")
        except Exception as e:
            print(f"⚠️  Warning: Could not load buster module: {e}")
            logging.warning(f"Could not load buster module: {e}")
        
        # Step 2: Check dependencies using buster-report
        print("Running dependency check...")
        logging.info("Running buster-report -checkdep")
        
        # Run buster-report -checkdep in a shell with module loaded
        checkdep_cmd = "source /etc/profile 2>/dev/null && module load buster 2>/dev/null && buster-report -checkdep"
        result = subprocess.run(
            checkdep_cmd,
            shell=True,
            executable='/bin/bash',
            capture_output=True,
            text=True,
            timeout=30
        )
        
        output = result.stderr.strip() if result.stderr else result.stdout.strip()
        
        if verbose:
            print("\n--- BUSTER Dependency Check Output ---")
            print(output)
            print("--- End of Output ---")
        
        if not verbose:
            logging.debug(f"buster-report -checkdep output:\n{output}")
        
        if result.returncode == 0:
            print("✓ All BUSTER dependencies are satisfied\n")
            logging.info("All BUSTER dependencies satisfied")
            return True
        else:
            if not verbose:
                print("\n--- BUSTER Dependency Check Output (FAILED) ---")
                print(output)
                print("--- End of Output ---\n")
            
            print("!"*70)
            print("! WARNING: BUSTER DEPENDENCY CHECK FAILED")
            print("!"*70)
            
            # Try to identify specific missing dependencies
            missing_deps = []
            for line in output.split('\n'):
                if 'not found' in line.lower() or 'missing' in line.lower() or 'error' in line.lower() or 'undefined' in line.lower():
                    missing_deps.append(line.strip())
            
            if missing_deps:
                print("\nMissing or problematic dependencies:")
                for dep in missing_deps:
                    print(f"  ! {dep}")
                logging.warning(f"Missing dependencies: {missing_deps}")
            
            print("\n" + "!"*70)
            print("! To resolve these issues:")
            print("!   1. Ensure you're on a system with module environment configured")
            print("!   2. Check that BUSTER module is available: 'module avail buster'")
            print("!   3. The cluster compute nodes will load BUSTER automatically")
            print("!"*70 + "\n")
            
            return False
            
    except FileNotFoundError:
        print("\n" + "!"*70)
        print("! ERROR: Could not run buster-report")
        print("! ")
        print("! This could mean:")
        print("!   1. BUSTER module failed to load")
        print("!   2. buster-report not in PATH after module load")
        print("!   3. Module system not available on this host")
        print("! ")
        print("! Note: Cluster compute nodes will load BUSTER automatically")
        print("!"*70 + "\n")
        logging.error("buster-report not found after attempting module load")
        return False
    except subprocess.TimeoutExpired:
        print("\n" + "!"*70)
        print("! ERROR: buster-report -checkdep timed out after 30 seconds")
        print("!"*70 + "\n")
        logging.error("buster-report -checkdep timed out")
        return False
    except Exception as e:
        print(f"\n" + "!"*70)
        print(f"! ERROR checking BUSTER dependencies: {e}")
        print("!"*70 + "\n")
        logging.error(f"Error checking BUSTER dependencies: {e}")
        return False
    

def submit_single_job(
    parent_dir: str,
    num_datasets: int,
    datasets_csv: str,
    processing_dir: str,
    params: Dict[str, Any],
    ssh_client: paramiko.SSHClient
) -> None:
    """
    Submit a single SLURM array job (for datasets <= MAX_ARRAY_SIZE).
    
    Args:
        parent_dir: Output directory
        num_datasets: Number of datasets to process
        datasets_csv: Path to CSV with dataset metadata
        processing_dir: Base processing directory
        params: Full parameters dict (includes refinement params and cluster config)
        ssh_client: Active SSH connection to wilson
    """
    script_path = os.path.join(parent_dir, "pipedream_array.sh")
    
    # Extract cluster configuration with defaults
    cluster_partition = params.get('Cluster_partition', DEFAULT_CLUSTER_PARTITION)
    job_priority = params.get('Job_priority', DEFAULT_JOB_PRIORITY)
    
    logging.info(f"Writing single job script for {num_datasets} datasets")
    write_array_job_script(
        datasets_csv=datasets_csv,
        output_dir=parent_dir,
        num_datasets=num_datasets,
        script_path=script_path,
        processing_dir=processing_dir,
        refinement_params=params.get('Refinement_parameters'),
        cluster_partition=cluster_partition,
        job_priority=job_priority
    )
    
    submit_sbatch_on_wilson(
        script_path=script_path,
        processing_dir=processing_dir,
        ssh_client=ssh_client
    )
    
    print(f"\n{format_separator()}")
    print(f"Summary:")
    print(f"Job submitted successfully")
    print(f"Datasets: {num_datasets}")
    print(f"Output directory: {parent_dir}")
    print(f"Script: {script_path}")
    print(f"{format_separator()}\n")
    
    logging.info(f"Single job submitted: {num_datasets} datasets")


def submit_chunked_jobs(
    parent_dir: str,
    num_datasets: int,
    datasets_csv: str,
    processing_dir: str,
    params: Dict[str, Any],
    ssh_client: paramiko.SSHClient
) -> None:
    """
    Submit multiple SLURM array jobs (for datasets > MAX_ARRAY_SIZE).
    
    Each chunk processes up to MAX_ARRAY_SIZE datasets. This is necessary
    because SLURM has limits on array job sizes.
    
    Args:
        parent_dir: Output directory
        num_datasets: Total number of datasets to process
        datasets_csv: Path to CSV with dataset metadata
        processing_dir: Base processing directory
        params: Full parameters dict (includes refinement params and cluster config)
        ssh_client: Active SSH connection to wilson
    """
    num_chunks = (num_datasets + MAX_ARRAY_SIZE - 1) // MAX_ARRAY_SIZE
    
    # Extract cluster configuration with defaults
    cluster_partition = params.get('Cluster_partition', DEFAULT_CLUSTER_PARTITION)
    job_priority = params.get('Job_priority', DEFAULT_JOB_PRIORITY)
    
    print(f"\n⚠️  Large dataset: splitting {num_datasets} datasets into {num_chunks} job chunks")
    print(f"   (Each chunk will process up to {MAX_ARRAY_SIZE} datasets)\n")
    
    logging.info(f"Chunking {num_datasets} datasets into {num_chunks} jobs")
    
    for chunk_idx in range(num_chunks):
        start_idx = chunk_idx * MAX_ARRAY_SIZE
        end_idx = min(start_idx + MAX_ARRAY_SIZE - 1, num_datasets - 1)
        chunk_size = end_idx - start_idx + 1
        
        print(f"Preparing chunk {chunk_idx + 1}/{num_chunks}:")
        print(f"  Datasets: {start_idx} to {end_idx} (total: {chunk_size})")
        
        script_path = os.path.join(parent_dir, f"pipedream_array_chunk_{chunk_idx}.sh")
        
        # Write script with array_offset for this chunk
        logging.info(f"Writing chunk {chunk_idx + 1} script: array_offset={start_idx}, size={chunk_size}")
        write_array_job_script(
            datasets_csv=datasets_csv,
            output_dir=parent_dir,
            num_datasets=chunk_size,
            script_path=script_path,
            processing_dir=processing_dir,
            refinement_params=params.get('Refinement_parameters'),
            array_offset=start_idx,
            cluster_partition=cluster_partition,
            job_priority=job_priority
        )
        
        print(f"Submitting chunk {chunk_idx + 1}...")
        submit_sbatch_on_wilson(
            script_path=script_path,
            processing_dir=processing_dir,
            ssh_client=ssh_client
        )
        print(f"✓ Chunk {chunk_idx + 1} submitted\n")
        
        logging.info(f"Chunk {chunk_idx + 1}/{num_chunks} submitted successfully")

    print(f"\n{format_separator()}")
    print(f"Successfully submitted {num_chunks} job chunks")
    print(f"Total datasets: {num_datasets}")
    print(f"Output directory: {parent_dir}")
    print(f"Chunk size: {MAX_ARRAY_SIZE} datasets/chunk")
    print(f"{format_separator()}\n")
    
    logging.info(f"All {num_chunks} chunks submitted successfully")


def print_dry_run_summary(num_datasets: int, json_output_path: str) -> None:
    """
    Print summary of dry-run results.
    
    Args:
        num_datasets: Number of datasets that would be processed
        json_output_path: Path to generated JSON metadata file
    """
    print(f"\nDry run complete.")
    
    if num_datasets > 0:
        print(f"✓ {num_datasets} datasets would be submitted for processing")
        if num_datasets > MAX_ARRAY_SIZE:
            num_chunks = (num_datasets + MAX_ARRAY_SIZE - 1) // MAX_ARRAY_SIZE
            print(f"  Would be split into {num_chunks} job chunks")
            print(f"  Each chunk processes up to {MAX_ARRAY_SIZE} datasets")
    else:
        print("✗ No datasets would be submitted (all filtered out)")

    logging.debug(f"Dry run complete: {num_datasets} datasets prepared")
    print(f"{format_separator()}\n")


def log_and_print(message: str, level: str = "INFO") -> None:
    """Log message and print to console."""
    getattr(logging, level.lower())(message)
    print(message)

def format_separator(width: int = 76, char: str = '=') -> str:
    return char * width

@contextmanager
def wilson_ssh_connection():
    """Context manager for SSH connection lifecycle."""
    ssh = None
    try:
        ssh = establish_wilson_connection()
        yield ssh
    finally:
        if ssh:
            try:
                ssh.close()
                logging.info("SSH connection closed")
            except Exception as e:
                logging.warning(f"Error closing SSH connection: {e}")

if __name__ == "__main__":
    try:
        args = parse_args()
        
        # Read parameters first to get Output_directory
        params = read_yaml(args.parameters)
        validate_params(params)
        
        # Determine log directory from Output_directory or generate timestamp-based path
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        default_output_dir = f"{params['Processing_directory']}/analysis/Pipedream/Pipedream_{timestamp}"
        log_dir = params.get('Output_directory') or default_output_dir
        
        # Ensure log directory exists
        os.makedirs(log_dir, exist_ok=True)
        
        # Setup logging with the determined directory
        setup_logging(log_dir=log_dir, log_level=args.log_level, verbose=args.verbose)
        
        logging.info(f"Pipedream XChem Integration Script v{VERSION}")
        logging.info(f"Log directory: {log_dir}")
        logging.info(f"Parameters file: {args.parameters}")
        
        # Check BUSTER dependencies
        deps_ok = check_buster_dependencies(verbose=args.verbose)

        if not deps_ok:
            if args.dry_run:
                print("⚠️  BUSTER dependency check failed, but continuing with dry-run.")
                print("Note: Cluster compute nodes will load BUSTER modules automatically.\n")
            else:
                print("\n⚠️  BUSTER dependency check failed!")
                print("\nNote: Cluster compute nodes will load BUSTER modules automatically,")
                print("but this may indicate configuration issues that could affect job execution.")
                print("\nDo you want to continue anyway? (y/N): ", end='')
                response = input().strip().lower()
                if response not in ['y', 'yes']:
                    print("Exiting due to dependency issues.")
                    sys.exit(1)
                print("Continuing despite dependency issues...\n")
        
        print("Getting datasets from database...")
        datasets = get_datasets(params)
        
        if datasets is None or datasets.empty:
            print("No datasets found or retrieved from database.")
            logging.warning("No datasets available for processing")
            sys.exit(0)
        
        try:
            # Submit jobs based on mode
            if args.dry_run:
                parent_dir, num_valid_datasets, json_output_path = submit_jobs(
                    datasets, 
                    params, 
                    output_json_path=args.output_json,
                    output_csv_path=args.output_csv,
                    dry_run=args.dry_run,
                    log_level=args.log_level,
                    verbose=args.verbose
                )
                print_dry_run_summary(num_valid_datasets, json_output_path)
            else:
                with wilson_ssh_connection() as ssh_connection:
                    parent_dir, num_valid_datasets, json_output_path = submit_jobs(
                        datasets, 
                        params, 
                        output_json_path=args.output_json,
                        output_csv_path=args.output_csv,
                        dry_run=args.dry_run,
                        log_level=args.log_level,
                        verbose=args.verbose
                    )
                    
                    if num_valid_datasets > MAX_ARRAY_SIZE:
                        submit_chunked_jobs(
                            parent_dir=parent_dir,
                            num_datasets=num_valid_datasets,
                            datasets_csv=params['Filtered_dataset_csv_path'],
                            processing_dir=params['Processing_directory'],
                            params=params,
                            ssh_client=ssh_connection
                        )
                    else:
                        submit_single_job(
                            parent_dir=parent_dir,
                            num_datasets=num_valid_datasets,
                            datasets_csv=params['Filtered_dataset_csv_path'],
                            processing_dir=params['Processing_directory'],
                            params=params,
                            ssh_client=ssh_connection
                        )
        
        except KeyboardInterrupt:
            print("\n\nOperation cancelled by user.")
            logging.info("Operation cancelled by user (KeyboardInterrupt)")
            sys.exit(1)
        except Exception as e:
            print(f"\n{'='*70}")
            print(f"ERROR: {e}")
            print(f"{'='*70}")
            import traceback
            traceback.print_exc()
            logging.error(f"Fatal error: {e}", exc_info=True)
            sys.exit(1)

    except KeyboardInterrupt:
        print("\n\nProcess interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"Fatal error: {e}")
        sys.exit(1)
