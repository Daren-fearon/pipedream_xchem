"""
XChem Pipedream Results Export Script

This script exports completed Pipedream refinement results to XChem Explorer (XCE) 
format and updates the associated SQLite database. It processes validated results 
from the collation script, integrates them into the XCE file structure, and handles 
special cases like chirality inversions by updating restraint files and molecular 
diagrams.

Features:
- Exports Pipedream results to XCE-compatible directory structure
- Updates SQLite database with refinement statistics and quality metrics
- Applies traffic-light scoring for validation metrics
- Handles chirality inversions by copying updated restraint files (CIF/PDB)
- Generates molecular structure diagrams (PNG) from SMILES using RDKit
- Creates proper symlinks for PDB, MTZ, and map files in XCE format
- Copies BUSTER output files for downstream analysis

Author: DFearon
Date: November 2025
"""

import os
import json
import shutil
import sqlite3
import yaml
import csv
import argparse
import getpass
import logging
from datetime import datetime

# Version information
VERSION = "1.0.2"

# Import RDKit for PNG generation
try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not available - PNG generation will be disabled")

# Configure logging to file and console
def setup_logging():
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"export_{timestamp}.log"
    
    # Create logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    
    # Create file handler for detailed logging
    file_handler = logging.FileHandler(log_filename)
    file_handler.setLevel(logging.INFO)
    file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(file_formatter)
    
    # Create console handler for minimal output
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console_handler.setFormatter(console_formatter)
    
    # Add handlers to logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    
    return logger, log_filename

# Get logger for module-level use
logger, log_filename = setup_logging()

# Function to log only to file (not console)
def log_to_file_only(message):
    # Create a temporary logger that only has the file handler
    file_logger = logging.getLogger('file_only')
    file_logger.setLevel(logging.INFO)
    # Clear any existing handlers
    file_logger.handlers.clear()
    # Add only the file handler
    file_handler = logging.getLogger().handlers[0]  # First handler is file handler
    file_logger.addHandler(file_handler)
    file_logger.propagate = False  # Don't propagate to parent logger
    file_logger.info(message)

# Helper function for crystal-specific logging
def log_crystal_progress(crystal, message):
    """Log a message with crystal context to file only."""
    log_to_file_only(f"Processing crystal {crystal}: {message}")

# Helper function to safely create symlinks
def safe_symlink(src, dst):
    try:
        if os.path.islink(dst) or os.path.exists(dst):
            os.remove(dst)
        os.symlink(src, dst)
    except Exception as e:
        logging.error(f"Error creating symlink from {src} to {dst}: {e}")

# Helper function to check if a column exists in a table
def column_exists(cursor, table, column):
    cursor.execute(f"PRAGMA table_info({table})")
    columns = [info[1] for info in cursor.fetchall()]
    return column in columns

# Traffic light function
def traffic_light(value, green, orange=None, reverse=False):
    try:
        if value is None or value == "" or value == "NA":
            return None
        val = float(value)
        if orange:
            if reverse:
                # For metrics where higher is better (e.g., Ramachandran Favored)
                if val > green:
                    return "green"
                elif val > orange:
                    return "orange"
                else:
                    return "red"
            else:
                # For metrics where lower is better (e.g., R-factor, resolution)
                if val < green:
                    return "green"
                elif val < orange:
                    return "orange"
                else:
                    return "red"
        else:
            if reverse:
                return "green" if val > green else "red"
            else:
                return "green" if val < green else "red"
    except (ValueError, TypeError) as e:
        logging.warning(f"Invalid value for traffic light calculation: {value} (error: {e})")
        return None

# Load dataset list from CSV
def load_dataset_list(csv_path):
    with open(csv_path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        return [row[0] for row in reader if row]

# Function to generate PNG from SMILES using RDKit
def generate_png_from_smiles(smiles, output_path):
    """Generate a 2D structure diagram PNG from SMILES string."""
    if not RDKIT_AVAILABLE:
        logging.warning("RDKit not available - cannot generate PNG")
        return False
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logging.error(f"Invalid SMILES string: {smiles}")
            return False
        
        # Generate 2D coordinates
        Chem.rdDepictor.Compute2DCoords(mol)
        
        # Create image (approximately 3 inch x 3 inch at 100 DPI)
        img = Draw.MolToImage(mol, size=(300, 300), kekulize=True)
        img.save(output_path)
        
        log_to_file_only(f"Generated PNG structure diagram: {output_path}")
        return True
        
    except Exception as e:
        logging.error(f"Error generating PNG from SMILES '{smiles}': {e}")
        return False

# Helper function to safely copy files
def safe_copy_file(src, dst, crystal_name):
    """Copy file to dst using a temp-file-and-rename strategy to avoid partial writes
    and to provide clearer logging/fallbacks when permission errors occur."""
    try:
        if not src or not os.path.exists(src):
            raise FileNotFoundError(src)

        dst_dir = os.path.dirname(dst)
        os.makedirs(dst_dir, exist_ok=True)

        # Copy to a temporary file in the destination directory then rename
        tmp_dst = dst + ".tmp"
        if os.path.exists(tmp_dst):
            try:
                os.remove(tmp_dst)
            except Exception:
                pass

        shutil.copy2(src, tmp_dst)
        # Ensure file permissions are reasonable (rw-r--r--)
        try:
            os.chmod(tmp_dst, 0o644)
        except Exception:
            # Not critical; may fail on some filesystems
            pass

        # Atomic move
        os.replace(tmp_dst, dst)
        log_crystal_progress(crystal_name, f"copied file {src} -> {dst} (atomic)")
        return True

    except PermissionError as pe:
        logging.error(f"Permission error copying {src} to {dst}: {pe}")
        # Fallback: try a direct copy without temp if allowed
        try:
            shutil.copy2(src, dst)
            log_crystal_progress(crystal_name, f"copied file {src} -> {dst} (fallback)")
            return True
        except Exception as e:
            logging.error(f"Fallback copy also failed for {src} -> {dst}: {e}")
            return False
    except FileNotFoundError:
        logging.warning(f"Source file not found: {src}")
        return False
    except Exception as e:
        logging.error(f"Error copying {src} to {dst}: {e}")
        return False

# Helper function to safely write text files
def safe_write_text_file(path, text, crystal_name):
    """Write text to file using atomic write strategy."""
    try:
        dst_dir = os.path.dirname(path)
        os.makedirs(dst_dir, exist_ok=True)
        tmp_path = path + ".tmp"
        with open(tmp_path, 'w', encoding='utf-8') as f:
            f.write(text)
        try:
            os.chmod(tmp_path, 0o644)
        except Exception:
            pass
        os.replace(tmp_path, path)
        log_crystal_progress(crystal_name, f"wrote text file {path} (atomic)")
        return True
    except PermissionError as pe:
        logging.error(f"Permission error writing {path}: {pe}")
        try:
            with open(path, 'w', encoding='utf-8') as f:
                f.write(text)
            log_crystal_progress(crystal_name, f"wrote text file {path} (fallback)")
            return True
        except Exception as e:
            logging.error(f"Fallback write failed for {path}: {e}")
            return False
    except Exception as e:
        logging.error(f"Error writing file {path}: {e}")
        return False

# Function to export restraints files to compound directory
def export_restraints_files(crystal, compound_code, entry, yaml_params, pipedream_dir):
    """Export CIF and PDB restraint files to the compound directory for all datasets.
    
    Args:
        crystal: Crystal name
        compound_code: Compound code
        entry: JSON entry with refinement information
        yaml_params: YAML parameters dictionary
        pipedream_dir: Path to Pipedream output directory
        
    Returns:
        tuple: (cif_path, pdb_path) of exported files, or (None, None) on failure
    """
    
    # Construct compound directory path
    compound_dir = os.path.join(yaml_params["Processing_directory"], "analysis", "model_building", crystal, "compound")
    crystal_dir = os.path.join(yaml_params["Processing_directory"], "analysis", "model_building", crystal)
    
    # Create compound directory if it doesn't exist
    os.makedirs(compound_dir, exist_ok=True)
    log_crystal_progress(crystal, f"created compound directory {compound_dir}")
    
    # Define target file paths in compound directory
    target_cif = os.path.join(compound_dir, f"{compound_code}.cif")
    target_pdb = os.path.join(compound_dir, f"{compound_code}.pdb")
    
    # Define symlink path for CIF in crystal directory
    symlink_cif = os.path.join(crystal_dir, f"{compound_code}.cif")
    
    try:
        # Get source files from Pipedream output
        output_pdb_file = entry.get("PDB File")  # This is the output PDB
        
        # Find output CIF and PDB files from rhofit directory
        # The CIF and PDB should be the refined restraints from Pipedream rhofit output
        output_cif_file = None
        refined_pdb_file = None
        
        # Primary location: rhofit-{CompoundCode} subdirectory in Pipedream output
        if pipedream_dir:
            rhofit_dir = os.path.join(pipedream_dir, f"rhofit-{compound_code}")
            primary_cif_path = os.path.join(rhofit_dir, "best.cif")
            primary_pdb_path = os.path.join(rhofit_dir, "best.pdb")
            
            if os.path.exists(primary_cif_path):
                output_cif_file = primary_cif_path
                log_crystal_progress(crystal, f"found CIF file at {primary_cif_path}")
            
            if os.path.exists(primary_pdb_path):
                refined_pdb_file = primary_pdb_path
                log_crystal_progress(crystal, f"found refined PDB file at {primary_pdb_path}")
        
        # Fallback: try other locations if rhofit files not found
        if not output_cif_file:
            # First, check if there's explicit CIF information in the entry
            if entry.get("Input Ligand Structure") and entry["Input Ligand Structure"] != "NA":
                input_cif_candidate = entry["Input Ligand Structure"].replace(".svg", ".cif").replace(".png", ".cif")
                if os.path.exists(input_cif_candidate):
                    output_cif_file = input_cif_candidate
                    log_crystal_progress(crystal, f"found CIF file from Input Ligand Structure: {output_cif_file}")
            
            # Try common fallback locations
            if not output_cif_file and output_pdb_file:
                pdb_dir = os.path.dirname(output_pdb_file)
                pipedream_base = os.path.dirname(os.path.dirname(pipedream_dir)) if pipedream_dir else ""
                
                possible_cif_paths = [
                    # In pipedream input_files (original restraints)
                    os.path.join(pipedream_dir, "input_files", f"{compound_code}.cif") if pipedream_dir else None,
                    # In original dataset input directory
                    os.path.join(pipedream_base, crystal, "input_files", f"{compound_code}.cif") if pipedream_base else None,
                    # In the same directory as output PDB
                    os.path.join(pdb_dir, f"{compound_code}.cif"),
                    # In restraints subdirectory
                    os.path.join(pdb_dir, "restraints", f"{compound_code}.cif")
                ]
                
                for cif_path in possible_cif_paths:
                    if cif_path and os.path.exists(cif_path):
                        output_cif_file = cif_path
                        log_crystal_progress(crystal, f"found fallback CIF file at {cif_path}")
                        break
        
        # Use refined PDB if available, otherwise fall back to original output PDB
        if not refined_pdb_file:
            refined_pdb_file = output_pdb_file
        
        if not output_cif_file:
            logging.warning(f"Could not find output CIF file for {compound_code} in {crystal}")
            log_crystal_progress(crystal, f"WARNING - no output CIF file found for {compound_code}")
        
        # Copy restraint files to compound directory
        cif_copied = False
        pdb_copied = False
        
        if refined_pdb_file and os.path.exists(refined_pdb_file):
            if safe_copy_file(refined_pdb_file, target_pdb, crystal):
                pdb_copied = True
            else:
                logging.warning(f"Refined PDB file copy failed: {refined_pdb_file} -> {target_pdb}")
        else:
            logging.warning(f"Refined PDB file not found: {refined_pdb_file}")

        if output_cif_file and os.path.exists(output_cif_file):
            if safe_copy_file(output_cif_file, target_cif, crystal):
                cif_copied = True
                # Create/update symlink to CIF in crystal directory
                safe_symlink(target_cif, symlink_cif)
                log_crystal_progress(crystal, f"updated CIF symlink {target_cif} -> {symlink_cif}")
            else:
                logging.warning(f"Failed to copy CIF file: {output_cif_file} -> {target_cif}")
        else:
            logging.warning(f"Output CIF file not found: {output_cif_file}")
        
        return (target_cif if cif_copied else None, target_pdb if pdb_copied else None)
        
    except Exception as e:
        logging.error(f"Error exporting restraints files for crystal {crystal}: {e}")
        log_crystal_progress(crystal, f"ERROR exporting restraints files: {e}")
        return (None, None)

# Function to handle chirality changes
def handle_chirality_change(crystal, compound_code, entry, yaml_params, pipedream_dir, new_smiles, cursor, user, db_timestamp):
    """Handle SMILES update, PNG regeneration, and database updates for chirality changes.
    Note: Restraints files are now exported separately for all datasets via export_restraints_files().
    """
    
    # Construct compound directory path
    compound_dir = os.path.join(yaml_params["Processing_directory"], "analysis", "model_building", crystal, "compound")
    crystal_dir = os.path.join(yaml_params["Processing_directory"], "analysis", "model_building", crystal)
    
    # Create compound directory if it doesn't exist
    os.makedirs(compound_dir, exist_ok=True)
    log_crystal_progress(crystal, f"ensured compound directory exists {compound_dir}")
    
    # Define target file paths in compound directory
    target_png = os.path.join(compound_dir, f"{compound_code}.png")
    target_smiles = os.path.join(compound_dir, f"{compound_code}.smiles")
    
    # Define symlink path for PNG in crystal directory
    symlink_png = os.path.join(crystal_dir, f"{compound_code}.png")
    
    try:
        # Handle SMILES update and PNG regeneration for chirality changes
        if new_smiles and new_smiles != "NA":
            if safe_write_text_file(target_smiles, new_smiles, crystal):
                log_crystal_progress(crystal, f"saved new SMILES to {target_smiles}")
            else:
                log_crystal_progress(crystal, f"failed to save SMILES to {target_smiles}")
            
            # Generate PNG from new SMILES
            if generate_png_from_smiles(new_smiles, target_png):
                log_crystal_progress(crystal, f"generated PNG diagram {target_png}")
                
                # Create/update symlink to PNG in crystal directory
                safe_symlink(target_png, symlink_png)
                log_crystal_progress(crystal, f"updated PNG symlink {target_png} -> {symlink_png}")
            else:
                log_crystal_progress(crystal, f"failed to generate PNG for {compound_code}")
        
        # Update SMILES in database
        if new_smiles and new_smiles != "NA":
            cursor.execute("UPDATE mainTable SET CompoundSMILES = ? WHERE CrystalName = ?", (new_smiles, crystal))
            log_crystal_progress(crystal, f"updated CompoundSMILES in database to: {new_smiles}")
        
        return True
        
    except Exception as e:
        logging.error(f"Error handling chirality change for crystal {crystal}: {e}")
        log_crystal_progress(crystal, f"ERROR handling chirality change: {e}")
        return False

# Main function
def main():
    parser = argparse.ArgumentParser(
        description="Export Pipedream results to XCE and update SQLite database.",
        epilog="Example: python export_pipedream_to_xce.py --input results.json --parameters parameters.yaml",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Version argument
    parser.add_argument(
        '--version', action='version', version=f'%(prog)s {VERSION}',
        help="Show version number"
    )
    
    # Required arguments
    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument(
        "--input", "--json", "-i",
        required=True,
        help="Path to the Pipedream results JSON file"
    )
    required_group.add_argument(
        "--parameters", "-p",
        required=True, 
        help="Path to the pipedream_parameters.yaml file"
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

    # Log input parameters to file
    logging.info(f"Starting export with parameters:")
    logging.info(f"  Input JSON file: {args.input}")
    logging.info(f"  Parameters file: {args.parameters}")
    logging.info(f"  Log file: {log_filename}")
    logging.info(f"  Log level: {args.log_level}")
    logging.info(f"  Verbose: {args.verbose}")

    # Load YAML parameters
    with open(args.parameters, "r") as f:
        yaml_params = yaml.safe_load(f)
    
    # Log YAML parameters to file
    log_to_file_only(f"YAML parameters: {yaml_params}")

    db_path = yaml_params["Database_path"]
    mode = yaml_params.get("Mode", "").strip('"').strip("'")
    output_dir = yaml_params.get("Output_directory", None)
    csv_path = yaml_params.get("Dataset_csv_path", None)  # Get CSV path from YAML
    
    # Log CSV path from YAML
    if csv_path:
        logging.info(f"  CSV file (from YAML): {csv_path}")
    if not output_dir:
        timestamp = datetime.now().strftime("Pipedream_%Y%m%d_%H%M%S")
        output_dir = os.path.join(yaml_params["Processing_directory"], "analysis", "Pipedream", timestamp)

    # Load dataset list if in specific_datasets mode
    dataset_list = None
    if mode == "specific_datasets":
        if csv_path:
            dataset_list = load_dataset_list(csv_path)
        else:
            raise ValueError("Dataset_csv_path must be specified in YAML file for specific_datasets mode")

    # Load JSON results
    with open(args.input, "r") as f:
        results = json.load(f)

    # Connect to the SQLite database
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    log_to_file_only(f"Connected to database: {db_path}")

    # Check required columns exist
    required_columns = [
        "RefinementResolution", "RefinementResolutionTL", "RefinementRcryst", "RefinementRcrystTraficLight", 
        "RefinementRfree", "RefinementRfreeTraficLight", "RefinementOutcome",
        "RefinementOutcomePerson", "RefinementOutcomeDate", "RefinementPDB_latest",
        "RefinementMTZ_latest", "RefinementMMCIFmodel_latest", "RefinementMMCIFreflections_latest",
        "RefinementLigandCC", "RefinementLigandConfidence", "RefinementBoundConformation", 
        "RefinementMolProbityScore", "RefinementMolProbityScoreTL",
        "RefinementRamachandranOutliers", "RefinementRamachandranOutliersTL", 
        "RefinementRamachandranFavored", "RefinementRamachandranFavoredTL",
        "RefinementStatus", "RefinementBusterReportHTML", "RefinementRefiner",
        "RefinementDate", "LastUpdated", "LastUpdated_by", "CrystalName",
        "DataProcessingResolutionHigh", "RefinementRmsdBonds", "RefinementRmsdBondsTL",
        "RefinementRmsdAngles", "RefinementRmsdAnglesTL", "CompoundSMILES", "RefinementCIFprogram",
        "RefinementComment"
    ]
    missing_columns = [col for col in required_columns if not column_exists(cursor, "mainTable", col)]
    if missing_columns:
        logging.warning(f"Missing columns in database table: {missing_columns}")
        log_to_file_only(f"Missing columns in database table: {missing_columns}")
    else:
        log_to_file_only("All required database columns are present")

    user = getpass.getuser()
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Special timestamp format for specific database columns: YYYY-MM-DD_hh-mm-ss.ms
    db_timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S.%f")[:-4]  # Remove last 4 digits to get milliseconds (ms)
    
    log_to_file_only(f"Processing user: {user}")
    log_to_file_only(f"Processing timestamp: {now}")
    log_to_file_only(f"Database timestamp format: {db_timestamp}")
    
    # Count total entries and entries to process
    total_entries = len(results)
    entries_to_process = len([e for e in results if e.get("Export to XCE") == "True"])
    log_to_file_only(f"Total entries in JSON: {total_entries}")
    log_to_file_only(f"Entries marked for export: {entries_to_process}")
    
    processed_count = 0
    error_count = 0

    for entry in results:
        if entry.get("Export to XCE") != "True":
            log_to_file_only(f"Skipping crystal {entry.get('Crystal Name', 'UNKNOWN')}: Export to XCE = {entry.get('Export to XCE')}")
            continue
        crystal = entry.get("Crystal Name")
        if dataset_list and crystal not in dataset_list:
            log_to_file_only(f"Skipping crystal {crystal}: not in dataset list")
            continue

        log_to_file_only(f"Starting processing for crystal: {crystal}")
        
        # Validate critical fields exist in JSON entry
        required_fields = ["Crystal Name", "Compound Code", "Pipedream Directory"]
        missing_fields = [field for field in required_fields if not entry.get(field)]
        if missing_fields:
            logging.error(f"Missing required fields for crystal {crystal}: {missing_fields}")
            log_crystal_progress(crystal, f"ERROR: Missing required fields: {missing_fields}")
            error_count += 1
            continue
        
        # Extract compound code once for this entry
        compound_code = entry.get("Compound Code", "UNK")
        
        # Check for chirality changes
        chirality_status = entry.get("Chiral inversion", "Not checked")
        has_chirality_change = (chirality_status != "stereochemistry matches input" and 
                              chirality_status != "Output stereochemistry matches input" and
                              chirality_status != "Not checked")
        
        if has_chirality_change:
            logging.info(f"Chirality change detected for crystal {crystal}: {chirality_status}")
            log_crystal_progress(crystal, f"chirality change detected - {chirality_status}")
            
            # Extract new SMILES from the {CompoundCode}_output_smiles.smiles file
            pipedream_dir = entry.get("Pipedream Directory", "")
            new_smiles = "NA"
            
            if pipedream_dir and compound_code:
                smiles_file = os.path.join(pipedream_dir, f"{compound_code}_output_smiles.smiles")
                log_crystal_progress(crystal, f"looking for SMILES file at {smiles_file}")
                
                if os.path.exists(smiles_file):
                    try:
                        with open(smiles_file, 'r') as f:
                            lines = f.readlines()
                            # Get the last non-empty line
                            for line in reversed(lines):
                                line = line.strip()
                                if line and not line.startswith('#'):
                                    new_smiles = line
                                    log_crystal_progress(crystal, f"extracted SMILES from last line: {new_smiles}")
                                    break
                        
                        if new_smiles == "NA":
                            log_crystal_progress(crystal, f"no valid SMILES found in file {smiles_file}")
                    except Exception as e:
                        log_crystal_progress(crystal, f"could not read SMILES file {smiles_file}: {e}")
                        new_smiles = "NA"
                else:
                    log_crystal_progress(crystal, f"SMILES file not found: {smiles_file}")
                    new_smiles = "NA"
            else:
                log_crystal_progress(crystal, f"missing Pipedream directory or compound code for SMILES extraction")
            
            # Handle chirality change (SMILES update, PNG regeneration, database update)
            if handle_chirality_change(crystal, compound_code, entry, yaml_params, entry.get("Pipedream Directory", ""), new_smiles, cursor, user, db_timestamp):
                log_crystal_progress(crystal, f"chirality change handling completed successfully")
            else:
                log_crystal_progress(crystal, f"chirality change handling failed")
        else:
            log_crystal_progress(crystal, f"no chirality change detected - {chirality_status}")

        try:
            pipedream_dir = entry["Pipedream Directory"]
            target_dir = os.path.join(yaml_params["Processing_directory"], "analysis", "model_building", crystal)
            
            # Export restraints files (CIF and PDB) to compound directory for ALL datasets
            logging.info(f"Exporting restraints files for crystal {crystal}")
            log_crystal_progress(crystal, f"exporting restraints files to compound directory")
            export_restraints_files(crystal, compound_code, entry, yaml_params, pipedream_dir)
            log_crystal_progress(crystal, f"restraints files export completed")
            os.makedirs(target_dir, exist_ok=True)
            
            log_crystal_progress(crystal, f"creating target directory {target_dir}")

            dest_pipedream = os.path.join(target_dir, os.path.basename(pipedream_dir))
            if not os.path.exists(dest_pipedream):
                log_crystal_progress(crystal, f"copying pipedream directory from {pipedream_dir} to {dest_pipedream}")
                shutil.copytree(pipedream_dir, dest_pipedream, symlinks=True, ignore_dangling_symlinks=True)
            else:
                log_crystal_progress(crystal, f"pipedream directory already exists at {dest_pipedream}")
            
            # Construct file paths directly from Pipedream directory structure
            # Standard Pipedream output locations in postrefine directory
            mtz_file_dest = os.path.join(dest_pipedream, f"postrefine-{compound_code}", "refine.mtz")
            postrefine_pdb = os.path.join(dest_pipedream, f"postrefine-{compound_code}", "refine.pdb")
            
            # Verify files exist and log warnings if they don't
            if not os.path.exists(mtz_file_dest):
                log_crystal_progress(crystal, f"WARNING: MTZ file not found at expected location: {mtz_file_dest}")
                logging.warning(f"MTZ file not found for crystal {crystal}: {mtz_file_dest}")
            
            if not os.path.exists(postrefine_pdb):
                log_crystal_progress(crystal, f"WARNING: postrefine PDB file not found at expected location: {postrefine_pdb}")
                logging.warning(f"Postrefine PDB not found for crystal {crystal}: {postrefine_pdb}")
            
            # Use postrefine PDB as main destination
            main_pdb_dest = postrefine_pdb
            log_crystal_progress(crystal, f"using postrefine files: PDB={postrefine_pdb}, MTZ={mtz_file_dest}")

            log_crystal_progress(crystal, f"creating symlinks")
            log_crystal_progress(crystal, f"  - {main_pdb_dest} -> {os.path.join(target_dir, 'refine.pdb')}")
            log_crystal_progress(crystal, f"  - {mtz_file_dest} -> {os.path.join(target_dir, 'refine.mtz')}")
            log_crystal_progress(crystal, f"  - {main_pdb_dest} -> {os.path.join(target_dir, 'refine.split.bound-state.pdb')}")

            safe_symlink(main_pdb_dest, os.path.join(target_dir, "refine.pdb"))
            safe_symlink(mtz_file_dest, os.path.join(target_dir, "refine.mtz"))
            safe_symlink(main_pdb_dest, os.path.join(target_dir, "refine.split.bound-state.pdb"))
            
            # Copy map files from the collation script results
            map_2fofc = entry.get("2Fo-Fc Map File")
            map_fofc = entry.get("Fo-Fc Map File")
            
            if map_2fofc and map_2fofc != "NA":
                map_2fofc_dest = map_2fofc.replace(pipedream_dir, dest_pipedream)
                if os.path.exists(map_2fofc_dest):
                    safe_symlink(map_2fofc_dest, os.path.join(target_dir, f"{crystal}_2fofc.map"))
                    log_crystal_progress(crystal, f"created 2Fo-Fc map symlink {map_2fofc_dest} -> {os.path.join(target_dir, f'{crystal}_2fofc.map')}")
                else:
                    log_crystal_progress(crystal, f"2Fo-Fc map file not found: {map_2fofc_dest}")
            else:
                log_crystal_progress(crystal, f"no 2Fo-Fc map file specified in JSON")
            
            if map_fofc and map_fofc != "NA":
                map_fofc_dest = map_fofc.replace(pipedream_dir, dest_pipedream)
                if os.path.exists(map_fofc_dest):
                    safe_symlink(map_fofc_dest, os.path.join(target_dir, f"{crystal}_fofc.map"))
                    log_crystal_progress(crystal, f"created Fo-Fc map symlink {map_fofc_dest} -> {os.path.join(target_dir, f'{crystal}_fofc.map')}")
                else:
                    log_crystal_progress(crystal, f"Fo-Fc map file not found: {map_fofc_dest}")
            else:
                log_crystal_progress(crystal, f"no Fo-Fc map file specified in JSON")

            # Construct summary path from Pipedream Directory
            summary_path = os.path.join(pipedream_dir, "pipedream_summary.json")
            with open(summary_path, "r") as f:
                summary = json.load(f)

            ligands = summary.get("ligandfitting", {}).get("ligands", [])
            
            # Process first ligand for molprobity and postref data (assuming same for all)
            if ligands:
                ligand = ligands[0]
                molprobity = ligand.get("validationstatistics", {}).get("molprobity", {})
                postref = summary.get("pipedream_outputs", {}).get("ligandfitting", {}).get("ligands", [{}])[0].get("postrefinement", [])
            else:
                molprobity = {}
                postref = []

            # Set mmCIF file paths to BUSTER files in the copied Pipedream directory
            # They should be in postrefine-{compound_code} subdirectory
            mmcif_model = os.path.join(dest_pipedream, f"postrefine-{compound_code}", "BUSTER_model.cif")
            mmcif_reflections = os.path.join(dest_pipedream, f"postrefine-{compound_code}", "BUSTER_refln.cif")
            
            # Check if BUSTER CIF files exist and log warnings if they don't
            if not os.path.exists(mmcif_model):
                log_crystal_progress(crystal, f"WARNING: BUSTER model CIF file not found at {mmcif_model}")
            if not os.path.exists(mmcif_reflections):
                log_crystal_progress(crystal, f"WARNING: BUSTER reflections CIF file not found at {mmcif_reflections}")

            rama_out = molprobity.get("ramaoutlierpercent")
            rama_fav = molprobity.get("ramafavoredpercent")
            molprob = molprobity.get("molprobityscore")
            r = entry.get("R")
            rfree = entry.get("Rfree")
            
            # Extract comment from JSON entry
            comment = entry.get("Comments", "")
            
            # Process multiple ligands for ligandcc
            ligandcc_entries = []
            overall_ligandcc = None
            
            for ligand in ligands:
                ligand_stats = ligand.get("validationstatistics", {}).get("ligandstatistics", [])
                for stat in ligand_stats:
                    ligand_id = stat.get("ligandid", "")
                    ligand_cc = stat.get("ligandcc", 0)
                    if ligand_id and ligand_cc:
                        ligandcc_entries.append(f"{ligand_id}: {ligand_cc}")
                        if overall_ligandcc is None or ligand_cc > overall_ligandcc:
                            overall_ligandcc = ligand_cc
            
            # Format ligandcc output
            if ligandcc_entries:
                ligandcc = "\n".join(ligandcc_entries)
            else:
                ligandcc = entry.get("Ligand CC")
                overall_ligandcc = ligandcc
            
            # Update report URL to use copied file path
            report_url = entry.get("Buster Report HTML")
            if report_url:
                report_url_dest = report_url.replace(pipedream_dir, dest_pipedream)
                report = report_url_dest
            else:
                report = None
            
            resolution = entry.get("High resolution (A)")
            
            # Extract RMSD values directly from the results JSON
            rmsd_bonds = entry.get("RMS bonds")
            rmsd_angles = entry.get("RMS angles")
            
            # Determine ligand confidence based on overall ligandcc value
            ligand_confidence = None
            refinement_outcome = "3 - In Refinement"  # Default value
            if overall_ligandcc is not None:
                try:
                    cc_value = float(overall_ligandcc)
                    if cc_value >= 0.8:
                        ligand_confidence = "4 - High Confidence"
                        refinement_outcome = "4 - CompChem ready"
                    elif cc_value >= 0.7:
                        ligand_confidence = "2 - Correct ligand, weak density"
                        refinement_outcome = "3 - In Refinement"
                    elif cc_value < 0.7:
                        ligand_confidence = "1 - Low Confidence"
                        refinement_outcome = "3 - In Refinement"
                except (ValueError, TypeError):
                    pass
            
            # Debug: Log the resolution value (file only)
            log_crystal_progress(crystal, f"resolution = {resolution}")
            log_crystal_progress(crystal, f"rmsd_bonds = {rmsd_bonds}, rmsd_angles = {rmsd_angles}")
            log_crystal_progress(crystal, f"r = {r}, rfree = {rfree}")
            log_crystal_progress(crystal, f"molprob = {molprob}, rama_out = {rama_out}, rama_fav = {rama_fav}")
            log_crystal_progress(crystal, f"overall_ligandcc = {overall_ligandcc}")
            log_crystal_progress(crystal, f"ligand_confidence = {ligand_confidence}, refinement_outcome = {refinement_outcome}")
            
            # Log all database columns to be updated
            log_crystal_progress(crystal, f"Database update values:")
            log_crystal_progress(crystal, f"  - RefinementResolution = {resolution}")
            log_crystal_progress(crystal, f"  - RefinementResolutionTL = {traffic_light(resolution, 2.0, 2.5) if resolution else None}")
            log_crystal_progress(crystal, f"  - RefinementOutcome = {refinement_outcome}")
            log_crystal_progress(crystal, f"  - RefinementRcryst = {r}")
            log_crystal_progress(crystal, f"  - RefinementRcrystTraficLight = {traffic_light(r, 0.20, 0.25)}")
            log_crystal_progress(crystal, f"  - RefinementRfree = {rfree}")
            log_crystal_progress(crystal, f"  - RefinementRfreeTraficLight = {traffic_light(rfree, 0.25, 0.30)}")
            log_crystal_progress(crystal, f"  - RefinementOutcomePerson = {user}")
            log_crystal_progress(crystal, f"  - RefinementOutcomeDate = {db_timestamp}")
            log_crystal_progress(crystal, f"  - RefinementPDB_latest = {main_pdb_dest}")
            log_crystal_progress(crystal, f"  - RefinementMTZ_latest = {mtz_file_dest}")
            log_crystal_progress(crystal, f"  - RefinementMMCIFmodel_latest = {mmcif_model}")
            log_crystal_progress(crystal, f"  - RefinementMMCIFreflections_latest = {mmcif_reflections}")
            log_crystal_progress(crystal, f"  - RefinementLigandCC = {ligandcc}")
            log_crystal_progress(crystal, f"  - RefinementLigandConfidence = {ligand_confidence}")
            log_crystal_progress(crystal, f"  - RefinementBoundConformation = {main_pdb_dest}")
            log_crystal_progress(crystal, f"  - RefinementMolProbityScore = {molprob}")
            log_crystal_progress(crystal, f"  - RefinementMolProbityScoreTL = {traffic_light(molprob, 2, 3)}")
            log_crystal_progress(crystal, f"  - RefinementRamachandranOutliers = {rama_out}")
            log_crystal_progress(crystal, f"  - RefinementRamachandranOutliersTL = {traffic_light(rama_out, 0.3, 1)}")
            log_crystal_progress(crystal, f"  - RefinementRamachandranFavored = {rama_fav}")
            log_crystal_progress(crystal, f"  - RefinementRamachandranFavoredTL = {traffic_light(rama_fav, 98, 95, reverse=True)}")
            log_crystal_progress(crystal, f"  - RefinementRmsdBonds = {rmsd_bonds}")
            log_crystal_progress(crystal, f"  - RefinementRmsdBondsTL = {traffic_light(rmsd_bonds, 0.012, 0.018)}")
            log_crystal_progress(crystal, f"  - RefinementRmsdAngles = {rmsd_angles}")
            log_crystal_progress(crystal, f"  - RefinementRmsdAnglesTL = {traffic_light(rmsd_angles, 1.5, 2.0)}")
            log_crystal_progress(crystal, f"  - RefinementStatus = 'finished'")
            log_crystal_progress(crystal, f"  - RefinementBusterReportHTML = {report}")
            log_crystal_progress(crystal, f"  - RefinementRefiner = {user}")
            log_crystal_progress(crystal, f"  - RefinementDate = {db_timestamp}")
            log_crystal_progress(crystal, f"  - LastUpdated = {db_timestamp}")
            log_crystal_progress(crystal, f"  - LastUpdated_by = {user}")
            log_crystal_progress(crystal, f"  - RefinementComment = {comment}")
            log_crystal_progress(crystal, f"  - WHERE CrystalName = {crystal}")
            
            update_query = f"""
            UPDATE mainTable SET
                RefinementResolution = ?,
                RefinementResolutionTL = ?,
                RefinementOutcome = ?,
                RefinementRcryst = ?,
                RefinementRcrystTraficLight = ?,
                RefinementRfree = ?,
                RefinementRfreeTraficLight = ?,
                RefinementOutcomePerson = ?,
                RefinementOutcomeDate = ?,
                RefinementPDB_latest = ?,
                RefinementMTZ_latest = ?,
                RefinementMMCIFmodel_latest = ?,
                RefinementMMCIFreflections_latest = ?,
                RefinementLigandCC = ?,
                RefinementLigandConfidence = ?,
                RefinementBoundConformation = ?,
                RefinementMolProbityScore = ?,
                RefinementMolProbityScoreTL = ?,
                RefinementRamachandranOutliers = ?,
                RefinementRamachandranOutliersTL = ?,
                RefinementRamachandranFavored = ?,
                RefinementRamachandranFavoredTL = ?,
                RefinementRmsdBonds = ?,
                RefinementRmsdBondsTL = ?,
                RefinementRmsdAngles = ?,
                RefinementRmsdAnglesTL = ?,
                RefinementStatus = 'finished',
                RefinementBusterReportHTML = ?,
                RefinementRefiner = ?,
                RefinementDate = ?,
                LastUpdated = ?,
                LastUpdated_by = ?,
                RefinementCIFprogram = 'Grade2',
                RefinementComment = ?
            WHERE CrystalName = ?
            """

            cursor.execute(update_query, (
                resolution,
                traffic_light(resolution, 2.0, 2.5) if resolution else None,
                refinement_outcome,
                r,
                traffic_light(r, 0.20, 0.25),
                rfree,
                traffic_light(rfree, 0.25, 0.30),
                user,
                db_timestamp,  # RefinementOutcomeDate
                main_pdb_dest,
                mtz_file_dest,
                mmcif_model,
                mmcif_reflections,
                ligandcc,
                ligand_confidence,
                main_pdb_dest,
                molprob,
                traffic_light(molprob, 2, 3),
                rama_out,
                traffic_light(rama_out, 0.3, 1,),
                rama_fav,
                traffic_light(rama_fav, 98, 95, reverse=True),
                rmsd_bonds,
                traffic_light(rmsd_bonds, 0.012, 0.018),
                rmsd_angles,
                traffic_light(rmsd_angles, 1.5, 2.0),
                report,
                user,
                db_timestamp,  # RefinementDate
                db_timestamp,  # LastUpdated
                user,
                comment,  # RefinementComment
                crystal
            ))

            conn.commit()
            log_crystal_progress(crystal, f"Database update completed successfully")
            logging.info(f"Updated database for crystal: {crystal}")
            processed_count += 1

        except Exception as e:
            logging.error(f"Error processing crystal {entry.get('Crystal Name', 'UNKNOWN')}: {e}")
            log_crystal_progress(crystal, f"ERROR processing: {e}")
            log_crystal_progress(crystal, f"Exception details: {str(e)}")
            error_count += 1

    conn.close()
    log_to_file_only(f"Database connection closed")
    log_to_file_only(f"Processing summary: {processed_count} successful, {error_count} errors")
    logging.info("Processing complete.")
    log_to_file_only(f"Export completed. Log saved to: {log_filename}")

if __name__ == "__main__":
    main()

