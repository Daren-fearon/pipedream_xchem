"""
Pipedream Post-Processing Script

This script generates electron density maps and runs edstats analysis on completed
Pipedream refinement outputs. Can be run standalone for backwards compatibility
or automatically as part of the Pipedream pipeline.

Features:
- Converts MTZ files to map files (2fofc and fofc) using gemmi
- Runs edstats to calculate ligand RSR values
- Saves metadata for downstream collation
- Can process single dataset or batch process multiple datasets
- Parallel processing support for batch mode

Author: DFearon
Date: November 2025
"""

import os
import json
import logging
import argparse
import subprocess
import sys
import shutil
from datetime import datetime
from pathlib import Path
from typing import Optional, Dict, List, Tuple, Any
from dataclasses import dataclass
import traceback

# Try to import tqdm for progress bar
try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False
    logging.debug("tqdm not available, progress bar will not be shown")

VERSION = "1.0.0"

@dataclass
class ProcessingConfig:
    """Configuration for processing."""
    gemmi_timeout_seconds: int = 30
    edstats_timeout_seconds: int = 60
    map_sample_rate: int = 5  # Sample rate for gemmi (ensures edstats compatibility)

CONFIG = ProcessingConfig()

# Software paths - rely on PATH like pipedream_xchem.py does
# Users should load appropriate modules (e.g., 'module load ccp4') before running
DEFAULT_GEMMI_PATHS = ['gemmi']
DEFAULT_EDSTATS_PATHS = ['edstats']


def setup_logging(log_file: str = None, log_level: str = "INFO", verbose: bool = False) -> None:
    """Configure logging."""
    handlers = []
    
    if log_file:
        log_dir = os.path.dirname(log_file)
        if log_dir:  # Only create directory if path contains a directory component
            os.makedirs(log_dir, exist_ok=True)
        handlers.append(logging.FileHandler(log_file))
    
    # Console handler
    console_handler = logging.StreamHandler()
    if verbose:
        console_handler.setLevel(getattr(logging, log_level))
    else:
        console_handler.setLevel(logging.WARNING)
    handlers.append(console_handler)
    
    logging.basicConfig(
        level=getattr(logging, log_level),
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=handlers,
        force=True
    )


def _load_module(module_name: str) -> bool:
    """
    Attempt to load an environment module (like 'module load').
    
    Args:
        module_name: Name of module to load (e.g., 'buster', 'ccp4')
        
    Returns:
        True if successful or not applicable, False if failed
    """
    try:
        # Check if running on a system with module command
        if shutil.which('modulecmd') is None:
            logging.debug(f"modulecmd not available, skipping module load {module_name}")
            return True  # Not an error, just not on a module system
        
        # Try to load the module using bash
        cmd = f"source /etc/profile 2>/dev/null && module load {module_name}"
        result = subprocess.run(
            ['bash', '-c', cmd],
            capture_output=True,
            text=True,
            timeout=10
        )
        
        if result.returncode == 0:
            logging.debug(f"Successfully loaded module: {module_name}")
            return True
        else:
            logging.debug(f"Could not load module {module_name}: {result.stderr}")
            return False
            
    except subprocess.TimeoutExpired:
        logging.debug(f"Timeout loading module: {module_name}")
        return False
    except Exception as e:
        logging.debug(f"Error loading module {module_name}: {e}")
        return False


def check_gemmi_availability(gemmi_paths: Optional[List] = None) -> Tuple[bool, Optional[List[str]]]:
    """Check if gemmi is available and return the command to use."""
    if gemmi_paths is None:
        gemmi_paths = DEFAULT_GEMMI_PATHS
    
    for gemmi_path in gemmi_paths:
        # Use shutil.which() to find command in PATH (like bash's 'command -v')
        found_path = shutil.which(str(gemmi_path))
        if found_path:
            logging.debug(f"gemmi found at: {found_path}")
            return True, [found_path]
    
    # Try loading buster module (contains gemmi)
    logging.info("gemmi not found in PATH, attempting to load buster module...")
    if _load_module('buster'):
        # Check again after loading module
        for gemmi_path in gemmi_paths:
            found_path = shutil.which(str(gemmi_path))
            if found_path:
                logging.info(f"gemmi found after loading buster module: {found_path}")
                return True, [found_path]
    
    logging.warning("gemmi not found in PATH")
    logging.info("Manually load buster module if gemmi is needed: module load buster")
    return False, None


def check_edstats_availability(edstats_paths: Optional[List] = None) -> Tuple[bool, Optional[List[str]]]:
    """Check if edstats is available and return the command to use."""
    if edstats_paths is None:
        edstats_paths = DEFAULT_EDSTATS_PATHS
    
    for edstats_path in edstats_paths:
        # Use shutil.which() to find command in PATH (like bash's 'command -v')
        found_path = shutil.which(str(edstats_path))
        if found_path:
            logging.debug(f"edstats found at: {found_path}")
            return True, [found_path]
    
    # Try loading ccp4 module (contains edstats)
    logging.info("edstats not found in PATH, attempting to load ccp4 module...")
    if _load_module('ccp4'):
        # Check again after loading module
        for edstats_path in edstats_paths:
            found_path = shutil.which(str(edstats_path))
            if found_path:
                logging.info(f"edstats found after loading ccp4 module: {found_path}")
                return True, [found_path]
    
    logging.warning("edstats not found in PATH")
    logging.info("Manually load CCP4 module if edstats is needed: module load ccp4")
    return False, None


def convert_mtz_to_map(mtz_file_path: str, map_type: str = '2fofc', 
                       gemmi_cmd: Optional[List[str]] = None) -> str:
    """
    Convert an MTZ file to a map file using gemmi with appropriate sampling for edstats.
    
    Args:
        mtz_file_path: Path to MTZ file
        map_type: '2fofc' or 'fofc'
        gemmi_cmd: gemmi command to use (if None, will auto-detect)
        
    Returns:
        Path to generated map file or 'NA' if failed
    """
    if not os.path.exists(mtz_file_path):
        logging.warning(f"MTZ file not found: {mtz_file_path}")
        return 'NA'
    
    map_file_path = mtz_file_path.replace('.mtz', f'_{map_type}.map')
    
    # Remove existing map file to ensure correct sampling
    if os.path.exists(map_file_path):
        try:
            os.remove(map_file_path)
            logging.debug(f"Removed existing map file: {map_file_path}")
        except (OSError, PermissionError) as e:
            logging.warning(f"Could not remove existing map file {map_file_path}: {e}")
            return map_file_path  # Try to use existing file
    
    # Auto-detect gemmi if not provided
    if gemmi_cmd is None:
        gemmi_available, gemmi_cmd = check_gemmi_availability()
        if not gemmi_available:
            logging.error("gemmi command not available")
            return 'NA'
    
    try:
        # Use --sample to ensure voxel_size < high_resolution/4 for edstats compatibility
        if map_type == '2fofc':
            cmd = gemmi_cmd + ['sf2map', '--sample', str(CONFIG.map_sample_rate), 
                              mtz_file_path, map_file_path]
        elif map_type == 'fofc':
            cmd = gemmi_cmd + ['sf2map', '--sample', str(CONFIG.map_sample_rate), 
                              '-d', mtz_file_path, map_file_path]
        else:
            raise ValueError(f"Unknown map type: {map_type}")
        
        logging.debug(f"Running gemmi command: {' '.join(str(c) for c in cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, 
                              timeout=CONFIG.gemmi_timeout_seconds, check=True)
        
        if os.path.exists(map_file_path):
            logging.info(f"Successfully created map file: {map_file_path}")
            return map_file_path
        else:
            logging.error(f"gemmi succeeded but map file not created: {map_file_path}")
            return 'NA'
            
    except subprocess.TimeoutExpired:
        logging.error(f"Timeout converting MTZ to map: {mtz_file_path}")
        return 'NA'
    except subprocess.CalledProcessError as e:
        logging.error(f"gemmi conversion failed: {e.stderr}")
        return 'NA'
    except Exception as e:
        logging.error(f"Error during map conversion: {e}")
        return 'NA'


def run_edstats(pdb_file: str, map_2fofc_file: str, map_fofc_file: str,
               output_file: str, low_resolution: Optional[float] = None,
               high_resolution: Optional[float] = None,
               edstats_cmd: Optional[List[str]] = None) -> bool:
    """
    Run edstats on PDB and map files.
    
    Args:
        pdb_file: Path to PDB file
        map_2fofc_file: Path to 2fofc map file
        map_fofc_file: Path to fofc map file
        output_file: Path to save edstats output
        low_resolution: Low resolution limit (optional)
        high_resolution: High resolution limit (optional)
        edstats_cmd: edstats command to use (if None, will auto-detect)
        
    Returns:
        True if successful, False otherwise
    """
    # Validate inputs
    if not os.path.exists(pdb_file):
        logging.error(f"PDB file not found: {pdb_file}")
        return False
    
    if not os.path.exists(map_2fofc_file):
        logging.warning(f"2fofc map file not found: {map_2fofc_file}")
        return False
        
    if not os.path.exists(map_fofc_file):
        logging.warning(f"fofc map file not found: {map_fofc_file}")
        return False
    
    # Auto-detect edstats if not provided
    if edstats_cmd is None:
        edstats_available, edstats_cmd = check_edstats_availability()
        if not edstats_available:
            logging.error("edstats command not available")
            return False
    
    # Check resolution limits are available (edstats requires them)
    if low_resolution is None or high_resolution is None:
        logging.error("Resolution limits are required for edstats but not available")
        return False
    
    try:
        # Prepare edstats input with resolution limits
        edstats_input = f"RESLO={low_resolution:.2f}\nRESHI={high_resolution:.2f}\nEND\n"
        logging.debug(f"Using resolution limits: RESLO={low_resolution:.2f}, RESHI={high_resolution:.2f}")
        
        # Build command using CCP4 style keywords
        cmd = edstats_cmd + [
            'XYZIN', pdb_file,
            'MAPIN1', map_2fofc_file,
            'MAPIN2', map_fofc_file,
            'OUT', output_file
        ]
        
        logging.debug(f"Running edstats command: {' '.join(str(c) for c in cmd)}")
        logging.debug(f"With input: {edstats_input.strip()}")
        
        result = subprocess.run(
            cmd,
            input=edstats_input,
            capture_output=True,
            text=True,
            timeout=CONFIG.edstats_timeout_seconds
        )
        
        if result.returncode != 0:
            logging.error(f"edstats failed with exit code {result.returncode}")
            if result.stderr:
                logging.error(f"edstats stderr: {result.stderr}")
            return False
        
        logging.info(f"Successfully ran edstats, output saved to: {output_file}")
        return True
        
    except subprocess.TimeoutExpired:
        logging.error(f"edstats timed out after {CONFIG.edstats_timeout_seconds} seconds")
        return False
    except Exception as e:
        logging.error(f"Error running edstats: {e}")
        return False


def find_pipedream_files(pipedream_dir: str, compound_code: str) -> Dict[str, str]:
    """
    Find required files in Pipedream output directory.
    
    Args:
        pipedream_dir: Pipedream output directory
        compound_code: Compound code
        
    Returns:
        Dictionary with file paths
    """
    files = {
        'mtz': 'NA',
        'pdb': 'NA',
        'postrefine_dir': 'NA'
    }
    
    # Try multiple possible directory structures
    possible_dirs = [
        # Standard structure: pipedream_dir/postrefine-{compound_code}
        os.path.join(pipedream_dir, f'postrefine-{compound_code}'),
        # Already in postrefine directory
        pipedream_dir,
        # Alternative: pipedream_dir might contain subdirectories
    ]
    
    # Also check for any postrefine-* directories in pipedream_dir
    if os.path.exists(pipedream_dir) and os.path.isdir(pipedream_dir):
        try:
            for item in os.listdir(pipedream_dir):
                item_path = os.path.join(pipedream_dir, item)
                if os.path.isdir(item_path) and item.startswith('postrefine-'):
                    possible_dirs.insert(1, item_path)  # Insert after first option
        except (PermissionError, OSError) as e:
            logging.debug(f"Could not list directory {pipedream_dir}: {e}")
    
    # Try each possible directory
    for postrefine_dir in possible_dirs:
        if not os.path.exists(postrefine_dir):
            continue
            
        # Look for MTZ and PDB files
        mtz_file = os.path.join(postrefine_dir, 'refine.mtz')
        pdb_file = os.path.join(postrefine_dir, 'refine.pdb')
        
        if os.path.exists(mtz_file) and os.path.exists(pdb_file):
            files['postrefine_dir'] = postrefine_dir
            files['mtz'] = mtz_file
            files['pdb'] = pdb_file
            logging.debug(f"Found files in: {postrefine_dir}")
            break
    
    return files


def extract_resolution_from_summary(pipedream_dir: str, compound_code: str) -> Tuple[Optional[float], Optional[float]]:
    """
    Extract resolution limits from Pipedream summary JSON file.
    
    Args:
        pipedream_dir: Pipedream output directory (base output directory, not postrefine subdirectory)
        compound_code: Compound code
        
    Returns:
        Tuple of (low_resolution, high_resolution) or (None, None) if failed
    """
    # Look for pipedream_summary.json in the base output directory
    summary_file = os.path.join(pipedream_dir, 'pipedream_summary.json')
    
    if not os.path.exists(summary_file):
        logging.debug(f"Summary file not found: {summary_file}")
        return None, None
    
    try:
        with open(summary_file, 'r', encoding='utf-8') as f:
            summary = json.load(f)
        
        # Extract resolution from dataprocessing.inputdata
        low_res = summary.get('dataprocessing', {}).get('inputdata', {}).get('reslo')
        high_res = summary.get('dataprocessing', {}).get('inputdata', {}).get('reshigh')
        
        if low_res is not None and high_res is not None:
            low_res = float(low_res)
            high_res = float(high_res)
            logging.debug(f"Extracted resolution from summary: low={low_res:.2f}, high={high_res:.2f}")
            return low_res, high_res
        else:
            logging.debug("Resolution not found in summary file")
            return None, None
            
    except Exception as e:
        logging.warning(f"Could not extract resolution from summary: {e}")
        return None, None


def process_single_dataset(pipedream_dir: str, compound_code: str, 
                          crystal_name: str = None,
                          edstats_path: str = None,
                          gemmi_path: str = None) -> Dict[str, Any]:
    """
    Process a single Pipedream output: generate maps and run edstats.
    
    Args:
        pipedream_dir: Pipedream output directory
        compound_code: Compound code
        crystal_name: Crystal name (optional, for logging)
        edstats_path: Custom path to edstats executable (optional)
        gemmi_path: Custom path to gemmi executable (optional)
        
    Returns:
        Dictionary with processing results
    """
    display_name = crystal_name or compound_code
    logging.info(f"Processing {display_name}")
    
    result = {
        'crystal_name': crystal_name,
        'compound_code': compound_code,
        'pipedream_dir': pipedream_dir,
        'success': False,
        'map_2fofc': 'NA',
        'map_fofc': 'NA',
        'edstats_output': 'NA',
        'timestamp': datetime.now().isoformat(),
        'errors': []
    }
    
    # Find required files
    files = find_pipedream_files(pipedream_dir, compound_code)
    
    if files['mtz'] == 'NA':
        error = f"MTZ file not found in {pipedream_dir}"
        logging.error(error)
        result['errors'].append(error)
        return result
    
    if files['pdb'] == 'NA':
        error = f"PDB file not found in {pipedream_dir}"
        logging.error(error)
        result['errors'].append(error)
        return result
    
    # Check gemmi availability
    if gemmi_path:
        gemmi_available, gemmi_cmd = check_gemmi_availability([gemmi_path])
    else:
        gemmi_available, gemmi_cmd = check_gemmi_availability()
    
    if not gemmi_available:
        error = "gemmi not available"
        logging.error(error)
        result['errors'].append(error)
        return result
    
    # Generate maps
    logging.info(f"Generating maps for {display_name}...")
    map_2fofc = convert_mtz_to_map(files['mtz'], '2fofc', gemmi_cmd)
    map_fofc = convert_mtz_to_map(files['mtz'], 'fofc', gemmi_cmd)
    
    result['map_2fofc'] = map_2fofc
    result['map_fofc'] = map_fofc
    
    if map_2fofc == 'NA' or map_fofc == 'NA':
        error = "Failed to generate maps"
        logging.warning(error)
        result['errors'].append(error)
        # Continue anyway - we can still try edstats if one map exists
    
    # Check edstats availability
    if edstats_path:
        edstats_available, edstats_cmd = check_edstats_availability([edstats_path])
    else:
        edstats_available, edstats_cmd = check_edstats_availability()
    
    if not edstats_available:
        error = "edstats not available"
        logging.warning(error)
        result['errors'].append(error)
        return result
    
    # Extract resolution limits from summary JSON
    logging.info(f"Extracting resolution limits for {display_name}...")
    low_res, high_res = extract_resolution_from_summary(pipedream_dir, compound_code)
    
    if low_res is None or high_res is None:
        logging.error(f"Could not extract resolution limits from summary JSON for {display_name} - edstats cannot run without resolution")
    
    # Run edstats
    if map_2fofc != 'NA' and map_fofc != 'NA':
        logging.info(f"Running edstats for {display_name}...")
        edstats_output = os.path.join(files['postrefine_dir'], 'edstats.out')
        
        edstats_success = run_edstats(
            files['pdb'],
            map_2fofc,
            map_fofc,
            edstats_output,
            low_res,
            high_res,
            edstats_cmd
        )
        
        if edstats_success:
            result['edstats_output'] = edstats_output
            result['success'] = True
            logging.info(f"Successfully processed {display_name}")
        else:
            error = "edstats failed"
            logging.error(error)
            result['errors'].append(error)
    else:
        error = "Cannot run edstats without both maps"
        logging.error(error)
        result['errors'].append(error)
    
    return result


def save_metadata(results: List[Dict[str, Any]], output_file: str) -> None:
    """Save processing results to JSON file."""
    try:
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2)
        logging.info(f"Saved metadata to: {output_file}")
    except Exception as e:
        logging.error(f"Failed to save metadata: {e}")


def process_batch(json_file: str, output_file: str = None, 
                  edstats_path: str = None, gemmi_path: str = None) -> None:
    """
    Process multiple datasets from JSON file (batch mode).
    
    Args:
        json_file: JSON file with dataset metadata from pipedream_xchem.py
        output_file: Output file for results (default: auto-generated)
        edstats_path: Custom path to edstats executable (optional)
        gemmi_path: Custom path to gemmi executable (optional)
    """
    logging.info(f"Processing batch from: {json_file}")
    
    # Read input JSON
    try:
        with open(json_file, 'r', encoding='utf-8') as f:
            data = json.load(f)
    except Exception as e:
        logging.error(f"Failed to read input JSON: {e}")
        return
    
    # Extract datasets
    datasets = []
    for crystal_name, metadata in data.items():
        if isinstance(metadata, dict) and 'PipedreamDirectory' in metadata:
            datasets.append({
                'crystal_name': crystal_name,
                'pipedream_dir': metadata['PipedreamDirectory'],
                'compound_code': metadata.get('CompoundCode', 'Unknown')
            })
    
    if not datasets:
        logging.error("No datasets found in JSON file")
        return
    
    logging.info(f"Found {len(datasets)} datasets to process")
    
    # Process datasets
    results = []
    
    if HAS_TQDM:
        iterator = tqdm(datasets, desc="Post-processing", unit="dataset")
    else:
        iterator = datasets
        print(f"Processing {len(datasets)} datasets...")
    
    for dataset in iterator:
        try:
            result = process_single_dataset(
                dataset['pipedream_dir'],
                dataset['compound_code'],
                dataset['crystal_name'],
                edstats_path=edstats_path,
                gemmi_path=gemmi_path
            )
            results.append(result)
        except Exception as e:
            logging.error(f"Error processing {dataset['crystal_name']}: {e}")
            results.append({
                'crystal_name': dataset['crystal_name'],
                'compound_code': dataset['compound_code'],
                'success': False,
                'errors': [str(e)]
            })
    
    # Save results
    if output_file is None:
        base_dir = os.path.dirname(json_file)
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        output_file = os.path.join(base_dir, f'post_process_results_{timestamp}.json')
    
    save_metadata(results, output_file)
    
    # Print summary
    successful = sum(1 for r in results if r.get('success', False))
    print(f"\n{'='*70}")
    print(f"Post-processing complete:")
    print(f"  Total datasets: {len(results)}")
    print(f"  Successful: {successful}")
    print(f"  Failed: {len(results) - successful}")
    print(f"  Results saved to: {output_file}")
    print(f"{'='*70}\n")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Post-process Pipedream outputs: generate maps and run edstats',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process single dataset
  python pipedream_post_process.py --pipedream-dir /path/to/output --compound-code LIG
  
  # Process batch from JSON
  python pipedream_post_process.py --json /path/to/Pipedream_output.json
        """
    )
    
    parser.add_argument('--version', action='version', version=f'%(prog)s {VERSION}')
    
    # Mode selection
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument(
        '--pipedream-dir',
        help='Pipedream output directory (single dataset mode)'
    )
    mode_group.add_argument(
        '--json', '--input', '-i',
        dest='input',
        help='JSON file with dataset metadata (batch mode)'
    )
    
    # Single dataset options
    parser.add_argument(
        '--compound-code',
        help='Compound code (required for single dataset mode)'
    )
    parser.add_argument(
        '--crystal-name',
        help='Crystal name (optional, for logging in single dataset mode)'
    )
    
    # Output options
    parser.add_argument(
        '--output',
        help='Output file for results (default: auto-generated with timestamp)'
    )
    
    # Logging options
    parser.add_argument(
        '--log-file',
        help='Log file path (default: post_process.log in current directory)'
    )
    parser.add_argument(
        '--log-level',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        default='INFO',
        help='Logging level (default: INFO)'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose console output'
    )
    
    # Software paths
    parser.add_argument(
        '--edstats-path',
        help='Path to edstats executable (default: auto-detect from PATH or CCP4 installations)'
    )
    parser.add_argument(
        '--gemmi-path',
        help='Path to gemmi executable (default: auto-detect from PATH or BUSTER installations)'
    )
    
    args = parser.parse_args()
    
    # Setup logging
    log_file = args.log_file or 'post_process.log'
    setup_logging(log_file, args.log_level, args.verbose)
    
    logging.info(f"Pipedream Post-Processing Script v{VERSION}")
    
    try:
        if args.input:
            # Batch mode
            process_batch(args.input, args.output, args.edstats_path, args.gemmi_path)
        else:
            # Single dataset mode
            if not args.compound_code:
                parser.error("--compound-code is required when using --pipedream-dir")
            
            result = process_single_dataset(
                args.pipedream_dir,
                args.compound_code,
                args.crystal_name,
                edstats_path=args.edstats_path,
                gemmi_path=args.gemmi_path
            )
            
            # Save result
            if args.output:
                save_metadata([result], args.output)
            
            # Print summary
            if result['success']:
                print(f"\n✓ Successfully processed {args.crystal_name or args.compound_code}")
                print(f"  Maps: {result['map_2fofc']}, {result['map_fofc']}")
                print(f"  Edstats: {result['edstats_output']}")
            else:
                print(f"\n✗ Failed to process {args.crystal_name or args.compound_code}")
                for error in result.get('errors', []):
                    print(f"  Error: {error}")
                sys.exit(1)
    
    except KeyboardInterrupt:
        print("\n\nOperation cancelled by user")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Fatal error: {e}")
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
