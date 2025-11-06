"""
Pipedream Results Collation and Analysis Script

This script collates and analyzes results from completed Pipedream crystallographic 
refinement runs. It processes multiple datasets, performs chirality analysis by 
comparing input and output ligand structures, extracts validation statistics, and 
generates comprehensive interactive HTML and JSON reports for quality assessment.

Features:
- Chirality inversion detection comparing input SMILES/PDB with output PDB structures
- Ligand validation statistics extraction (correlation coefficients, B-factors, etc.)
- MolProbity and structural quality metrics analysis
- SMILES extraction from refined ligand structures
- Interactive HTML reports with file links and data filtering
- Map file generation from MTZ files using gemmi

Author: DFearon
Date: October 2025
"""

import os
import json
import pandas as pd
import logging
import argparse
import sys
import traceback
from datetime import datetime
import webbrowser
import subprocess
from collections import OrderedDict
from typing import Any, Dict, List, Tuple, Optional, Protocol, TypedDict, Union, Literal
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D
import gemmi
from dataclasses import dataclass
from pathlib import Path
from contextlib import contextmanager
from enum import Enum, auto

# Version information
VERSION = "1.0.1"

# Constants
@dataclass
class ProcessingConfig:
    """Configuration for molecular processing."""
    min_ligand_atoms: int = 5
    mcs_timeout_seconds: int = 30
    min_mcs_atoms: int = 5
    max_debug_files_display: int = 5
    gemmi_timeout_seconds: int = 30

@dataclass
class OutputConfig:
    """Configuration for output generation."""
    default_page_length: int = 25
    scroll_height: str = '60vh'
    image_max_width: int = 150
    density_gif_width: int = 180
    density_gif_height: int = 120

PROCESSING = ProcessingConfig()
OUTPUT = OutputConfig()

# File name constants
class FileType(str, Enum):
    """File types used in Pipedream analysis."""
    BEST_PDB = "best.pdb"
    REFINE_PDB = "refine.pdb"
    BEST_CIF = "best.cif"
    REFINE_MTZ = "refine.mtz"
    DENSITY_GIF = "electrondensity_movie.gif"
    OUTPUT_SMILES = "_output_smiles.smiles"

BEST_PDB_FILENAME = FileType.BEST_PDB.value
REFINE_PDB_FILENAME = FileType.REFINE_PDB.value
BEST_CIF_FILENAME = FileType.BEST_CIF.value
REFINE_MTZ_FILENAME = FileType.REFINE_MTZ.value
ELECTRON_DENSITY_GIF_SUFFIX = FileType.DENSITY_GIF.value
OUTPUT_SMILES_SUFFIX = FileType.OUTPUT_SMILES.value

class SubdirPrefix(str, Enum):
    """Subdirectory prefixes in Pipedream output."""
    RHOFIT = "rhofit"
    POSTREFINE = "postrefine"
    REPORT = "report"

class MapType(str, Enum):
    """Map file types for electron density."""
    TWOFOFC = "2fofc"
    FOFC = "fofc"

def check_gemmi_availability():
    """
    Check if gemmi is available and return the command to use.
    
    Returns:
        Tuple[bool, Optional[List[str]]]: (is_available, command)
    """
    # Try direct gemmi command first
    try:
        result = subprocess.run(['gemmi', '--version'], 
                              capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            logging.debug(f"gemmi command available: {result.stdout.strip()}")
            return True, ['gemmi']
    except (FileNotFoundError, subprocess.TimeoutExpired, Exception) as e:
        logging.debug(f"Direct gemmi command not available: {e}")
    
    # Try python -m gemmi as fallback
    try:
        result = subprocess.run([sys.executable, '-m', 'gemmi', '--version'], 
                              capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            logging.debug(f"python -m gemmi available: {result.stdout.strip()}")
            return True, [sys.executable, '-m', 'gemmi']
    except Exception as e:
        logging.debug(f"python -m gemmi not available: {e}")
    
    # Try BUSTER gemmi path as final fallback
    buster_gemmi_path = '/dls_sw/apps/GPhL/BUSTER/20240123/autoBUSTER/bin/linux64/gemmi'
    try:
        result = subprocess.run([buster_gemmi_path, '--version'], 
                              capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            logging.debug(f"BUSTER gemmi available: {result.stdout.strip()}")
            return True, [buster_gemmi_path]
    except Exception as e:
        logging.debug(f"BUSTER gemmi not available: {e}")
    
    logging.debug("No gemmi installation found")
    return False, None

# CIF/PDB handling functions (from compare_chirality.py)
def load_mol_from_cif(cif_path: str) -> Chem.Mol:
    """Load a molecule from a CIF file with topology and coordinates."""
    if not os.path.exists(cif_path):
        raise FileNotFoundError(f"CIF file not found: {cif_path}")
    
    try:
        doc = gemmi.cif.read_file(cif_path)
    except (IOError, OSError) as e:
        raise ValueError(f"Cannot read CIF file {cif_path}: {e}")
    
    try:
        block = doc.find_block("comp_LIG") or doc.sole_block()
    except (RuntimeError, ValueError, AttributeError) as e:
        raise ValueError(f"Could not find valid CIF block in {cif_path}: {e}")
    
    BOND_TYPES = {
        'SING': Chem.BondType.SINGLE, 'SINGLE': Chem.BondType.SINGLE,
        'DOUB': Chem.BondType.DOUBLE, 'DOUBLE': Chem.BondType.DOUBLE,
        'TRIP': Chem.BondType.TRIPLE, 'TRIPLE': Chem.BondType.TRIPLE,
        'AROM': Chem.BondType.AROMATIC, 'AROMATIC': Chem.BondType.AROMATIC,
    }

    def strip_quotes(s): return s.strip('"').strip("'")

    mol = Chem.RWMol()
    conf = Chem.Conformer()

    comp_ids = block.find_loop("_chem_comp_atom.comp_id")
    atom_ids = block.find_loop("_chem_comp_atom.atom_id")
    atom_symbols = block.find_loop("_chem_comp_atom.type_symbol")
    x = block.find_loop("_chem_comp_atom.x") or block.find_loop("_chem_comp_atom.model_Cartn_x")
    y = block.find_loop("_chem_comp_atom.y") or block.find_loop("_chem_comp_atom.model_Cartn_y")
    z = block.find_loop("_chem_comp_atom.z") or block.find_loop("_chem_comp_atom.model_Cartn_z")
    charges = [0] * len(atom_ids)
    if block.find_loop("_chem_comp_atom.charge"):
        charges = list(block.find_loop("_chem_comp_atom.charge"))
    elif block.find_loop("_chem_comp_atom.partial_charge"):
        charges = list(block.find_loop("_chem_comp_atom.partial_charge"))

    atoms = {}
    for name, s, id, px, py, pz, charge in zip(comp_ids, atom_symbols, atom_ids, x, y, z, charges):
        id = strip_quotes(id)
        if len(s) == 2: s = s[0] + s[1].lower()
        atom = Chem.Atom(s)
        atom.SetFormalCharge(round(float(charge)))
        atom.SetProp("atom_id", id)
        idx = mol.AddAtom(atom)
        atom.SetIntProp("idx", idx)
        atoms[id] = atom
        conf.SetAtomPosition(idx, Point3D(float(px), float(py), float(pz)))

    atom1 = block.find_loop("_chem_comp_bond.atom_id_1")
    atom2 = block.find_loop("_chem_comp_bond.atom_id_2")
    bond_type = block.find_loop("_chem_comp_bond.type") or block.find_loop("_chem_comp_bond.value_order")

    for a1, a2, bt in zip(atom1, atom2, bond_type):
        mol.AddBond(atoms[strip_quotes(a1)].GetIntProp("idx"),
                    atoms[strip_quotes(a2)].GetIntProp("idx"),
                    BOND_TYPES[bt.upper()])

    mol.AddConformer(conf)
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistryFrom3D(mol)
    return mol

def assign_coordinates_from_pdb(mol: Chem.Mol, pdb_path: str, dataset_name: str = "Unknown", structure_type: str = "Unknown") -> Chem.Mol:
    """Assign 3D coordinates from a PDB file to a molecule from CIF."""
    with open(pdb_path, 'r') as f:
        pdb_block = f.read()
    pdb_mol = Chem.MolFromPDBBlock(pdb_block, removeHs=False)
    if pdb_mol is None:
        raise ValueError("Failed to parse PDB file.")
    
    mol_with_h = Chem.AddHs(mol)
    
    # Try constrained embedding first
    try:
        AllChem.ConstrainedEmbed(mol_with_h, pdb_mol)
    except (ValueError, RuntimeError) as e:
        # ConstrainedEmbed failed - this is normal when CIF and PDB have 
        # structural differences (e.g., after stereochemical changes)
        logging.debug(f"Dataset {dataset_name}: ConstrainedEmbed failed for {structure_type}, using manual coordinate assignment: {e}")
        conf = Chem.Conformer(mol_with_h.GetNumAtoms())
        
        # Get the conformers safely
        try:
            pdb_conf = pdb_mol.GetConformer()
            mol_conf = mol_with_h.GetConformer()
        except (ValueError, RuntimeError) as e:
            # If there's no conformer, we can't proceed with coordinate assignment
            logging.debug(f"Dataset {dataset_name}: No valid conformer available for {structure_type} coordinate assignment from PDB, using original CIF coordinates")
            Chem.AssignStereochemistryFrom3D(mol_with_h)
            return mol_with_h
        
        # Assign coordinates by atom index
        for i in range(mol_with_h.GetNumAtoms()):
            if i < pdb_mol.GetNumAtoms():
                try:
                    pos = pdb_conf.GetAtomPosition(i)
                    conf.SetAtomPosition(i, pos)
                except (ValueError, IndexError, RuntimeError) as e:
                    # If we can't get position from PDB, use CIF position
                    try:
                        pos = mol_conf.GetAtomPosition(i)
                        conf.SetAtomPosition(i, pos)
                    except (ValueError, IndexError, RuntimeError):
                        # Last resort: set to origin
                        conf.SetAtomPosition(i, Point3D(0.0, 0.0, 0.0))
            else:
                # Keep original CIF position for atoms not in PDB
                try:
                    pos = mol_conf.GetAtomPosition(i)
                    conf.SetAtomPosition(i, pos)
                except (ValueError, IndexError, RuntimeError):
                    # Last resort: set to origin
                    conf.SetAtomPosition(i, Point3D(0.0, 0.0, 0.0))
        
        # Remove existing conformer and add the new one
        try:
            mol_with_h.RemoveConformer(0)
        except (ValueError, RuntimeError):
            # No conformer to remove
            pass
        mol_with_h.AddConformer(conf)
    
    Chem.AssignStereochemistryFrom3D(mol_with_h)
    return mol_with_h

def compare_chiral_centers_cif_pdb(input_cif, input_pdb, output_cif, output_pdb, smiles_file, dataset_name="Unknown"):
    """Compare stereochemistry using CIF files for topology and PDB files for coordinates."""
    try:
        # Load reference SMILES
        with open(smiles_file, 'r') as f:
            smiles = f.readline().strip()

        # Process input molecule
        logging.debug(f"Loading input CIF: {input_cif}")
        mol_input = load_mol_from_cif(input_cif)
        
        # Assign coordinates from input PDB if available, otherwise use CIF coordinates
        if input_pdb and input_pdb != 'NA' and os.path.isfile(input_pdb):
            logging.debug(f"Assigning coordinates from input PDB: {input_pdb}")
            mol_input = assign_coordinates_from_pdb(mol_input, input_pdb, dataset_name, "input")
        else:
            logging.debug(f"No input PDB available, using CIF coordinates for input molecule")
            # Ensure stereochemistry is assigned from the CIF 3D coordinates
            Chem.AssignStereochemistryFrom3D(mol_input)

        # Process output molecule - use CIF topology + PDB coordinates
        logging.debug(f"Loading output CIF: {output_cif}")
        mol_output = load_mol_from_cif(output_cif)
        logging.debug(f"Assigning coordinates from output PDB: {output_pdb}")
        mol_output = assign_coordinates_from_pdb(mol_output, output_pdb, dataset_name, "output")

        # Find chiral centers
        chiral_input = Chem.FindMolChiralCenters(mol_input, includeUnassigned=True)
        chiral_output = Chem.FindMolChiralCenters(mol_output, includeUnassigned=True)

        # Generate SMILES without explicit hydrogens - handle potential errors
        try:
            mol_input_no_h = Chem.RemoveHs(mol_input)
            smiles_input = Chem.MolToSmiles(mol_input_no_h, isomericSmiles=True)
        except Exception as e:
            logging.warning(f"Failed to generate input SMILES: {e}")
            smiles_input = "Failed to generate"
            
        try:
            mol_output_no_h = Chem.RemoveHs(mol_output)
            smiles_output = Chem.MolToSmiles(mol_output_no_h, isomericSmiles=True)
        except Exception as e:
            logging.warning(f"Failed to generate output SMILES: {e}")
            smiles_output = "Failed to generate"

        # Compare chiral centers
        differences = []
        input_dict = {idx: config for idx, config in chiral_input}
        output_dict = {idx: config for idx, config in chiral_output}
        
        all_chiral_indices = set(input_dict.keys()) | set(output_dict.keys())
        
        for idx in all_chiral_indices:
            input_config = input_dict.get(idx, 'Not found')
            output_config = output_dict.get(idx, 'Not found')
            if input_config != output_config:
                differences.append(f"Atom {idx}: {input_config} → {output_config}")

        if differences:
            result = f"Chiral centre inverted:\n{'; '.join(differences)}"
        else:
            result = "Output stereochemistry matches input"
            
        logging.debug(f"Input chiral centers: {chiral_input}")
        logging.debug(f"Output chiral centers: {chiral_output}")
        logging.debug(f"Reference SMILES: {smiles}")
        logging.debug(f"Input SMILES: {smiles_input}")
        logging.debug(f"Output SMILES: {smiles_output}")
        
        return result, smiles_output
        
    except FileNotFoundError as e:
        logging.error(f"File not found in CIF/PDB chirality comparison: {e}")
        return f"File not found: {e}", "NA"
    except (IOError, OSError) as e:
        logging.error(f"File I/O error in CIF/PDB chirality comparison: {e}")
        return f"File error: {e}", "NA"
    except ValueError as e:
        logging.error(f"Value error in CIF/PDB chirality comparison: {e}")
        logging.debug(f"Full traceback: {traceback.format_exc()}")
        return f"Structure error: {e}", "NA"
    except Exception as e:
        # Keep generic catch for truly unexpected errors
        logging.error(f"Unexpected error in CIF/PDB chirality comparison: {e}")
        logging.debug(f"Full traceback: {traceback.format_exc()}")
        return f"Error: {str(e)}", "NA"

# Column order for consistent output across HTML and JSON
COLUMN_ORDER = [
    'Export to XCE',
    'Comments',
    'Crystal Name',
    'Compound Code',
    'Input Ligand Structure',
    'Ligand Density',
    'Chiral inversion',
    'Pipedream Directory',
    'Buster Report HTML',
    'Ligand Report HTML',
    'Ligand CC',
    'Ligand occupancy',
    'Ligand avg B factor',
    'Mean B factor',
    'High resolution (A)',
    'R',
    'Rfree',
    'Mogul Z angle',
    'Mogul Z bond',
    'c beta deviations',
    'Rama outlier percent',
    'Rama favored percent',
    'Poor rotamer percent',
    'Clash score',
    'Mol probity score',
    'RMS bonds',
    'RMS angles',
    'Pipedream Summary',
    'PDB File',
    'MTZ File',
    '2Fo-Fc Map File',
    'Fo-Fc Map File',
    'Output SMILES'
]

def _extract_ligand_from_pdb(mol: Chem.Mol) -> Chem.Mol:
    """Extract ligand fragment from a PDB molecule that may contain protein + ligand."""
    try:
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
        if len(frags) > 1:
            # Sort fragments by size and take the smallest (likely the ligand)
            frags_with_size = [(frag, frag.GetNumAtoms()) for frag in frags]
            frags_with_size.sort(key=lambda x: x[1])
            
            # Take the smallest fragment that's not too small (at least 5 atoms)
            for frag, size in frags_with_size:
                if size >= PROCESSING.min_ligand_atoms:
                    logging.debug(f"Extracted ligand fragment with {size} atoms from {len(frags)} fragments")
                    return frag
            
            logging.warning("Could not identify ligand fragment, using full structure")
        else:
            logging.debug(f"Output PDB contains single molecule with {mol.GetNumAtoms()} atoms")
    except (ValueError, RuntimeError, AttributeError) as e:
        logging.warning(f"Error extracting ligand from PDB, using full structure: {e}")
    
    return mol

def _save_smiles_file(pipedream_dir: str, compound_code: str, 
                     chirality_flip: str, updated_smiles: str, 
                     method: str, dataset_name: str) -> str:
    """Save SMILES analysis file and return its path."""
    smiles_file = "NA"
    
    # Only save if we have valid paths and meaningful SMILES
    if pipedream_dir != 'NA' and compound_code != 'NA' and updated_smiles not in ['NA', 'Failed to generate SMILES']:
        try:
            output_path = get_output_smiles_path(pipedream_dir, compound_code)
            
            with open(output_path, "w", encoding="utf-8") as f:
                f.write(f"# Output SMILES analysis for compound {compound_code}\n")
                f.write(f"# Dataset: {dataset_name}\n")
                f.write(f"# Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"# Analysis method: {method}\n")
                f.write(f"# Chirality analysis: {chirality_flip}\n")
                f.write(f"# Output SMILES: {updated_smiles}\n")
                f.write(f"#\n")
                f.write(f"# SMILES string:\n")
                f.write(f"{updated_smiles}\n")
            
            smiles_file = output_path
            logging.debug(f"Saved output SMILES analysis to: {output_path}")
            
        except (IOError, OSError, PermissionError) as e:
            logging.error(f"Error saving output SMILES for {compound_code}: {e}")
            smiles_file = "NA"
    else:
        # Construct expected path even if not saving
        if pipedream_dir != 'NA' and compound_code != 'NA':
            expected_path = get_output_smiles_path(pipedream_dir, compound_code)
            if os.path.exists(expected_path):
                smiles_file = expected_path
    
    return smiles_file

def save_output_smiles_with_analysis(output_mol, compound_code, pipedream_dir, chirality_result="Not analyzed", input_smiles=""):
    """Save output SMILES with chirality analysis results to the specific pipedream run directory."""
    try:
        # Remove hydrogens and assign stereochemistry
        mol_no_h = Chem.RemoveHs(output_mol)
        Chem.AssignStereochemistry(mol_no_h, force=True, cleanIt=True)
        # Generate canonical isomeric SMILES
        output_smiles = Chem.MolToSmiles(mol_no_h, isomericSmiles=True, canonical=True)
        
        # Create output path using helper function
        output_path = get_output_smiles_path(pipedream_dir, compound_code)
            
        # Write file with comprehensive metadata
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(f"# Output SMILES analysis for compound {compound_code}\n")
            f.write(f"# Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"# Chirality analysis: {chirality_result}\n")
            if input_smiles:
                f.write(f"# Input SMILES:  {input_smiles}\n")
            f.write(f"# Output SMILES: {output_smiles}\n")
            f.write(f"#\n")
            f.write(f"# SMILES string:\n")
            f.write(f"{output_smiles}\n")
            
        logging.debug(f"Saved output SMILES analysis to: {output_path}")
        return output_path
    except (ValueError, RuntimeError, AttributeError) as e:
        logging.error(f"RDKit error saving output SMILES for {compound_code}: {e}")
        return f"Error: RDKit processing failed - {e}"
    except (IOError, OSError, PermissionError) as e:
        logging.error(f"File I/O error saving output SMILES for {compound_code}: {e}")
        return f"Error: Cannot write file - {e}"

def detect_chiral_inversion(input_smi_path: str, output_pdb_path: str, input_pdb_path: Optional[str] = None, pipedream_dir: str = 'NA', compound_code: str = 'NA') -> Tuple[str, str]:
    """
    Detect chirality differences between input and output structures.
    Uses input SMILES (or optionally input PDB) and output PDB files.
    
    Args:
        input_smi_path: Path to input SMILES file
        output_pdb_path: Path to output PDB file
        input_pdb_path: Optional path to input PDB file (preferred over SMILES)
        pipedream_dir: Pipedream directory path for saving canonical SMILES output files
        compound_code: Compound code used for naming output SMILES files
    
    Returns:
        Tuple of (chirality comparison result, updated SMILES from output structure)
    """
    # File existence checks
    if not os.path.isfile(input_smi_path):
        return f"Input SMILES not found: {input_smi_path}", "NA"
    if not os.path.isfile(output_pdb_path):
        return f"Output PDB not found: {output_pdb_path}", "NA"

    try:
        # Try to use input PDB if available, otherwise fall back to SMILES
        if input_pdb_path and os.path.isfile(input_pdb_path):
            logging.debug(f"Using input PDB for comparison: {input_pdb_path}")
            input_mol = Chem.MolFromPDBFile(input_pdb_path, removeHs=False)
            if input_mol is None:
                logging.warning(f"Failed to parse input PDB, falling back to SMILES")
                input_mol = None
            else:
                Chem.AssignStereochemistry(input_mol, force=True, cleanIt=True)
                source = "PDB"
        
        # Fall back to SMILES if PDB not available or failed
        if input_mol is None:
            logging.debug(f"Using input SMILES for comparison: {input_smi_path}")
            with open(input_smi_path, 'r', encoding='utf-8') as f:
                smiles_line = f.readline().strip()
            smiles = smiles_line.split()[0] if smiles_line else ""
            if not smiles:
                return f"Empty or invalid SMILES file: {input_smi_path}", "NA"

            logging.debug(f"Processing SMILES: {smiles}")
            input_mol = Chem.MolFromSmiles(smiles)
            if input_mol is None:
                return f"Failed to parse SMILES: {smiles}", "NA"
            source = "SMILES"

        # Read output structure from PDB file
        logging.debug(f"Reading output structure from PDB file: {output_pdb_path}")
        output_mol = Chem.MolFromPDBFile(output_pdb_path, removeHs=False)
        if output_mol is None:
            return f"Failed to parse output PDB file: {output_pdb_path}", "NA"

        # Extract only the ligand from the output PDB
        output_mol = _extract_ligand_from_pdb(output_mol)
        
        # CRITICAL: Use input SMILES as template to preserve bond orders
        # PDB files lack bond order information, so we use the input structure as a template
        template_mol = None
        try:
            # Get the input SMILES string to use as template
            if source == "SMILES":
                # We already have the input molecule from SMILES
                template_mol = input_mol
                logging.debug("Using input SMILES molecule as template for bond order assignment")
            else:
                # If input was from PDB, try to get SMILES from the original SMILES file
                logging.debug(f"Reading SMILES template from: {input_smi_path}")
                with open(input_smi_path, 'r', encoding='utf-8') as f:
                    smiles_line = f.readline().strip()
                template_smiles = smiles_line.split()[0] if smiles_line else ""
                if template_smiles:
                    template_mol = Chem.MolFromSmiles(template_smiles)
                    logging.debug(f"Using SMILES template: {template_smiles}")
                
            if template_mol is not None:
                # Use template to assign correct bond orders to output molecule
                
                # Remove hydrogens from both molecules for template matching
                template_no_h = Chem.RemoveHs(template_mol)
                output_no_h = Chem.RemoveHs(output_mol)
                
                # Assign bond orders from template - this preserves the correct connectivity
                output_mol_corrected = AllChem.AssignBondOrdersFromTemplate(template_no_h, output_no_h)
                
                if output_mol_corrected is not None:
                    output_mol = output_mol_corrected
                    logging.debug("Successfully assigned bond orders using input SMILES template")
                else:
                    logging.warning("Failed to assign bond orders from template, using standard sanitization")
                    raise Exception("Template assignment failed")
            else:
                raise Exception("No template molecule available")
                
        except Exception as e:
            logging.warning(f"Could not use template for bond order assignment: {e}, falling back to standard sanitization")
            # Fallback to standard sanitization if template approach fails
            try:
                Chem.SanitizeMol(output_mol)
            except Exception as e2:
                logging.warning(f"Standard sanitization failed: {e2}, trying alternative approach")
                try:
                    output_mol = Chem.AddHs(output_mol)
                    Chem.SanitizeMol(output_mol)
                    output_mol = Chem.RemoveHs(output_mol)
                except Exception as e3:
                    logging.warning(f"Alternative sanitization also failed: {e3}")
        
        # Now assign stereochemistry with correct bond orders
        Chem.AssignStereochemistry(output_mol, force=True, cleanIt=True)
        
        # Generate SMILES from output ligand (exclude hydrogens)
        output_mol_no_h = Chem.RemoveHs(output_mol)
        
        # Get input SMILES for comparison metadata
        input_smiles_for_metadata = ""
        if source == "SMILES":
            input_smiles_for_metadata = smiles
        else:
            # Try to read from SMILES file for metadata
            try:
                with open(input_smi_path, 'r', encoding='utf-8') as f:
                    input_smiles_for_metadata = f.readline().strip().split()[0]
            except (IOError, IndexError, OSError):
                pass
        
        updated_smiles = Chem.MolToSmiles(output_mol_no_h, isomericSmiles=True, canonical=True)
        if not updated_smiles:
            updated_smiles = "Failed to generate SMILES"

        # Find chiral centers
        input_chiral_centers = _get_chiral_centers(input_mol)
        output_chiral_centers = _get_chiral_centers(output_mol)
        
        logging.debug(f"Input chiral centers: {len(input_chiral_centers)}")
        logging.debug(f"Output chiral centers: {len(output_chiral_centers)}")

        # Additional molecular information
        input_atoms = input_mol.GetNumAtoms()
        output_atoms = output_mol.GetNumAtoms();
        
        logging.debug(f"Input ligand ({source}): {input_atoms} atoms")
        logging.debug(f"Output ligand (PDB): {output_atoms} atoms")
        
        # For both PDB and SMILES inputs, we need to handle atom mapping since 
        # atom indices don't correspond between input and output structures        
        try:
            # Align molecules and map chiral centers
            mapped_output_chiral = _align_molecules_and_map_chiral_centers(
                input_mol, output_mol, input_chiral_centers, output_chiral_centers, source
            )
            
            # Compare mapped chirality
            input_dict = {idx: config for idx, config in input_chiral_centers}
            output_dict = mapped_output_chiral;
                
        except Exception as e:
            return f"Error during structural alignment: {e}", updated_smiles

        # Check for differences in chiral centers
        differences = []
        for atom_idx in input_dict:
            if atom_idx in output_dict:
                if output_dict[atom_idx] is None:
                    differences.append(f"Atom {atom_idx}: {input_dict[atom_idx]} -> missing")
                elif input_dict[atom_idx] != output_dict[atom_idx]:
                    differences.append(f"Atom {atom_idx}: {input_dict[atom_idx]} -> {output_dict[atom_idx]}")
            else:
                differences.append(f"Atom {atom_idx}: {input_dict[atom_idx]} -> missing")

        if differences:
            result = "Chiral centre inverted: " + "; ".join(differences)
            logging.debug(f"Chirality comparison complete: {result}")
            logging.debug(f"Input source: {source}, Output source: PDB")
        else:
            result = "Stereochemistry matches input"
            logging.debug(f"Chirality comparison complete: {result}")
            logging.debug(f"Input source: {source}, Output source: PDB")
            # If no differences, don't provide updated SMILES
            updated_smiles = "NA"
        
        # Save comprehensive SMILES analysis file
        smiles_file_path = "NA"
        if pipedream_dir != 'NA' and compound_code != 'NA':
            smiles_file_path = save_output_smiles_with_analysis(
                output_mol, compound_code, pipedream_dir, result, input_smiles_for_metadata
            )
            if not smiles_file_path.startswith("Error"):
                logging.debug(f"Output SMILES analysis saved to: {smiles_file_path}")
        
        return result, updated_smiles

    except (IOError, OSError) as e:
        logging.error(f"File I/O error during chirality comparison for {input_smi_path} -> {output_pdb_path}: {e}")
        return f"File error: {e}", "NA"
    except (ValueError, RuntimeError) as e:
        logging.error(f"RDKit error during chirality comparison for {input_smi_path} -> {output_pdb_path}: {e}")
        return f"Structure error: {e}", "NA"
    except Exception as e:
        # Keep a generic catch for truly unexpected errors, but log them prominently
        logging.error(f"Unexpected error during chirality comparison for {input_smi_path} -> {output_pdb_path}: {e}")
        logging.error(f"Traceback: {traceback.format_exc()}")
        return f"Unexpected error: {e}", "NA"


def _get_chiral_centers(mol: Chem.Mol) -> List[Tuple[int, str]]:
    """
    Get all chiral centers from a molecule, including genuinely unassigned ones.
    
    Args:
        mol: RDKit molecule object
        
    Returns:
        List of (atom_index, configuration) tuples
    """
    # Get assigned chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
    
    # Get all potential chiral centers (including unassigned)
    all_chiral = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    
    # Add genuinely unassigned chiral centers (exclude false positives)
    assigned_indices = {idx for idx, _ in chiral_centers}
    
    for idx, config in all_chiral:
        if config == '?' and idx not in assigned_indices:
            # Verify this is actually a chiral center
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetDegree() == 4 and len(set(n.GetSymbol() for n in atom.GetNeighbors())) >= 3:
                chiral_centers.append((idx, config))
    
    return chiral_centers


def _align_molecules_and_map_chiral_centers(input_mol: Chem.Mol, output_mol: Chem.Mol, 
                                          input_chiral_centers: List[Tuple[int, str]], 
                                          output_chiral_centers: List[Tuple[int, str]], 
                                          source: str) -> Dict[int, Optional[str]]:
    """
    Align molecules using MCS and map chiral centers from input to output.
    
    Args:
        input_mol: Input molecule
        output_mol: Output molecule  
        input_chiral_centers: Input chiral centers
        output_chiral_centers: Output chiral centers
        source: Source type ("PDB" or "SMILES")
        
    Returns:
        Dictionary mapping input atom indices to output configurations
    """
    from rdkit.Chem import rdFMCS
    
    # Find maximum common substructure to establish atom mapping
    mcs = rdFMCS.FindMCS([input_mol, output_mol], 
                       bondCompare=rdFMCS.BondCompare.CompareOrderExact,
                       atomCompare=rdFMCS.AtomCompare.CompareElements,
                       timeout=PROCESSING.mcs_timeout_seconds)
    
    if mcs.numAtoms <= PROCESSING.min_mcs_atoms:
        raise ValueError(f"Insufficient structural similarity for mapping (MCS: {mcs.numAtoms} atoms)")
    
    # Get the match mapping
    mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
    input_match = input_mol.GetSubstructMatch(mcs_mol)
    output_match = output_mol.GetSubstructMatch(mcs_mol)
    
    if not input_match or not output_match:
        raise ValueError("Could not establish atom mapping between structures (no substructure match)")
    
    # Create mapping from input to output atom indices
    atom_mapping = dict(zip(input_match, output_match))
    
    # Map input chiral centers to output using the mapping
    mapped_output_chiral = {}
    for input_idx, config in input_chiral_centers:
        if input_idx in atom_mapping:
            output_idx = atom_mapping[input_idx]
            # Find the corresponding output chiral center
            for out_idx, out_config in output_chiral_centers:
                if out_idx == output_idx:
                    mapped_output_chiral[input_idx] = out_config
                    break
            else:
                mapped_output_chiral[input_idx] = None  # Lost chirality
    
    logging.debug(f"Successfully mapped {len(atom_mapping)} atoms via MCS ({source} input)")
    logging.debug(f"Input chiral centers: {input_chiral_centers}")
    logging.debug(f"Mapped output chiral: {mapped_output_chiral}")
    
    return mapped_output_chiral


def read_json(json_file: str) -> dict:
    """Read a JSON file and return its contents as a dictionary."""
    if not os.path.exists(json_file):
        raise FileNotFoundError(f"JSON file not found: {json_file}")
    
    try:
        with open(json_file, 'r', encoding='utf-8') as file:
            data = json.load(file)
            return data
    except json.JSONDecodeError as e:
        raise ValueError(f"Invalid JSON format in {json_file}: {e}")
    except PermissionError:
        raise PermissionError(f"Permission denied reading {json_file}")
    except Exception as e:
        raise IOError(f"Error reading {json_file}: {e}")

def validate_json_structure(json_data: dict) -> None:
    """Warn if required keys are missing in any dataset info."""
    required_keys = ['ExpectedSummary']
    
    # Handle both dict and list formats
    if isinstance(json_data, list):
        logging.warning("JSON data is a list format. Expected dictionary with dataset names as keys.")
        for idx, info in enumerate(json_data):
            if not isinstance(info, dict):
                logging.warning(f"List item {idx} is not a dictionary")
                continue
            for key in required_keys:
                if key not in info:
                    logging.warning(f"Missing required key '{key}' in list item {idx}")
    elif isinstance(json_data, dict):
        for dataset, info in json_data.items():
            if not isinstance(info, dict):
                logging.warning(f"Dataset '{dataset}' info is not a dictionary")
                continue
            for key in required_keys:
                if key not in info:
                    logging.warning(f"Missing required key '{key}' in dataset '{dataset}'")
    else:
        logging.error(f"JSON data is neither a list nor a dictionary, but {type(json_data)}")


def build_file_path(base_dir: str, relative_path: str, filename: str) -> str:
    """Build a file path from base directory, relative path, and filename."""
    if not base_dir or not relative_path or not filename:
        return 'NA'
    
    # Normalize paths
    base_dir = os.path.normpath(base_dir)
    relative_path = os.path.normpath(relative_path)
    
    # Handle cases where relative_path is actually an absolute path
    if os.path.isabs(relative_path):
        return os.path.join(relative_path, filename)
    
    # Construct the full path
    full_path = os.path.join(base_dir, relative_path, filename)
    
    return full_path if os.path.exists(full_path) else 'NA'


def safe_get(d: Any, keys: List[Any], default: Any = "NA") -> Any:
    """Safely get a nested value from a dict or list, or return default."""
    for i, key in enumerate(keys):
        if isinstance(d, dict):
            d = d.get(key, {})
        elif isinstance(d, list) and isinstance(key, int) and len(d) > key:
            d = d[key]
        else:
            return default if i == len(keys) - 1 else {}
    return d if d != {} else default


def safe_round(value: Any, digits: int = 3) -> Any:
    """Round a value to a given number of digits, or return 'NA' if not possible."""
    try:
        return round(float(value), digits)
    except (TypeError, ValueError):
        return 'NA'


def safe_dict_get(stats: Optional[dict], key: str, default: Any = 'NA') -> Any:
    """Safely get a value from a dictionary, returning default if dict is None or key missing."""
    return stats.get(key, default) if isinstance(stats, dict) else default


def get_file_link_html(cell_val: str) -> str:
    """Generate HTML for file links with appropriate data attributes."""
    if not cell_val or cell_val == 'NA':
        return '<td></td>'
    
    file_ext = os.path.splitext(cell_val)[1].lower()
    file_name = os.path.basename(cell_val)
    return f'<td><a href="#" class="file-link" data-file="file://{cell_val}" data-ext="{file_ext}">{file_name}</a></td>'


class PipedreamPaths:
    """Helper class for constructing Pipedream directory paths."""
    
    def __init__(self, pipedream_dir: str, compound_code: str):
        """Initialize with base directory and compound code."""
        self.base_dir = Path(pipedream_dir) if pipedream_dir != 'NA' else None
        self.compound_code = compound_code if compound_code != 'NA' else None
    
    def _build_subdir_path(self, subdir_prefix: str, filename: str = '') -> Optional[Path]:
        """Build path to a file in a subdirectory."""
        if not self.base_dir or not self.compound_code:
            return None
        
        path = self.base_dir / f"{subdir_prefix}-{self.compound_code}"
        if filename:
            path = path / filename
        return path if path.exists() else None
    
    def rhofit_file(self, filename: Union[str, FileType]) -> str:
        """Get path to a file in the rhofit directory."""
        if isinstance(filename, FileType):
            filename = filename.value
        path = self._build_subdir_path(SubdirPrefix.RHOFIT.value, filename)
        return str(path) if path else 'NA'
    
    def postrefine_file(self, filename: str) -> str:
        """Get path to a file in the postrefine directory."""
        path = self._build_subdir_path(SubdirPrefix.POSTREFINE.value, filename)
        return str(path) if path else 'NA'
    
    def report_file(self, *path_parts: str) -> str:
        """Get path to a file in the report directory."""
        if not self.base_dir or not self.compound_code:
            return 'NA'
        
        report_dir = self.base_dir / f"report-{self.compound_code}"
        if not report_dir.exists():
            return 'NA'
        
        full_path = report_dir.joinpath(*path_parts)
        return str(full_path) if full_path.exists() else 'NA'
    
    def output_smiles_file(self) -> str:
        """Get path to the output SMILES file."""
        if not self.base_dir or not self.compound_code:
            return 'NA'
        
        path = self.base_dir / f"{self.compound_code}{OUTPUT_SMILES_SUFFIX}"
        return str(path)
    
    def find_electron_density_gif(self) -> str:
        """Find the electron density movie GIF file."""
        if not self.base_dir or not self.compound_code:
            return 'NA'
        
        gif_dir = self.base_dir / f"report-{self.compound_code}" / 'ligand' / 'pictures'
        if not gif_dir.exists():
            return 'NA'
        
        for gif_file in gif_dir.glob(f"*{ELECTRON_DENSITY_GIF_SUFFIX}"):
            return str(gif_file)
        
        return 'NA'


class InputPaths:
    """Helper class for finding input structure files."""
    
    def __init__(self, input_dir: str, compound_code: str):
        """Initialize with input directory and compound code."""
        self.input_dir = Path(input_dir) if input_dir != 'NA' else None
        self.compound_code = compound_code if compound_code != 'NA' else None
    
    def _find_file(self, *extensions: str) -> str:
        """Find file with any of the given extensions."""
        if not self.input_dir or not self.compound_code or not self.input_dir.exists():
            return 'NA'
        
        for ext in extensions:
            file_path = self.input_dir / f"{self.compound_code}{ext}"
            if file_path.exists():
                return str(file_path)
        
        return 'NA'
    
    def smiles_file(self) -> str:
        """Get path to SMILES file."""
        return self._find_file('.smiles', '.smi')
    
    def pdb_file(self) -> str:
        """Get path to input PDB file."""
        return self._find_file('.pdb')
    
    def cif_file(self) -> str:
        """Get path to input CIF file."""
        return self._find_file('.cif')
    
    def ligand_diagram(self) -> str:
        """Get path to ligand diagram SVG."""
        return self._find_file('.diagram.svg')
    
    def all_files(self) -> dict:
        """Get dictionary of all input file paths."""
        return {
            'smiles': self.smiles_file(),
            'pdb': self.pdb_file(),
            'cif': self.cif_file(),
            'diagram': self.ligand_diagram()
        }


# Now update the existing functions to use these helpers:

def get_ligand_png_path(input_dir: str, compound_code: str) -> str:
    """Construct the path to the ligand diagram SVG file from input directory and compound code."""
    paths = InputPaths(input_dir, compound_code)
    return paths.ligand_diagram()


def get_output_smiles_path(pipedream_dir: str, compound_code: str) -> str:
    """Construct the standard path for output SMILES files."""
    paths = PipedreamPaths(pipedream_dir, compound_code)
    return paths.output_smiles_file()


def find_input_files(input_dir: str, compound_code: str) -> dict:
    """Find input structure files (SMILES, PDB, CIF) in the input directory."""
    paths = InputPaths(input_dir, compound_code)
    result = paths.all_files()
    
    # Keep the debug logging for SMILES file
    if result['smiles'] == 'NA' and paths.input_dir and paths.input_dir.exists():
        try:
            available_files = [f.name for f in paths.input_dir.iterdir() 
                             if f.suffix in ('.smiles', '.smi', '.pdb', '.cif')]
            if available_files:
                truncated_files = available_files[:PROCESSING.max_debug_files_display]
                ellipsis = '...' if len(available_files) > PROCESSING.max_debug_files_display else ''
                logging.warning(f"SMILES file not found. Available files: {truncated_files}{ellipsis}")
        except (OSError, PermissionError):
            pass
    
    # Add debug logging for found files
    if result['pdb'] != 'NA':
        logging.debug(f"Found input PDB: {result['pdb']}")
    if result['cif'] != 'NA':
        logging.debug(f"Found input CIF: {result['cif']}")
    
    return result


class MoleculeStats(TypedDict, total=False):
    """Type definition for ligand statistics."""
    ligandid: str
    ligandcc: float
    ligandbavg: float
    ligandomin: float
    mogulzangl: float
    mogulzbond: float

class MolprobityStats(TypedDict, total=False):
    """Type definition for MolProbity statistics."""
    cbetadeviations: int
    ramaoutlierpercent: float
    ramafavoredpercent: float
    poorrotamerspercent: float
    clashscore: float
    molprobityscore: float
    rmsbonds: float
    rmsangles: float

class PostrefinementStats(TypedDict, total=False):
    """Type definition for post-refinement statistics."""
    R: float
    Rfree: float
    MeanB: float

@dataclass
class ValidationThresholds:
    """Validation thresholds for crystallographic quality assessment."""
    ligand_cc: float = 0.8
    rfree_max: float = 0.3
    r_max: float = 0.3
    mogul_z_max: float = 2.0
    rama_outlier_max: float = 0.5
    poor_rotamer_max: float = 2.0
    clash_score_max: float = 20.0
    molprobity_score_max: float = 2.0
    b_factor_ratio: float = 1.5
    rms_angles_max: float = 3.0
    rms_bonds_max: float = 0.03

THRESHOLDS = ValidationThresholds()

def convert_mtz_to_map(mtz_file_path: str, map_type: MapType = '2fofc') -> str:
    """Convert an MTZ file to a map file using gemmi."""
    if mtz_file_path == 'NA':
        return 'NA'
    
    if not os.path.exists(mtz_file_path):
        logging.warning(f"MTZ file not found: {mtz_file_path}")
        return 'NA'
        
    map_file_path = mtz_file_path.replace('.mtz', f'_{map_type}.map')
    
    if os.path.exists(map_file_path):
        logging.debug(f"Map file already exists: {map_file_path}")
        return map_file_path
    
    gemmi_available, gemmi_cmd = check_gemmi_availability()
    if not gemmi_available:
        logging.warning("gemmi command not available")
        return 'NA'
    
    try:
        if map_type == '2fofc':
            cmd = gemmi_cmd + ['sf2map', mtz_file_path, map_file_path]
        elif map_type == 'fofc':
            cmd = gemmi_cmd + ['sf2map', mtz_file_path, '-d', map_file_path]
        else:
            raise ValueError(f"Unknown map type: {map_type}")
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30, check=True)
        
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
    except (OSError, PermissionError) as e:
        logging.error(f"File system error during conversion: {e}")
        return 'NA'


def build_result(
    dataset: str,
    info: dict,
    summary_file: str = 'NA',
    ligand_stats: Optional[dict] = None,
    molprobity_stats: Optional[dict] = None,
    postrefinement_stats: Optional[dict] = None,
    electron_density_gif: str = 'NA',
    chirality_flip: str = 'Not checked',
    updated_smiles_file: str = 'NA',
    pdb_file: str = 'NA',
    mtz_file: str = 'NA',
    map_2fofc_file: str = 'NA',
    map_fofc_file: str = 'NA',
    high_resolution: str = 'NA'
) -> Dict[str, Any]:
    """Build a result dictionary for a dataset."""
    # Construct ligand SVG diagram path from Input_dir and CompoundCode
    input_dir = info.get('Input_dir', 'NA')
    compound_code = info.get('CompoundCode', 'NA')
    ligand_png = get_ligand_png_path(input_dir, compound_code)
    
    return {
        'Crystal Name': dataset,
        'Compound Code': compound_code,
        'Input Ligand Structure': ligand_png,
        'Ligand Density': electron_density_gif,
        'Chiral inversion': chirality_flip,
        'Pipedream Directory': info.get('PipedreamDirectory', 'NA'),
        'Buster Report HTML': info.get('ReportHTML', 'NA'),
        'Ligand Report HTML': info.get('LigandReportHTML', 'NA'),
        'Pipedream Summary': summary_file,
        'PDB File': pdb_file,
        'MTZ File': mtz_file,
        '2Fo-Fc Map File': map_2fofc_file,
        'Fo-Fc Map File': map_fofc_file,
        'Output SMILES': updated_smiles_file,
        'Ligand ID': safe_dict_get(ligand_stats, 'ligandid'),
        'Ligand CC': safe_dict_get(ligand_stats, 'ligandcc'),
        'R': safe_round(safe_dict_get(postrefinement_stats, 'R')),
        'Rfree': safe_round(safe_dict_get(postrefinement_stats, 'Rfree')),
        'Ligand avg B factor': safe_dict_get(ligand_stats, 'ligandbavg'),
        'Ligand occupancy': safe_dict_get(ligand_stats, 'ligandomin'),
        'Mogul Z angle': safe_dict_get(ligand_stats, 'mogulzangl'),
        'Mogul Z bond': safe_dict_get(ligand_stats, 'mogulzbond'),
        'c beta deviations': safe_dict_get(molprobity_stats, 'cbetadeviations'),
        'Rama outlier percent': safe_dict_get(molprobity_stats, 'ramaoutlierpercent'),
        'Rama favored percent': safe_dict_get(molprobity_stats, 'ramafavoredpercent'),
        'Poor rotamer percent': safe_dict_get(molprobity_stats, 'poorrotamerspercent'),
        'Clash score': safe_dict_get(molprobity_stats, 'clashscore'),
        'Mol probity score': safe_dict_get(molprobity_stats, 'molprobityscore'),
        'RMS bonds': safe_dict_get(molprobity_stats, 'rmsbonds'),
        'RMS angles': safe_dict_get(molprobity_stats, 'rmsangles'),
        'Mean B factor': safe_dict_get(postrefinement_stats, 'MeanB'),
        'High resolution (A)': high_resolution,
        'Comments': '',
        'Export to XCE': 'True' if (isinstance(ligand_stats, dict) and ligand_stats.get('ligandcc', 0) > THRESHOLDS.ligand_cc) else 'False'
    }


def collect_results_from_json(json_data: dict) -> List[Dict[str, Any]]:
    """Collect results from the JSON data for all datasets."""
    # Input validation
    if not json_data:
        logging.warning("No data provided to collect_results_from_json")
        return []
    
    results = []
    total_datasets = len(json_data)
    logging.info(f"Processing {total_datasets} datasets")
    print(f"Processing {total_datasets} datasets...")  # Console progress
    
    # Check gemmi availability once for all datasets (not per dataset)
    gemmi_available, gemmi_cmd = check_gemmi_availability()
    if gemmi_available:
        logging.info("gemmi is available for map generation")
    else:
        logging.info("gemmi not available, map generation will be skipped")
    
    for i, (dataset, info) in enumerate(json_data.items(), 1):
        logging.info(f"Processing dataset {i}/{total_datasets}: {dataset}")
        
        if i % 10 == 0 or i == 1 or i == total_datasets:
            print(f"Progress: {i}/{total_datasets} datasets processed")
        
        summary_file = info.get('ExpectedSummary', 'NA')
        
        if not os.path.exists(summary_file):
            logging.warning(f"Summary file not found for {dataset}: {summary_file}")
            results.append(build_result(dataset, info, summary_file, updated_smiles_file="NA"))
            continue
        
        try:
            with open(summary_file, 'r', encoding='utf-8') as f:
                summary = json.load(f)
            logging.debug(f"Successfully loaded summary for {dataset}")
        except json.JSONDecodeError as e:
            logging.error(f"Error decoding JSON in {summary_file} for dataset {dataset}: {e}")
            results.append(build_result(dataset, info, summary_file, updated_smiles_file="NA"))
            continue
        
        # Extract high resolution from dataprocessing.inputdata.reshigh
        high_resolution = safe_get(summary, ['dataprocessing', 'inputdata', 'reshigh'], 'NA')
        if high_resolution != 'NA':
            try:
                high_resolution = f"{float(high_resolution):.2f}"
            except (ValueError, TypeError):
                high_resolution = 'NA'
        
        # Extract validation statistics
        ligand = safe_get(summary, ['ligandfitting', 'ligands', 0])
        ligand_stats = safe_get(ligand, ['validationstatistics', 'ligandstatistics', 0])
        molprobity_stats = safe_get(ligand, ['validationstatistics', 'molprobity'])
        postrefinement_stats = safe_get(ligand, ['postrefinement', 1])
        
        # Extract key directory and compound information
        pipedream_dir = info.get('PipedreamDirectory', 'NA')
        compound_code = info.get('CompoundCode', 'NA')
        input_dir = info.get('Input_dir', 'NA')
        
        # Use compound_code or fallback to dataset name
        if compound_code == 'NA':
            compound_code = dataset
        
        logging.debug(f"Dataset: {dataset}, Compound: {compound_code}")
        
        # Initialize path helpers
        pipedream_paths = PipedreamPaths(pipedream_dir, compound_code)
        input_paths = InputPaths(input_dir, compound_code)
        
        # Get all file paths using helper classes
        pdb_file = pipedream_paths.rhofit_file(REFINE_PDB_FILENAME)
        output_cif_file = pipedream_paths.rhofit_file(BEST_CIF_FILENAME)
        mtz_file = pipedream_paths.postrefine_file(REFINE_MTZ_FILENAME)
        electron_density_gif = pipedream_paths.find_electron_density_gif()
        
        # Get input files
        input_files = input_paths.all_files()
        input_smi_path = input_files['smiles']
        input_pdb_path = input_files['pdb']
        input_cif_path = input_files['cif']
        
        # Fallback to JSON parsing if needed
        if pdb_file == 'NA' or output_cif_file == 'NA':
            logging.debug(f"Using fallback JSON parsing for missing files")
            postrefinement = safe_get(summary, ['pipedream_outputs', 'ligandfitting', 'ligands', 0, 'postrefinement'], [])
            if isinstance(postrefinement, list):
                for entry in postrefinement:
                    if entry.get('description') == 'final':
                        if entry.get('type') == 'model' and entry.get('format') == 'PDB' and pdb_file == 'NA':
                            pdb_file = build_file_path(
                                pipedream_dir,
                                entry.get('relative_path', ''),
                                entry.get('filename', '')
                            )
                        elif entry.get('type') == 'model' and entry.get('format') == 'CIF' and output_cif_file == 'NA':
                            output_cif_file = build_file_path(
                                pipedream_dir,
                                entry.get('relative_path', ''),
                                entry.get('filename', '')
                            )
        
        # Generate map files
        if gemmi_available and mtz_file != 'NA':
            logging.debug(f"Generating map files for: {mtz_file}")
            map_2fofc_file = convert_mtz_to_map(mtz_file, '2fofc')
            map_fofc_file = convert_mtz_to_map(mtz_file, 'fofc')
        else:
            map_2fofc_file = 'NA'
            map_fofc_file = 'NA'
        
        # Use best.pdb for chirality analysis
        chirality_pdb_file = pipedream_paths.rhofit_file(BEST_PDB_FILENAME)
        
        # Perform chirality analysis
        # Always prioritize using CIF input file when available
        # Input: Must have CIF (PDB optional), Output: Must have both CIF and PDB
        if (input_cif_path != 'NA' and output_cif_file != 'NA' and 
            chirality_pdb_file != 'NA' and input_smi_path != 'NA'):
            
            logging.debug(f"Using CIF+PDB approach for chirality analysis of {dataset}")
            chirality_flip, updated_smiles = compare_chiral_centers_cif_pdb(
                input_cif_path, input_pdb_path, output_cif_file, chirality_pdb_file, input_smi_path, dataset)
            
            # Save SMILES file
            smiles_file = _save_smiles_file(pipedream_dir, compound_code, 
                                            chirality_flip, updated_smiles, "CIF+PDB", dataset)
            
        # Fallback to original SMILES+PDB approach if CIF requirements not met
        elif input_smi_path != 'NA' and chirality_pdb_file != 'NA':
            missing_cif_info = []
            if input_cif_path == 'NA':
                missing_cif_info.append("input CIF")
            if output_cif_file == 'NA':
                missing_cif_info.append("output CIF")
            
            logging.debug(f"CIF requirements not met (missing: {', '.join(missing_cif_info)}), using fallback SMILES+PDB approach for {dataset}")
            chirality_flip, updated_smiles = detect_chiral_inversion(input_smi_path, chirality_pdb_file, input_pdb_path, pipedream_dir, compound_code)
            
            # Construct expected path to analysis file - directly without condition
            expected_smiles_file = get_output_smiles_path(pipedream_dir, compound_code)
            smiles_file = expected_smiles_file if os.path.exists(expected_smiles_file) else 'NA'
        
        else:
            missing_info = []
            if input_smi_path == 'NA':
                missing_info.append("SMILES file")
            if chirality_pdb_file == 'NA':
                missing_info.append("output PDB file")
            chirality_flip = f"Missing: {', '.join(missing_info)}"
            smiles_file = "NA"
            logging.warning(f"Skipping chirality analysis for {dataset}: {chirality_flip}")
        
        # Build result
        result = build_result(
            dataset=dataset,
            info=info,
            summary_file=summary_file,
            ligand_stats=ligand_stats,
            molprobity_stats=molprobity_stats,
            postrefinement_stats=postrefinement_stats,
            electron_density_gif=electron_density_gif,
            chirality_flip=chirality_flip,
            updated_smiles_file=smiles_file,
            pdb_file=pdb_file,
            mtz_file=mtz_file,
            map_2fofc_file=map_2fofc_file,
            map_fofc_file=map_fofc_file,
            high_resolution=high_resolution
        )
        results.append(result)
    
    return results


def get_cell_style_and_value(col: str, cell_val: Any, row: pd.Series) -> tuple[str, str]:
    """Get the appropriate styling and display value for a table cell."""
    style = ''
    display_val = str(cell_val)
    
    # Try to format numbers to 2 decimal places if possible
    try:
        num_val = float(cell_val)
        display_val = f"{num_val:.2f}"
    except (ValueError, TypeError):
        pass
    
    # Apply numeric highlighting rules
    try:
        val = float(cell_val)
        if col == 'Rfree' and val > THRESHOLDS.rfree_max:
            style = 'background-color: #ffcccc; font-weight: bold;'
        elif col == 'R' and val > THRESHOLDS.r_max:
            style = 'background-color: #ffcccc; font-weight: bold;'
        elif col == 'Mogul Z angle' and val > THRESHOLDS.mogul_z_max:
            style = 'background-color: #ffcccc; font-weight: bold;'
        elif col == 'Mogul Z bond' and val > THRESHOLDS.mogul_z_max:
            style = 'background-color: #ffcccc; font-weight: bold;'
        elif col == 'Ligand avg B factor':
            try:
                mean_b = float(row['Mean B factor'])
                if mean_b != 0 and val / mean_b >= THRESHOLDS.b_factor_ratio:
                    style = 'background-color: #ffcccc; font-weight: bold;'
            except (ValueError, TypeError, ZeroDivisionError):
                pass
        elif col == 'Rama outlier percent' and val > THRESHOLDS.rama_outlier_max:
            style = 'background-color: #ffcccc; font-weight: bold;'
        elif col == 'Poor rotamer percent' and val > THRESHOLDS.poor_rotamer_max:
            style = 'background-color: #ffcccc; font-weight: bold;'
        elif col == 'Clash score' and val > THRESHOLDS.clash_score_max:
            style = 'background-color: #ffcccc; font-weight: bold;'
        elif col == 'Mol probity score' and val > THRESHOLDS.molprobity_score_max:
            style = 'background-color: #ffcccc; font-weight: bold;'
        elif col == 'RMS angles' and val > THRESHOLDS.rms_angles_max:
            style = 'background-color: #ffcccc; font-weight: bold;'
        elif col == 'RMS bonds' and val > THRESHOLDS.rms_bonds_max:
            style = 'background-color: #ffcccc; font-weight: bold;'
    except (ValueError, TypeError):
        pass
    
    return style, display_val


def get_cell_html(col: str, cell_val: Any, row: pd.Series) -> str:
    """Generate HTML for a specific table cell based on column type."""
    if col == 'Export to XCE':
        return f'<td><input type="checkbox" {"checked" if cell_val == "True" else ""}></td>'
    elif col == 'Comments':
        return f'<td><input type="text" value="{cell_val}"></td>'
    elif col == 'Input Ligand Structure':
        return f'<td><img src="file://{cell_val}" alt="Ligand Image" width="150"></td>' if cell_val and cell_val != 'NA' else '<td></td>'
    elif col == 'Ligand Density':
        return f'<td><img src="file://{cell_val}" alt="Ligand Density" width="180" style="max-width:180px;max-height:120px;"></td>' if cell_val and cell_val != 'NA' else '<td></td>'
    elif col == 'Chiral inversion':
        return f'<td style="color:red;font-weight:bold">{cell_val}</td>' if cell_val and cell_val not in ['Output stereochemistry matches input', 'Not checked'] else f'<td>{cell_val}</td>'
    elif col in [
        'Pipedream Directory', 'Buster Report HTML', 'Ligand Report HTML', 'Pipedream Summary',
        'PDB File', 'MTZ File', '2Fo-Fc Map File', 'Fo-Fc Map File', 'Output SMILES']:
        return get_file_link_html(cell_val)
    else:
        # Regular cell with potential highlighting
        style, display_val = get_cell_style_and_value(col, cell_val, row)
        return f'<td style="{style}">{display_val}</td>'


def save_results_to_html(results: List[Dict[str, Any]], output_file: str, open_browser: bool = True) -> None:
    """Save results to an interactive HTML file."""
    logging.info(f"Generating HTML output with {len(results)} results")
    print(f"Generating HTML report...")
    
    df = pd.DataFrame(results)

    # Ensure 'Ligand CC' is numeric for sorting, non-numeric values become NaN
    df['Ligand CC'] = pd.to_numeric(df['Ligand CC'], errors='coerce')
    # Auto sort by Ligand CC descending (NaNs will be at the end)
    df = df.sort_values(by='Ligand CC', ascending=False)
    logging.debug(f"Sorted results by Ligand CC")

    # Define the desired column order (ensure it's always defined)
    column_order = COLUMN_ORDER

    # Build table rows with actual cell values in the new order
    table_rows = []
    for _, row in df.iterrows():
        row_html = '<tr>'
        for col in column_order:
            row_html += get_cell_html(col, row[col], row)
        row_html += '</tr>'
        table_rows.append(row_html)

    # Render table header HTML outside JS template
    table_header_html = ''.join([f'<th>{col}</th>' for col in column_order])
    
    # Calculate Ligand CC column index for sorting
    ligand_cc_index = column_order.index('Ligand CC') if 'Ligand CC' in column_order else 0

    html_content = f"""
    <!DOCTYPE html>
    <html lang=\"en\">
    <head>
        <meta charset=\"UTF-8\">
        <title>Pipedream XChem Results</title>
        <link rel=\"stylesheet\" href=\"https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css\">
        <link rel=\"stylesheet\" href=\"https://cdn.datatables.net/1.13.4/css/dataTables.bootstrap5.min.css\">
        <script src=\"https://code.jquery.com/jquery-3.6.0.min.js\"></script>
        <script src=\"https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js\"></script>
        <script src=\"https://cdn.datatables.net/1.13.4/js/jquery.dataTables.min.js\"></script>
        <script src=\"https://cdn.datatables.net/1.13.4/js/dataTables.bootstrap5.min.js\"></script>
        <style>
            html, body {{
                width: 100vw;
                min-width: 100vw;
                margin: 0;
                padding: 0;
                background-color: white;
            }}
            .container {{
                width: 100vw !important;
                max-width: 100vw !important;
                padding: 2rem 0.5rem;
            }}
            .table-responsive {{
                width: 100vw !important;
                max-width: 100vw !important;
                overflow-x: auto;
            }}
            table {{
                width: 100% !important;
            }}
            table img {{
                max-width: 150px;
                height: auto;
            }}
            /* Move DataTables search box to the left */
            div.dataTables_filter {{
                float: left !important;
                text-align: left !important;
            }}
            /* Remove top scroll bar by disabling scroll on head */
            .dataTables_scrollHead {{
                overflow-x: hidden !important;
                overflow-y: hidden;
            }}
            .dataTables_scrollFoot {{
                overflow-x: auto !important;
                overflow-y: hidden;
            }}
            .dataTables_scrollBody {{
                overflow-x: auto !important;
                overflow-y: auto !important;
                max-height: 60vh;
            }}
            /* Highlighted cells style is inline */
            /* Make checkboxes in Export to XCE column bigger */
            td input[type=checkbox] {{
                transform: scale(1.5);
                width: 20px;
                height: 20px;
                margin: 4px;
            }}
            /* File link styling */
            .file-link {{
                color: #0066cc;
                text-decoration: none;
                cursor: pointer;
            }}
            .file-link:hover {{
                color: #0052a3;
                text-decoration: underline;
            }}
            /* Modal styling */
            .modal-body iframe {{
                border: none;
                width: 100%;
                height: 70vh;
            }}
            .modal-xl {{
                max-width: 90vw;
            }}
        </style>
        <script>
            // Embed the original JSON filename for export
            const originalJsonFilename = "{os.path.basename(output_file).replace('.html', '.json')}";
            
            // Function to attach modal event handlers
            function attachModalHandlers() {{
                // Modal popup for file links
                $(document).off('click', '.file-link').on('click', '.file-link', function(e) {{
                    e.preventDefault();
                    var fileUrl = $(this).data('file');
                    var ext = $(this).data('ext');
                    var fileName = $(this).text();
                    var modalBody = $('#fileModal .modal-body');
                    var modalTitle = $('#fileModalLabel');
                    
                    modalBody.empty();
                    modalTitle.text('File Viewer - ' + fileName);
                    
                    if(ext === '.html') {{
                        modalBody.append(`<iframe src="${{fileUrl}}" style="width:100%;height:70vh;border:none;"></iframe>`);
                    }} else if(ext === '.pdb' || ext === '.mtz' || ext === '.map') {{
                        modalBody.append(`
                            <div class="alert alert-info">
                                <h6>File Information:</h6>
                                <p><strong>File:</strong> ${{fileName}}</p>
                                <p><strong>Path:</strong> ${{fileUrl.replace('file://', '')}}</p>
                                <p><strong>Type:</strong> ${{ext.toUpperCase()}} file</p>
                            </div>
                            <div class="d-grid gap-2">
                                <a href="${{fileUrl}}" class="btn btn-primary" download>Download File</a>
                                <button class="btn btn-secondary" onclick="navigator.clipboard.writeText('${{fileUrl.replace('file://', '')}}'); alert('Path copied to clipboard!')">Copy Path</button>
                            </div>
                        `);
                    }} else if(ext === '.png' || ext === '.jpg' || ext === '.jpeg' || ext === '.gif') {{
                        modalBody.append(`<img src="${{fileUrl}}" class="img-fluid" alt="${{fileName}}" style="max-width:100%;height:auto;">`);
                    }} else {{
                        modalBody.append(`
                            <div class="alert alert-info">
                                <h6>File Information:</h6>
                                <p><strong>File:</strong> ${{fileName}}</p>
                                <p><strong>Path:</strong> ${{fileUrl.replace('file://', '')}}</p>
                            </div>
                            <div class="d-grid gap-2">
                                <a href="${{fileUrl}}" class="btn btn-primary" target="_blank">Open File</a>
                                <button class="btn btn-secondary" onclick="navigator.clipboard.writeText('${{fileUrl.replace('file://', '')}}'); alert('Path copied to clipboard!')">Copy Path</button>
                            </div>
                        `);
                    }}
                    
                    var modal = new bootstrap.Modal(document.getElementById('fileModal'));
                    modal.show();
                }});
            }}
            
            
            $(document).ready(function() {{
                $('#resultsTable').DataTable({{
                    responsive: true,

                    pageLength: 25,
                    order: [[{ligand_cc_index}, 'desc']], // Ligand CC column index (0-based)
                    columnDefs: [
                        {{ targets: '_all', type: 'num' }}
                    ],
                    scrollX: true,
                    scrollY: '60vh',
                    scrollCollapse: true,
                    fixedHeader: true,
                    dom: 'Bflrt<"bottom-scrollbar"ip>',
                }});
                
                // Attach modal handlers initially
                attachModalHandlers();
                
                // Only sync scroll for bottom
                function syncScroll() {{
                    var scrollBody = $('.dataTables_scrollBody')[0];
                    var scrollFoot = $('.dataTables_scrollFoot')[0];
                    if(scrollBody && scrollFoot) {{
                        scrollFoot.scrollLeft = scrollBody.scrollLeft;
                    }}
                }}
                $('.dataTables_scrollBody').on('scroll', syncScroll);
            }});
            // --- Add exportToJSON function ---
            function exportToJSON() {{
                var table = $('#resultsTable').DataTable();
                var data = [];
                var headers = [];
                // Get headers
                $('#resultsTable thead th').each(function() {{
                    headers.push($(this).text());
                }});
                // Get data from visible rows
                table.rows({{ search: 'applied' }}).every(function(rowIdx, tableLoop, rowLoop) {{
                }});
                // Get data from visible rows
                table.rows({{ search: 'applied' }}).every(function(rowIdx, tableLoop, rowLoop) {{
                    var rowNode = this.node();
                    var rowData = {{}};
                    $(rowNode).find('td').each(function(i, cell) {{
                        var col = headers[i];
                        if (col === 'Export to XCE') {{
                            var checked = $(cell).find('input[type=checkbox]').is(':checked');
                            rowData[col] = checked ? 'True' : 'False';
                        }} else if (col === 'Comments') {{
                            var val = $(cell).find('input[type=text]').val();
                            rowData[col] = val !== undefined ? val : '';
                        }} else if (col === 'Ligand Structure' || col === 'Ligand Density') {{
                            var img = $(cell).find('img');
                            if (img.length > 0) {{
                                var src = img.attr('src') || '';
                                // Remove file:// prefix if present
                                rowData[col] = src.startsWith('file://') ? src.substring(7) : src;
                            }} else {{
                                rowData[col] = '';
                            }}
                        }} else if ([
                            'Pipedream Directory', 'Buster Report HTML', 'Ligand Report HTML', 'Pipedream Summary',
                            'PDB File', 'MTZ File', '2Fo-Fc Map File', 'Fo-Fc Map File', 'Output SMILES File'
                                               ].includes(col)) {{
                            var link = $(cell).find('a.file-link');
                            if (link.length > 0) {{
                                var file = link.data('file') || '';
                                // Remove file:// prefix if present
                                rowData[col] = file.startsWith('file://') ? file.substring(7) : file;
                            }} else {{
                                rowData[col] = '';
                            }}
                        }} else {{
                            // Remove HTML tags from cell
                            var cellText = $(cell).text();
                            rowData[col] = cellText;
                        }}
                    }});
                    data.push(rowData);
                }});
                var jsonStr = JSON.stringify(data, null, 2);
                var blob = new Blob([jsonStr], {{type: 'application/json'}});
                var a = document.createElement('a');
                a.href = URL.createObjectURL(blob);
                a.download = originalJsonFilename;
                document.body.appendChild(a);
                a.click();
                document.body.removeChild(a);
            }}
            // --- Add loadJsonAndUpdateTable function ---
            function loadJsonAndUpdateTable(event) {{
                var file = event.target.files[0];
                if (!file) return;
                var reader = new FileReader();
                reader.onload = function(e) {{
                    try {{
                        var jsonData = JSON.parse(e.target.result);
                        // Remove current table
                        $('#resultsTable').DataTable().destroy();
                        $('#resultsTable tbody').empty();
                        // Rebuild table body
                        jsonData.forEach(function(row) {{
                            var rowHtml = '<tr>';
                            Object.values(row).forEach(function(cell, i) {{
                                var col = Object.keys(row)[i];
                                if (col === 'Export to XCE') {{
                                    rowHtml += '<td><input type="checkbox" ' + (cell === 'True' ? 'checked' : '') + '></td>';
                                }} else if (col === 'Comments') {{
                                    rowHtml += '<td><input type="text" value="' + (cell || '') + '"></td>';
                                }} else if (col === 'Ligand Structure') {{
                                    rowHtml += (cell && cell !== 'NA') ? '<td><img src="file://' + cell + '" alt="Ligand Image" width="150"></td>' : '<td></td>';
                                }} else if (col === 'Ligand Density') {{
                                    rowHtml += (cell && cell !== 'NA') ? '<td><img src="file://' + cell + '" alt="Ligand Density" width="180" style="max-width:180px;max-height:120px;"></td>' : '<td></td>';                        }} else if (col === 'Chiral inversion') {{
                            rowHtml += (cell && cell !== 'Output stereochemistry matches input' && cell !== 'Not checked') ? '<td style="color:red;font-weight:bold">' + cell + '</td>' : '<td>' + cell + '</td>';
                        }} else if (col === 'Pipedream Directory' || col === 'Buster Report HTML' || col === 'Ligand Report HTML' || col === 'Pipedream Summary' || col === 'PDB File' || col === 'MTZ File' || col === '2Fo-Fc Map File' || col === 'Fo-Fc Map File' || col === 'Output SMILES') {{
                                    if (cell && cell !== 'NA') {{
                                        var fileExt = '.' + cell.split('.').pop().toLowerCase();
                                        var fileName = cell.split(/[\\/]/).pop();
                                        rowHtml += '<td><a href="#" class="file-link" data-file="file://' + cell + '" data-ext="' + fileExt + '">' + fileName + '</a></td>';
                                    }} else {{
                                        rowHtml += '<td></td>';
                                    }}
                                }} else {{
                                    rowHtml += '<td>' + (cell || '') + '</td>';
                                }}
                            }});
                            rowHtml += '</tr>';
                            $('#resultsTable tbody').append(rowHtml);
                        }});
                        // Re-initialize DataTable
                        $('#resultsTable').DataTable({{
                            responsive: true,
                            pageLength: 25,
                            order: [[{ligand_cc_index}, 'desc']],
                            columnDefs: [{{ targets: '_all', type: 'num' }}],
                            scrollX: true,
                            scrollY: '60vh',
                            scrollCollapse: true,
                            fixedHeader: true,
                            dom: 'Bflrt<"bottom-scrollbar"ip>',
                        }});
                        
                        // Reattach modal handlers after table rebuild
                        attachModalHandlers();
                        
                        alert('Table updated successfully with ' + jsonData.length + ' records.');
                    }} catch (err) {{
                        alert('Error parsing JSON file: ' + err.message);
                        console.error('JSON parsing error:', err);
                    }}
                }};
                reader.readAsText(file);
            }}
            
            function selectAllExport() {{
                $('#resultsTable tbody input[type=checkbox]').prop('checked', true);
            }}
            function unselectAllExport() {{
                $('#resultsTable tbody input[type=checkbox]').prop('checked', false);
            }}
        </script>
    </head>
    <body>
        <div class="container">
            <div class="mb-3 d-flex flex-wrap gap-2 align-items-center">
                <input type="file" id="jsonFileInput" accept="application/json" style="display:none" onchange="loadJsonAndUpdateTable(event)">
                <button class="btn btn-warning" onclick="document.getElementById('jsonFileInput').click()">Load JSON and Update Table</button>
                <button class="btn btn-success" onclick="exportToJSON()">Export Table to JSON</button>
                <button class="btn btn-primary" onclick="selectAllExport()">Select All</button>
                <button class="btn btn-secondary" onclick="unselectAllExport()">Unselect All</button>
            </div>
            <h1 class="text-center">Pipedream XChem Results</h1>
            <div class="table-responsive">
                <table id="resultsTable" class="table table-striped table-bordered">
                    <thead>
                        <tr>
                            {table_header_html}
                        </tr>
                    </thead>
                    <tbody>
                        {''.join(table_rows)}
                    </tbody>
                </table>
            </div>
        </div>
        <!-- Modal HTML -->
        <div class="modal fade" id="fileModal" tabindex="-1" aria-labelledby="fileModalLabel" aria-hidden="true">
          <div class="modal-dialog modal-xl">
            <div class="modal-content">
              <div class="modal-header">
                <h5 class="modal-title" id="fileModalLabel">File Viewer</h5>
                <button type="button" class="btn-close" data-bs-dismiss="modal" aria-label="Close"></button>
              </div>
              <div class="modal-body" id="fileModalBody">
                <!-- Content will be dynamically loaded here -->
              </div>
              <div class="modal-footer">
                <button type="button" class="btn btn-secondary" data-bs-dismiss="modal">Close</button>
              </div>
            </div>
          </div>
        </div>
    </body>
    </html>
    """

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_content)

    logging.info(f"HTML report saved to {os.path.basename(output_file)}")
    print(f"HTML report saved: {os.path.basename(output_file)}")
    
    if open_browser:
        webbrowser.open(f'file://{os.path.abspath(output_file)}')


def save_results_to_json(results: List[Dict[str, Any]], output_file: str) -> None:
    """Save results to a JSON file with a consistent column order."""
    # Use the same column order as for HTML
    column_order = COLUMN_ORDER
    ordered_results = [OrderedDict((col, row.get(col, '')) for col in column_order) for row in results]
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(ordered_results, f, indent=2)
    logging.info(f"Results saved to {os.path.basename(output_file)}")



def main() -> None:
    """Main entry point for collating Pipedream results."""
    parser = argparse.ArgumentParser(
        description='Collate Pipedream analysis results into interactive reports.',
        epilog='Example: python collate_pipedream_results.py --input results.json --output-dir /path/to/output --format both'
    )
    
    # Version argument
    parser.add_argument(
        '--version', action='version', version=f'%(prog)s {VERSION}',
        help="Show version number"
    )
    parser.add_argument(
        '--input', '--json', '-i',
        dest='input_json',
        required=True, 
        help='Path to the JSON file with dataset metadata'
    )
    parser.add_argument(
        '--output-dir',
        help='Output directory for generated reports (default: same as input JSON)'
    )
    parser.add_argument(
        '--output-name',
        help='Base name for output files without extension (default: timestamped Pipedream_results_YYYYMMDD_HHMMSS)'
    )
    parser.add_argument(
        '--format',
        choices=['json', 'html', 'both'],
        default='both',
        help='Output format(s) to generate (default: both)'
    )
    parser.add_argument(
        '--no-browser',
        action='store_true',
        help='Don\'t automatically open HTML report in browser'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose console output'
    )
    parser.add_argument(
        '--log-level',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        default='DEBUG',
        help='Set logging level (default: DEBUG)'
    )
    
    args = parser.parse_args()
    
    # Validate input file exists
    if not os.path.isfile(args.input_json):
        parser.error(f"Input JSON file not found: {args.input_json}")
    
    # Validate input file extension
    if not args.input_json.lower().endswith('.json'):
        parser.error("Input file must have a .json extension")
    
    # Determine output directory
    if args.output_dir:
        output_dir = args.output_dir
        try:
            os.makedirs(output_dir, exist_ok=True)
        except PermissionError:
            parser.error(f"Cannot create output directory: {output_dir} (permission denied)")
        except Exception as e:
            parser.error(f"Cannot create output directory: {output_dir} ({e})")
    else:
        output_dir = os.path.dirname(args.input_json)
    
    # Determine output file names
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    if args.output_name:
        base_name = args.output_name
    else:
        base_name = f'Pipedream_results_{timestamp}'
    
    output_file_json = os.path.join(output_dir, f'{base_name}.json')
    output_file_html = os.path.join(output_dir, f'{base_name}.html')
    
    log_file = os.path.join(output_dir, f'{base_name}.log')
    
    # Configure logging with different levels for file vs console
    logger = logging.getLogger()
    logger.setLevel(getattr(logging, args.log_level))
    
    # File handler - captures everything based on log level
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(getattr(logging, args.log_level))
    file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(file_formatter)
    
    # Console handler - verbose mode or essential messages only
    console_handler = logging.StreamHandler()
    if args.verbose:
        console_handler.setLevel(logging.INFO)
    else:
        console_handler.setLevel(logging.WARNING)
    console_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console_handler.setFormatter(console_formatter)
    
    # Add a custom filter for console when not in verbose mode
    if not args.verbose:
        class ConsoleFilter(logging.Filter):
            def filter(self, record):
                # Allow WARNING and above
                if record.levelno >= logging.WARNING:
                    return True
                # Allow specific INFO messages for terminal
                if record.levelno == logging.INFO:
                    essential_messages = [
                        "JSON file read successfully",
                        "Processing dataset:",
                        "Results saved to",
                        "Result collection completed"
                    ]
                    # Check for exact matches or "Processing X datasets" pattern
                    message = record.getMessage()
                    if any(msg in message for msg in essential_messages):
                        return True
                    # Also allow "Processing X datasets" pattern
                    if message.startswith("Processing ") and message.endswith(" datasets"):
                        return True
                    return False
        
        console_handler.addFilter(ConsoleFilter())
    
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    
    logging.info("Starting Pipedream results collation")
    print("Starting Pipedream results collation...")
    
    logging.info(f"Reading JSON file: {args.input_json}")
    json_data = read_json(args.input_json)
    
    logging.info("Validating JSON structure")
    validate_json_structure(json_data)
    
    logging.info("Collecting results from datasets")
    results = collect_results_from_json(json_data)
    
    # Generate outputs based on format selection
    if args.format in ['json', 'both']:
        logging.info(f"Saving results to JSON: {output_file_json}")
        print(f"Saving JSON results...")
        save_results_to_json(results, output_file_json)
    
    if args.format in ['html', 'both']:
        logging.info(f"Generating HTML report: {output_file_html}")
        # Pass the no_browser flag to the HTML generation function
        save_results_to_html(results, output_file_html, open_browser=not args.no_browser)
    
    logging.info("Result collection completed.")
    print("Pipedream results collation completed!")
    
    # Print summary of what was generated
    print(f"\nOutput summary:")
    if args.format in ['json', 'both']:
        print(f"  JSON: {output_file_json}")
    if args.format in ['html', 'both']:
        print(f"  HTML: {output_file_html}")
    print(f"  Log:  {log_file}")


if __name__ == "__main__":
    main()
