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

Author: DFearon
Date: November 2025
"""

import os
import json
import pandas as pd
import logging
import argparse
import traceback
from datetime import datetime
import webbrowser
from collections import OrderedDict
from typing import Any, Dict, List, Tuple, Optional, TypedDict, Union
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D
import gemmi
from dataclasses import dataclass
from pathlib import Path
from enum import Enum
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Import extracted modules
from pipedream_plots import (
    generate_pca_plots,
    generate_ligand_quality_plot,
    generate_pdb_percentile_quality_plot,
    generate_individual_pdb_percentile_plot,
    generate_spider_plot,
    generate_individual_batch_percentile_plot
)
from pipedream_thresholds import ValidationThresholds, THRESHOLDS

# Try to import plotly for interactive plots
try:
    import plotly.graph_objects as go
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False
    logging.debug("plotly not available, will use static matplotlib plots")

# Try to import tqdm for progress bar
try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False
    logging.debug("tqdm not available, progress bar will not be shown")

# Import Jinja2 for HTML templating
from jinja2 import Environment, FileSystemLoader, select_autoescape

# Version information
VERSION = "1.0.2"

# Constants
@dataclass
class ProcessingConfig:
    """Configuration for molecular processing."""
    min_ligand_atoms: int = 5
    mcs_timeout_seconds: int = 30
    min_mcs_atoms: int = 5
    max_debug_files_display: int = 5

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

class Config:
    """Configuration class for external tool paths and settings."""
    
    def __init__(self):
        """Initialize configuration with default paths."""
        # Default paths to search for PDB PCA reference data
        self.pca_reference_paths = [
            Path(__file__).parent / 'Data' / 'PCA_PDB_references.csv',  # Next to script
            Path(__file__).parent / 'PCA_PDB_references.csv',  # In script directory
            Path('Data/PCA_PDB_references.csv'),  # Relative path
            Path('PCA_PDB_references.csv')  # Current directory
        ]

class PCAReference:
    """
    Class to load PDB archive reference data and calculate PCA percentile scores.
    Based on calculateCompositeScore.py by Chenghua Shao (2020).
    
    This provides absolute quality assessment by comparing against the entire PDB archive,
    in contrast to batch-relative scoring.
    """
    
    def __init__(self):
        """Initialize empty reference data structures."""
        self.d_ref = {}  # Reference data from PDB archive
        self.d_par = {}  # Calculated PCA parameters (means, stds, loadings)
        self.loaded = False
    
    def load_reference(self, filepath: str) -> bool:
        """
        Load reference PDB archive data from CSV file.
        
        Args:
            filepath: Path to PCA_PDB_references.csv file
            
        Returns:
            True if loaded successfully, False otherwise
        """
        try:
            self._read_pca_references(filepath)
            self._calculate_pca_parameters()
            self.loaded = True
            logging.info(f"Loaded PDB PCA reference data from: {filepath}")
            return True
        except Exception as e:
            logging.warning(f"Failed to load PDB PCA reference data from {filepath}: {e}")
            self.loaded = False
            return False
    
    def _read_pca_references(self, filepath: str):
        """Load reference data as dictionary with columns as keys."""
        import csv
        
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Reference file not found: {filepath}")
        
        with open(filepath, 'r') as f:
            reader = csv.reader(f)
            header = next(reader, None)
            
            if header is None:
                raise ValueError("Empty CSV file")
            
            # Initialize lists for each column (skip first column - usually ID)
            for item in header[1:]:
                self.d_ref[item] = []
            
            # Read data rows
            for line in reader:
                for item, value in zip(header[1:], line[1:]):
                    self.d_ref[item].append(float(value))
        
        logging.debug(f"Loaded {len(self.d_ref.get('rsr', []))} reference structures from PDB archive")
    
    def _calculate_pca_parameters(self):
        """Calculate PCA parameters (mean, std, loading) for each variable."""
        from math import sqrt
        import numpy as np
        
        for var in ["rsr", "rscc", "mogul_bonds_rmsz", "mogul_angles_rmsz"]:
            if var not in self.d_ref:
                raise ValueError(f"Missing required column '{var}' in reference data")
            
            self.d_par[var] = {}
            l_value = self.d_ref[var]
            self.d_par[var]["mean"] = np.mean(l_value)
            self.d_par[var]["std"] = np.std(l_value)
            self.d_par[var]["loading"] = sqrt(2) / 2.0
        
        # RSCC loading must be negative (opposite to RSR)
        self.d_par["rscc"]["loading"] = -self.d_par["rscc"]["loading"]
        
        logging.debug(f"PDB archive parameters - RSR: mean={self.d_par['rsr']['mean']:.3f}, std={self.d_par['rsr']['std']:.3f}")
        logging.debug(f"PDB archive parameters - RSCC: mean={self.d_par['rscc']['mean']:.3f}, std={self.d_par['rscc']['std']:.3f}")
    
    def _get_ranking_percentile(self, score: float, ref_list: list) -> float:
        """
        Calculate ranking percentile of a score within a reference list.
        
        Args:
            score: PC1 score to rank
            ref_list: Reference list of PC1 scores from PDB archive
            
        Returns:
            Percentile (0-1) where 1.0 = best quality, 0.0 = worst quality
        """
        # Create a copy and add new score
        ref_copy = ref_list.copy()
        ref_copy.append(score)
        
        # Sort the list
        sorted_list = sorted(ref_copy)
        
        # Count how many values are less than this score (handles ties correctly)
        rank = sum(1 for x in sorted_list if x < score)
        
        # Calculate percentile
        if len(ref_copy) > 1:
            percentile = rank / float(len(ref_copy) - 1)
        else:
            percentile = 0.5
        
        # Invert so best (highest rank) = 1.0
        return 1 - percentile
    
    def calculate_fit_score(self, ligand_cc: float, ligand_rsr: float) -> Optional[float]:
        """
        Calculate fit quality percentile score vs PDB archive.
        
        Args:
            ligand_cc: Ligand correlation coefficient (equivalent to RSCC)
            ligand_rsr: Ligand real space R-factor
            
        Returns:
            Percentile score (0-100) or None if calculation fails
        """
        if not self.loaded:
            return None
        
        try:
            # Calculate fit PC1
            fit_pc1 = ((ligand_rsr - self.d_par["rsr"]["mean"]) / self.d_par["rsr"]["std"]) * \
                      self.d_par["rsr"]["loading"] + \
                      ((ligand_cc - self.d_par["rscc"]["mean"]) / self.d_par["rscc"]["std"]) * \
                      self.d_par["rscc"]["loading"]
            
            # Get percentile vs PDB archive
            percentile = self._get_ranking_percentile(fit_pc1, self.d_ref["fit_pc1"])
            
            return percentile * 100  # Convert to 0-100 scale
            
        except (KeyError, ZeroDivisionError, TypeError, ValueError) as e:
            logging.debug(f"Failed to calculate PDB fit score: {e}")
            return None
    
    def calculate_geometry_score(self, mogul_z_bond: float, mogul_z_angle: float) -> Optional[float]:
        """
        Calculate geometry quality percentile score vs PDB archive.
        
        Args:
            mogul_z_bond: Mogul Z-score for bonds
            mogul_z_angle: Mogul Z-score for angles
            
        Returns:
            Percentile score (0-100) or None if calculation fails
        """
        if not self.loaded:
            return None
        
        try:
            # Calculate geometry PC1
            geo_pc1 = ((mogul_z_bond - self.d_par["mogul_bonds_rmsz"]["mean"]) / 
                       self.d_par["mogul_bonds_rmsz"]["std"]) * \
                      self.d_par["mogul_bonds_rmsz"]["loading"] + \
                      ((mogul_z_angle - self.d_par["mogul_angles_rmsz"]["mean"]) / 
                       self.d_par["mogul_angles_rmsz"]["std"]) * \
                      self.d_par["mogul_angles_rmsz"]["loading"]
            
            # Get percentile vs PDB archive
            percentile = self._get_ranking_percentile(geo_pc1, self.d_ref["geo_pc1"])
            
            return percentile * 100  # Convert to 0-100 scale
            
        except (KeyError, ZeroDivisionError, TypeError, ValueError) as e:
            logging.debug(f"Failed to calculate PDB geometry score: {e}")
            return None

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
    """
    Compare stereochemistry between input and output structures.
    
    Prefers CIF files for topology (authoritative stereochemistry from refinement restraints),
    but falls back to SMILES template approach when CIF files are unavailable.
    
    Args:
        input_cif: Path to input CIF file (or 'NA' if unavailable)
        input_pdb: Path to input PDB file (or 'NA' if unavailable)
        output_cif: Path to output CIF file (or 'NA' if unavailable)
        output_pdb: Path to output PDB file
        smiles_file: Path to input SMILES file
        dataset_name: Name of dataset for logging
    
    Returns:
        Tuple of (chirality_result, output_smiles, heavy_atom_count)
        Result will be prefixed with "[CIF]" or "[SMILES fallback]" to indicate method used
    """
    # Check which files are available
    has_input_cif = input_cif and input_cif != 'NA' and os.path.isfile(input_cif)
    has_output_cif = output_cif and output_cif != 'NA' and os.path.isfile(output_cif)
    has_smiles = smiles_file and smiles_file != 'NA' and os.path.isfile(smiles_file)
    has_output_pdb = output_pdb and output_pdb != 'NA' and os.path.isfile(output_pdb)
    
    if not has_output_pdb:
        return "Output PDB not found", "NA", 0
    if not has_smiles:
        return "Input SMILES not found", "NA", 0
    
    # Determine which method to use
    use_cif_method = has_input_cif and has_output_cif
    
    try:
        # Load reference SMILES
        with open(smiles_file, 'r', encoding='utf-8') as f:
            smiles_line = f.readline().strip()
        smiles = smiles_line.split()[0] if smiles_line else ""
        
        if not smiles:
            return f"Empty or invalid SMILES file: {smiles_file}", "NA", 0
        
        # ============ METHOD 1: CIF TOPOLOGY (PREFERRED) ============
        if use_cif_method:
            logging.debug(f"Using CIF topology method for {dataset_name}")
            
            # Process input molecule
            logging.debug(f"Loading input CIF: {input_cif}")
            mol_input = load_mol_from_cif(input_cif)
            
            # Assign coordinates from input PDB if available, otherwise use CIF coordinates
            if input_pdb and input_pdb != 'NA' and os.path.isfile(input_pdb):
                logging.debug(f"Assigning coordinates from input PDB: {input_pdb}")
                mol_input = assign_coordinates_from_pdb(mol_input, input_pdb, dataset_name, "input")
            else:
                logging.debug(f"No input PDB available, using CIF coordinates for input molecule")
                Chem.AssignStereochemistryFrom3D(mol_input)

            # Process output molecule - use CIF topology + PDB coordinates
            logging.debug(f"Loading output CIF: {output_cif}")
            mol_output = load_mol_from_cif(output_cif)
            logging.debug(f"Assigning coordinates from output PDB: {output_pdb}")
            mol_output = assign_coordinates_from_pdb(mol_output, output_pdb, dataset_name, "output")

            # Extract heavy atom count
            heavy_atoms = get_heavy_atom_count(mol_output)
            logging.info(f"Ligand contains {heavy_atoms} atoms")

            # Find chiral centers
            chiral_input = Chem.FindMolChiralCenters(mol_input, includeUnassigned=True)
            chiral_output = Chem.FindMolChiralCenters(mol_output, includeUnassigned=True)

            # Generate output SMILES
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
                    differences.append(f"Atom {idx}: {input_config} -> {output_config}")

            if differences:
                result = f"Chiral centre inverted: {'; '.join(differences)} [CIF]"
            else:
                result = "Stereochemistry matches input [CIF]"
                
            logging.debug(f"CIF method - Input chiral centers: {chiral_input}")
            logging.debug(f"CIF method - Output chiral centers: {chiral_output}")
            
            return result, smiles_output, heavy_atoms
        
        # ============ METHOD 2: SMILES TEMPLATE (FALLBACK) ============
        else:
            missing_cif = []
            if not has_input_cif:
                missing_cif.append("input CIF")
            if not has_output_cif:
                missing_cif.append("output CIF")
            
            logging.debug(f"CIF files missing ({', '.join(missing_cif)}), using SMILES template fallback for {dataset_name}")
            
            # Load input from SMILES
            mol_input = Chem.MolFromSmiles(smiles)
            if mol_input is None:
                return f"[SMILES fallback] Failed to parse SMILES: {smiles}", "NA", 0
            
            # Read output structure from PDB
            mol_output = Chem.MolFromPDBFile(output_pdb, removeHs=False)
            if mol_output is None:
                return f"[SMILES fallback] Failed to parse output PDB: {output_pdb}", "NA", 0
            
            # Extract only the ligand from output PDB
            mol_output = _extract_ligand_from_pdb(mol_output)
            
            # Use SMILES as template to assign bond orders (PDB lacks this info)
            try:
                template_no_h = Chem.RemoveHs(mol_input)
                output_no_h = Chem.RemoveHs(mol_output)
                mol_output_corrected = AllChem.AssignBondOrdersFromTemplate(template_no_h, output_no_h)
                
                if mol_output_corrected is not None:
                    mol_output = mol_output_corrected
                    logging.debug("Successfully assigned bond orders using SMILES template")
                else:
                    logging.warning("Template assignment failed, using standard sanitization")
                    Chem.SanitizeMol(mol_output)
            except Exception as e:
                logging.warning(f"Template approach failed: {e}, using standard sanitization")
                try:
                    Chem.SanitizeMol(mol_output)
                except Exception as e2:
                    logging.warning(f"Sanitization failed: {e2}")
            
            # Assign stereochemistry
            Chem.AssignStereochemistry(mol_output, force=True, cleanIt=True)
            
            # Extract atom count
            heavy_atoms = get_heavy_atom_count(mol_output)
            logging.info(f"Ligand contains {heavy_atoms} atoms")
            
            # Generate output SMILES
            try:
                mol_output_no_h = Chem.RemoveHs(mol_output)
                smiles_output = Chem.MolToSmiles(mol_output_no_h, isomericSmiles=True, canonical=True)
            except Exception as e:
                logging.warning(f"Failed to generate output SMILES: {e}")
                smiles_output = "Failed to generate"
            
            # Find chiral centers
            input_chiral_centers = _get_chiral_centers(mol_input)
            output_chiral_centers = _get_chiral_centers(mol_output)
            
            # Align molecules and map chiral centers (needed because atom indices differ)
            try:
                mapped_output_chiral = _align_molecules_and_map_chiral_centers(
                    mol_input, mol_output, input_chiral_centers, output_chiral_centers, "SMILES"
                )
                
                input_dict = {idx: config for idx, config in input_chiral_centers}
                output_dict = mapped_output_chiral
                
            except Exception as e:
                return f"[SMILES fallback] Error during structural alignment: {e}", smiles_output, heavy_atoms
            
            # Compare chiral centers
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
                result = f"Chiral centre inverted: {'; '.join(differences)} [SMILES fallback]"
            else:
                result = "Stereochemistry matches input [SMILES fallback]"
            
            logging.debug(f"SMILES fallback method - Input chiral centers: {len(input_chiral_centers)}")
            logging.debug(f"SMILES fallback method - Output chiral centers: {len(output_chiral_centers)}")
            
            return result, smiles_output, heavy_atoms
        
    except FileNotFoundError as e:
        logging.error(f"File not found in chirality comparison: {e}")
        return f"File not found: {e}", "NA", 0
    except (IOError, OSError) as e:
        logging.error(f"File I/O error in chirality comparison: {e}")
        return f"File error: {e}", "NA", 0
    except ValueError as e:
        logging.error(f"Value error in chirality comparison: {e}")
        logging.debug(f"Full traceback: {traceback.format_exc()}")
        return f"Structure error: {e}", "NA", 0
    except Exception as e:
        logging.error(f"Unexpected error in chirality comparison: {e}")
        logging.debug(f"Full traceback: {traceback.format_exc()}")
        return f"Error: {str(e)}", "NA", 0

# Column order for consistent output across HTML and JSON
COLUMN_ORDER = [
    'Export to XCE',
    'Comments',
    'Crystal Name',
    'Ligand ID',
    'Compound Code',
    'Ligand Structure',
    'Chiral inversion',
    'Ligand Density',
    'Ligand Report HTML',
    'Composite Quality Score',
    'Spider Plot',
    'PCA Percentile Plot (Batch)',
    'PCA Percentile Plot (PDB)',
    'Ligand CC',
    'Ligand RSR',
    'Ligand Fit Quality Score',
    'Ligand Fit PCA Percentile (Batch)',
    'Ligand Fit PCA Percentile (PDB)',
    'Mogul Z Angle',
    'Mogul Z Bond',
    'Ligand Geometry Quality Score',
    'Ligand Geometry PCA Percentile (Batch)',
    'Ligand Geometry PCA Percentile (PDB)',
    'Ligand Occupancy',
    'Ligand Mean B Factor',
    'Ligand Atoms',
    'Ligand Clashes',
    'Ligand Clashes per Atom',
    'Buster Report HTML',
    'High Resolution (Å)',
    'Low Resolution (Å)',
    'R',
    'Rfree',
    'Mol Probity Score',
    'Mean B Factor',
    'Clash Score',
    'Rama Favored Percent',
    'Rama Outlier Percent',
    'C-beta Deviations',
    'Poor Rotamer Percent',
    'RMS Bonds',
    'RMS Angles',
    'Pipedream Directory',
    # 'Pipedream Summary',
    # 'PDB File',
    # 'MTZ File',
    # '2Fo-Fc Map File',
    # 'Fo-Fc Map File',
    # 'Output SMILES'
]

# Column tooltips for HTML table headers
COLUMN_TOOLTIPS = {
    'Crystal Name': 'Purple border indicates dataset contains multiple ligands',
    'Compound Code': 'Green=best quality, yellow=middle, red=worst quality; purple border indicates duplicate compound.',
    'Chiral inversion': 'Has stereochemistry changed from input SMILES?',
    'Ligand CC': 'Green ≥0.95 (very good), Light green ≥0.90 (good), Yellow ≥0.80 (OK), Red <0.80 (poor)',
    'Ligand RSR': 'Green <0.2 (very good), Light green <0.3 (good), Yellow <0.4 (OK), Red ≥0.4 (poor)',
    'Ligand Fit Quality Score': 'Combined score (0-100) from Ligand CC and RSR using piecewise linear conversion. CC scale emphasizes high quality (flat 100 for CC≥0.95, steep slope 0.90-0.95). RSR scale inverted (lower is better) with progressive penalties. Average of both component scores. Green ≥75 (excellent), Yellow ≥50 (good/OK), Red <50 (poor)',
    'Ligand Geometry Quality Score': 'Combined score (0-100) from Mogul Z-scores (bonds and angles) using piecewise linear conversion. Both Z-scores converted with thresholds: <1.5 (good, 80-100 pts), 1.5-2.5 (OK, 50-80 pts), 2.5-4.0 (poor, 20-50 pts), >4.0 (bad, 0-20 pts). Average of both component scores. Green ≥75 (excellent), Yellow ≥50 (good/OK), Red <50 (poor)',
    'Composite Quality Score': 'Weighted average (60% Fit, 40% Geometry) with penalties for high B-factor ratio (>1.2) and ligand clashes (≥1). Green ≥75 (excellent), Yellow ≥50 (good/OK), Red <50 (poor)',
    'Ligand Fit PCA Percentile (Batch)': 'Percentile rank (0-100%) within this batch based on PCA of Ligand CC and RSR. Score represents percentage of batch with worse fit quality. Uses z-score normalization and PC1 projection, then rank-based percentile calculation. 100% = best fit in batch, 0% = worst fit. Green ≥75%, Yellow ≥50%, Red <50%',
    'Ligand Fit PCA Percentile (PDB)': 'Percentile rank (0-100%) vs PDB based on PCA of Ligand CC and RSR. Score represents percentage of all PDB structures with worse fit quality. Uses 2020 PDB statistics (mean/std) for z-score normalization and PC1 projection. 100% = better than all PDB, 0% = worse than all PDB. Green ≥75%, Yellow ≥50%, Red <50%',
    'Ligand Geometry PCA Percentile (Batch)': 'Percentile rank (0-100%) within this batch based on PCA of Mogul Z-scores. Score represents percentage of batch with worse geometry. Uses z-score normalization and PC1 projection, then rank-based percentile calculation. 100% = best geometry in batch, 0% = worst. Green ≥75%, Yellow ≥50%, Red <50%',
    'Ligand Geometry PCA Percentile (PDB)': 'Percentile rank (0-100%) vs PDB based on PCA of Mogul Z-scores. Score represents percentage of all PDB structures with worse geometry. Uses 2020 PDB statistics for z-score normalization and PC1 projection. 100% = better than all PDB, 0% = worse than all PDB. Green ≥75%, Yellow ≥50%, Red <50%',
    'PCA Percentile Plot (Batch)': 'Interactive 2D scatter plot showing dataset position within this batch. X-axis: Ligand Fit PCA Percentile (Batch), Y-axis: Ligand Geometry PCA Percentile (Batch). Colored quadrants indicate relative quality within batch: green (excellent), yellow (mixed), red (poor). Reference lines at 50% and 75%. Click to view interactive plot with hover details',
    'PCA Percentile Plot (PDB)': 'Interactive 2D scatter plot showing dataset position vs PDB. X-axis: Ligand Fit PCA Percentile (PDB), Y-axis: Ligand Geometry PCA Percentile (PDB). Colored quadrants indicate quality regions: green (excellent), yellow (mixed), red (poor). Reference lines at 50% and 75%. Click to view interactive plot with hover details comparing to PDB structures',
    'Spider Plot': 'Interactive radar/spider plot showing 6 key quality metrics on a 0-100 scale: Ligand CC (higher better), Ligand Fit/RSR (inverted, lower better), Mogul Z Angle (inverted), Mogul Z Bond (inverted), Ligand Clashes (inverted), B-factor Ratio (closer to 1.0 better). All metrics normalized to 0-100 where 100=best quality. Click to view interactive plot with hover details for each metric',
    'Mogul Z Angle': 'Green <1.5 (good), Yellow 1.5-2.5 (OK), Light purple 2.5-4 (poor), Purple >4 (bad)',
    'Mogul Z Bond': 'Green <1.5 (good), Yellow 1.5-2.5 (OK), Light purple 2.5-4 (poor), Purple >4 (bad)',
    'Ligand Mean B Factor': f'Red indicates ligand B-factor ratio ≥{THRESHOLDS.b_factor_ratio} times mean B-factor',
    'Ligand Clashes': 'Number of ligand clashes (>0.4 Å overlap) from MolProbity',
    'Ligand Clashes per Atom': 'Ligand clashes normalized by total ligand atom count',
    'R': f'Red indicates R-factor >{THRESHOLDS.r_max}',
    'Rfree': f'Red indicates R-free >{THRESHOLDS.rfree_max}',
    'Rama Favored Percent': 'Green ≥98% (excellent Ramachandran geometry)',
    'Rama Outlier Percent': f'Red indicates Ramachandran outliers >{THRESHOLDS.rama_outlier_max}%',
    'Poor Rotamer Percent': f'Red indicates poor rotamers >{THRESHOLDS.poor_rotamer_max}%',
    'Clash Score': f'Red indicates clash score >{THRESHOLDS.clash_score_max}',
    'Mol Probity Score': f'Red indicates MolProbity score >{THRESHOLDS.molprobity_score_max}',
    'RMS Angles': f'Red indicates RMS angles >{THRESHOLDS.rms_angles_max}°',
    'RMS Bonds': f'Red indicates RMS bonds >{THRESHOLDS.rms_bonds_max} Å',
}


def calculate_quality_scores(ligand_cc: Any, ligand_rsr: Any, 
                            mogul_z_bond: Any, mogul_z_angle: Any,
                            ligand_mean_b: Any = None, global_mean_b: Any = None,
                            ligand_clashes: Any = None) -> Dict[str, Any]:
    """
    Calculate fit quality, geometry quality, and composite quality scores.
    
    Args:
        ligand_cc: Ligand correlation coefficient
        ligand_rsr: Ligand real space R-factor
        mogul_z_bond: Mogul Z score for bonds
        mogul_z_angle: Mogul Z score for angles
        ligand_mean_b: Mean ligand B-factor (optional, for composite score)
        global_mean_b: Global mean B-factor (optional, for composite score)
        ligand_clashes: Number of ligand clashes (optional, for composite score)
        
    Returns:
        Dictionary with 'Ligand Fit Quality Score', 'Ligand Geometry Quality Score', 'Composite Quality Score'
    """
    # Convert to numeric, handling 'NA' values
    try:
        cc = float(ligand_cc) if ligand_cc not in ['NA', None] else None
    except (ValueError, TypeError):
        cc = None
    
    try:
        rsr = float(ligand_rsr) if ligand_rsr not in ['NA', None] else None
    except (ValueError, TypeError):
        rsr = None
    
    try:
        z_bond = float(mogul_z_bond) if mogul_z_bond not in ['NA', None] else None
    except (ValueError, TypeError):
        z_bond = None
    
    try:
        z_angle = float(mogul_z_angle) if mogul_z_angle not in ['NA', None] else None
    except (ValueError, TypeError):
        z_angle = None
    
    # Calculate CC score
    def score_ligand_cc(cc_val):
        if cc_val >= 0.95:
            return 100
        elif cc_val >= 0.90:
            return 90 + ((cc_val - 0.90) / 0.05) * 10
        elif cc_val >= 0.80:
            return 50 + ((cc_val - 0.80) / 0.10) * 40
        else:
            return max(((cc_val - 0.60) / 0.20) * 50, 0)
    
    # Calculate RSR score
    def score_rsr(rsr_val):
        if rsr_val < 0.2:
            return 100 - (rsr_val / 0.2) * 10
        elif rsr_val < 0.3:
            return 90 - ((rsr_val - 0.2) / 0.1) * 20
        elif rsr_val < 0.4:
            return 70 - ((rsr_val - 0.3) / 0.1) * 30
        else:
            return max(40 - ((rsr_val - 0.4) / 0.2) * 40, 0)
    
    # Calculate Mogul Z score
    def score_mogul_z(z_val):
        if z_val < THRESHOLDS.mogul_z_good:
            return 100 - (z_val / THRESHOLDS.mogul_z_good) * 20
        elif z_val < THRESHOLDS.mogul_z_ok:
            return 80 - ((z_val - THRESHOLDS.mogul_z_good) / 
                          (THRESHOLDS.mogul_z_ok - THRESHOLDS.mogul_z_good) * 30)
        elif z_val < THRESHOLDS.mogul_z_poor:
            return 50 - ((z_val - THRESHOLDS.mogul_z_ok) / 
                         (THRESHOLDS.mogul_z_poor - THRESHOLDS.mogul_z_ok) * 30)
        else:
            return max(20 - ((z_val - THRESHOLDS.mogul_z_poor) / 2.0 * 20), 0)
    
    # Calculate component scores
    cc_score = score_ligand_cc(cc) if cc is not None else None
    rsr_score = score_rsr(rsr) if rsr is not None else None
    bond_score = score_mogul_z(z_bond) if z_bond is not None else None
    angle_score = score_mogul_z(z_angle) if z_angle is not None else None
    
    # Calculate fit quality (average of CC and RSR scores)
    if cc_score is not None and rsr_score is not None:
        fit_quality = (cc_score + rsr_score) / 2
    elif cc_score is not None:
        fit_quality = cc_score
    elif rsr_score is not None:
        fit_quality = rsr_score
    else:
        fit_quality = 'NA'
    
    # Calculate geometry quality (average of bond and angle scores)
    if bond_score is not None and angle_score is not None:
        geometry_quality = (bond_score + angle_score) / 2
    elif bond_score is not None:
        geometry_quality = bond_score
    elif angle_score is not None:
        geometry_quality = angle_score
    else:
        geometry_quality = 'NA'
    
    # Calculate base overall quality (weighted average: 60% fit, 40% geometry)
    if fit_quality != 'NA' and geometry_quality != 'NA':
        base_quality = (fit_quality * 0.6) + (geometry_quality * 0.4)
    elif fit_quality != 'NA':
        base_quality = fit_quality
    elif geometry_quality != 'NA':
        base_quality = geometry_quality
    else:
        base_quality = 'NA'
    
    # Calculate composite quality score with contextual penalties
    composite_quality = base_quality
    if composite_quality != 'NA':
        total_penalty = 0
        
        # B-factor ratio penalty (max 10)
        try:
            ligand_b = float(ligand_mean_b) if ligand_mean_b not in ['NA', None] else None
            global_b = float(global_mean_b) if global_mean_b not in ['NA', None] else None
            
            if ligand_b is not None and global_b is not None and global_b > 0:
                b_factor_ratio = ligand_b / global_b
                if b_factor_ratio > 1.5:
                    total_penalty += 10  # Bad - likely poorly ordered or artifactual
                elif b_factor_ratio > 1.2:
                    total_penalty += 5  # Caution - ligand more mobile
                # else: ratio <= 1.2, no penalty (OK or perfect)
        except (ValueError, TypeError):
            pass  # Skip penalty if values invalid
        
        # Ligand clashes penalty (max 10)
        try:
            clashes = int(ligand_clashes) if ligand_clashes not in ['NA', None] else None
            
            if clashes is not None and clashes >= 0:
                if clashes > 3:
                    total_penalty += 10  # Bad - serious stereochemical problem
                elif clashes >= 1:
                    total_penalty += 5  # Caution - some clashes present
                # else: clashes == 0, no penalty (good)
        except (ValueError, TypeError):
            pass  # Skip penalty if value invalid
        
        # Apply penalties (max combined penalty of 40)
        composite_quality = max(composite_quality - total_penalty, 0)
    
    return {
        'Ligand Fit Quality Score': safe_round(fit_quality, 1) if fit_quality != 'NA' else 'NA',
        'Ligand Geometry Quality Score': safe_round(geometry_quality, 1) if geometry_quality != 'NA' else 'NA',
        'Composite Quality Score': safe_round(composite_quality, 1) if composite_quality != 'NA' else 'NA'
    }


def calculate_pca_scores_batch(results: List[Dict[str, Any]], pca_ref: Optional[PCAReference] = None) -> List[Dict[str, Any]]:
    """
    Calculate PCA-based scores for fit and geometry quality relative to the current batch.
    This shows where each dataset sits within the group being analyzed.
    
    Args:
        results: List of result dictionaries with quality metrics
        pca_ref: Optional PCAReference object with PDB archive parameters for consistent PC1 calculation
        
    Returns:
        Updated results list with added PCA score columns
    """
    import numpy as np
    from math import sqrt
    
    # Extract fit metrics (Ligand CC and RSR)
    fit_data = []
    fit_indices = []
    for i, result in enumerate(results):
        cc = result.get('Ligand CC', 'NA')
        rsr = result.get('Ligand RSR', 'NA')
        try:
            cc_val = float(cc) if cc not in ['NA', None] else None
            rsr_val = float(rsr) if rsr not in ['NA', None] else None
            if cc_val is not None and rsr_val is not None:
                fit_data.append([rsr_val, cc_val])  # Note: [RSR, CC] order
                fit_indices.append(i)
        except (ValueError, TypeError):
            pass
    
    # Calculate fit PCA scores if we have enough data
    if len(fit_data) >= 2:
        fit_array = np.array(fit_data)
        
        # Use PDB parameters if available for consistent PC1 calculation
        # This ensures batch and PDB rankings have the same relative order
        if pca_ref is not None and pca_ref.loaded:
            mean_rsr = pca_ref.d_par["rsr"]["mean"]
            std_rsr = pca_ref.d_par["rsr"]["std"]
            loading_rsr = pca_ref.d_par["rsr"]["loading"]
            mean_cc = pca_ref.d_par["rscc"]["mean"]
            std_cc = pca_ref.d_par["rscc"]["std"]
            loading_cc = pca_ref.d_par["rscc"]["loading"]
        else:
            # Fallback to batch statistics if PDB ref not available
            mean_rsr = np.mean(fit_array[:, 0])
            std_rsr = np.std(fit_array[:, 0], ddof=1) if len(fit_array) > 1 else 1.0
            mean_cc = np.mean(fit_array[:, 1])
            std_cc = np.std(fit_array[:, 1], ddof=1) if len(fit_array) > 1 else 1.0
            
            # Avoid division by zero
            if std_rsr == 0:
                std_rsr = 1.0
            if std_cc == 0:
                std_cc = 1.0
            
            # PCA loadings (equal magnitude, opposite signs)
            loading_rsr = sqrt(2) / 2.0
            loading_cc = -sqrt(2) / 2.0
        
        # Calculate PC1 for each dataset using consistent parameters
        pc1_values = []
        for idx, data_idx in enumerate(fit_indices):
            rsr = fit_array[idx, 0]
            cc = fit_array[idx, 1]
            
            # Z-score normalization and PCA projection
            fit_pc1 = ((rsr - mean_rsr) / std_rsr) * loading_rsr + \
                     ((cc - mean_cc) / std_cc) * loading_cc
            
            pc1_values.append((data_idx, fit_pc1))
        
        # Calculate percentile ranking within batch
        # Lower PC1 = better quality (low RSR, high CC)
        if len(pc1_values) > 0:
            # Extract all PC1 values from the batch
            all_pc1_list = [pc1 for _, pc1 in pc1_values]
            
            for data_idx, fit_pc1 in pc1_values:
                # Sort all PC1 values in the batch
                sorted_list = sorted(all_pc1_list)
                
                # Count how many values are STRICTLY less than this one
                rank = sum(1 for x in sorted_list if x < fit_pc1)
                
                # Calculate percentile (0-1 range)
                if len(all_pc1_list) > 1:
                    percentile = rank / float(len(all_pc1_list) - 1)
                else:
                    percentile = 0.5  # Single structure gets middle score
                
                # Invert so lowest PC1 (best quality) gets highest percentile
                fit_pca_score = (1 - percentile) * 100
                results[data_idx]['Ligand Fit PCA Percentile (Batch)'] = round(fit_pca_score, 1)
    
    # Extract geometry metrics (Mogul Z scores)
    geom_data = []
    geom_indices = []
    for i, result in enumerate(results):
        z_bond = result.get('Mogul Z Bond', 'NA')
        z_angle = result.get('Mogul Z Angle', 'NA')
        try:
            bond_val = float(z_bond) if z_bond not in ['NA', None] else None
            angle_val = float(z_angle) if z_angle not in ['NA', None] else None
            if bond_val is not None and angle_val is not None:
                geom_data.append([bond_val, angle_val])
                geom_indices.append(i)
        except (ValueError, TypeError):
            pass
    
    # Calculate geometry PCA scores if we have enough data
    if len(geom_data) >= 2:
        geom_array = np.array(geom_data)
        
        # Use PDB parameters if available for consistent PC1 calculation
        # This ensures batch and PDB rankings have the same relative order
        if pca_ref is not None and pca_ref.loaded:
            mean_z_bond = pca_ref.d_par["mogul_bonds_rmsz"]["mean"]
            std_z_bond = pca_ref.d_par["mogul_bonds_rmsz"]["std"]
            loading_z_bond = pca_ref.d_par["mogul_bonds_rmsz"]["loading"]
            mean_z_angle = pca_ref.d_par["mogul_angles_rmsz"]["mean"]
            std_z_angle = pca_ref.d_par["mogul_angles_rmsz"]["std"]
            loading_z_angle = pca_ref.d_par["mogul_angles_rmsz"]["loading"]
        else:
            # Fallback to batch statistics if PDB ref not available
            mean_z_bond = np.mean(geom_array[:, 0])
            std_z_bond = np.std(geom_array[:, 0], ddof=1) if len(geom_array) > 1 else 1.0
            mean_z_angle = np.mean(geom_array[:, 1])
            std_z_angle = np.std(geom_array[:, 1], ddof=1) if len(geom_array) > 1 else 1.0
            
            # Avoid division by zero
            if std_z_bond == 0:
                std_z_bond = 1.0
            if std_z_angle == 0:
                std_z_angle = 1.0
            
            # PCA loadings (both positive for Mogul Z scores - higher is worse)
            loading_z_bond = sqrt(2) / 2.0
            loading_z_angle = sqrt(2) / 2.0
        
        # Calculate PC1 for each dataset using consistent parameters
        pc1_values = []
        for idx, data_idx in enumerate(geom_indices):
            z_bond = geom_array[idx, 0]
            z_angle = geom_array[idx, 1]
            
            # Z-score normalization and PCA projection
            geom_pc1 = ((z_bond - mean_z_bond) / std_z_bond) * loading_z_bond + \
                      ((z_angle - mean_z_angle) / std_z_angle) * loading_z_angle
            
            pc1_values.append((data_idx, geom_pc1))
        
        # Calculate percentile ranking within batch
        # Lower PC1 = better geometry (lower Mogul Z scores)
        if len(pc1_values) > 0:
            # Extract all PC1 values from the batch
            all_pc1_list = [pc1 for _, pc1 in pc1_values]
            
            for data_idx, geom_pc1 in pc1_values:
                # Sort all PC1 values in the batch
                sorted_list = sorted(all_pc1_list)
                
                # Count how many values are STRICTLY less than this one
                rank = sum(1 for x in sorted_list if x < geom_pc1)
                
                # Calculate percentile (0-1 range)
                if len(all_pc1_list) > 1:
                    percentile = rank / float(len(all_pc1_list) - 1)
                else:
                    percentile = 0.5  # Single structure gets middle score
                
                # Invert so lowest PC1 (best quality) gets highest percentile
                geom_pca_score = (1 - percentile) * 100
                results[data_idx]['Ligand Geometry PCA Percentile (Batch)'] = round(geom_pca_score, 1)
    
    return results


def calculate_pca_scores_pdb(results: List[Dict[str, Any]], pca_ref: Optional[PCAReference] = None) -> List[Dict[str, Any]]:
    """
    Calculate PCA-based scores relative to the entire PDB archive.
    This shows where each dataset sits compared to all known structures.
    
    Args:
        results: List of result dictionaries with quality metrics
        pca_ref: PCAReference object with loaded PDB archive data
        
    Returns:
        Updated results list with added PCA score columns
    """
    if pca_ref is None or not pca_ref.loaded:
        logging.debug("PDB PCA reference not available, skipping PDB archive scores")
        # Set all PDB scores to 'NA'
        for result in results:
            result['Ligand Fit PCA Percentile (PDB)'] = 'NA'
            result['Ligand Geometry PCA Percentile (PDB)'] = 'NA'
        return results
    
    logging.info(f"Calculating PDB archive percentile scores for {len(results)} structures")
    
    # Calculate fit scores against PDB archive
    for result in results:
        ligand_cc = result.get('Ligand CC', 'NA')
        ligand_rsr = result.get('Ligand RSR', 'NA')
        
        try:
            cc_val = float(ligand_cc) if ligand_cc not in ['NA', None] else None
            rsr_val = float(ligand_rsr) if ligand_rsr not in ['NA', None] else None
            
            if cc_val is not None and rsr_val is not None:
                fit_score = pca_ref.calculate_fit_score(cc_val, rsr_val)
                result['Ligand Fit PCA Percentile (PDB)'] = round(fit_score, 1) if fit_score is not None else 'NA'
            else:
                result['Ligand Fit PCA Percentile (PDB)'] = 'NA'
        except (ValueError, TypeError):
            result['Ligand Fit PCA Percentile (PDB)'] = 'NA'
    
    # Calculate geometry scores against PDB archive
    for result in results:
        mogul_z_bond = result.get('Mogul Z Bond', 'NA')
        mogul_z_angle = result.get('Mogul Z Angle', 'NA')
        
        try:
            bond_val = float(mogul_z_bond) if mogul_z_bond not in ['NA', None] else None
            angle_val = float(mogul_z_angle) if mogul_z_angle not in ['NA', None] else None
            
            if bond_val is not None and angle_val is not None:
                geom_score = pca_ref.calculate_geometry_score(bond_val, angle_val)
                result['Ligand Geometry PCA Percentile (PDB)'] = round(geom_score, 1) if geom_score is not None else 'NA'
            else:
                result['Ligand Geometry PCA Percentile (PDB)'] = 'NA'
        except (ValueError, TypeError):
            result['Ligand Geometry PCA Percentile (PDB)'] = 'NA'
    
    return results


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


def get_heavy_atom_count(mol: Chem.Mol) -> int:
    """
    Get the total number of atoms in a molecule.
    
    Args:
        mol: RDKit molecule object
        
    Returns:
        Total number of atoms (including hydrogens), or 0 if molecule is invalid
    """
    if mol is None:
        return 0
    
    try:
        return mol.GetNumAtoms()
    except (ValueError, RuntimeError, AttributeError) as e:
        logging.warning(f"Error counting atoms: {e}")
        return 0

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
        
        # Check write permissions before attempting
        output_dir = os.path.dirname(output_path)
        if not os.access(output_dir, os.W_OK):
            logging.error(f"No write permission for directory: {output_dir}")
            return f"Error: No write permission"
            
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
    try:
        from rdkit.Chem import rdFMCS
    except ImportError as e:
        raise ImportError(f"Failed to import rdFMCS: {e}")
    
    try:
        mcs = rdFMCS.FindMCS([input_mol, output_mol], 
                           bondCompare=rdFMCS.BondCompare.CompareOrderExact,
                           atomCompare=rdFMCS.AtomCompare.CompareElements,
                           timeout=PROCESSING.mcs_timeout_seconds)
    except RuntimeError as e:
        raise ValueError(f"MCS calculation failed or timed out: {e}")
    
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


def safe_round(value: Any, digits: int = 3) -> Union[float, str]:
    """Round a value to a given number of digits, or return 'NA' if not possible."""
    try:
        num_val = float(value)
        if pd.isna(num_val):
            return 'NA'
        return round(num_val, digits)
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
    
    # Return button styled consistently with percentile plot buttons (white with blue border)
    return f'''<td><a href="#" class="btn btn-sm btn-outline-primary file-link" 
        data-file="file://{cell_val}" data-ext="{file_ext}" 
        style="padding:5px 10px;text-decoration:none;">View</a></td>'''


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
    
    def find_electron_density_gif(self, ligand_id: str = None) -> str:
        """
        Find the electron density movie GIF file.
        
        Args:
            ligand_id: Optional ligand ID (e.g., "LIG A4000") to find ligand-specific GIF
        """
        if not self.base_dir or not self.compound_code:
            return 'NA'
        
        gif_dir = self.base_dir / f"report-{self.compound_code}" / 'ligand' / 'pictures'
        if not gif_dir.exists():
            return 'NA'
        
        # Get all available GIFs
        all_gifs = list(gif_dir.glob(f"*{ELECTRON_DENSITY_GIF_SUFFIX}"))
        if not all_gifs:
            return 'NA'
        
        logging.debug(f"find_electron_density_gif called with ligand_id={ligand_id}")
        
        # If ligand_id is provided, try to find ligand-specific GIF
        # Expected format: "LIG A4000" -> filename "A4000_electrondensity_movie.gif"
        if ligand_id and ligand_id != 'NA':
            # Parse ligand ID: "LIG A4000" -> restype="LIG", chain_resnum="A4000"
            parts = ligand_id.split()
            logging.debug(f"Parsed ligand_id '{ligand_id}' into {len(parts)} parts: {parts}")
            if len(parts) == 2:
                restype, chain_resnum = parts
                # chain_resnum is like "A4000" where first char is chain, rest is resnum
                # We need the GIF pattern: {CHAIN}{RESNUM}_electrondensity_movie.gif = A4000_electrondensity_movie.gif
                
                # Primary pattern: {CHAIN_RESNUM}_electrondensity_movie.gif
                pattern = f"{chain_resnum}_{ELECTRON_DENSITY_GIF_SUFFIX}"
                gif_path = gif_dir / pattern
                logging.debug(f"Looking for ligand-specific GIF: {gif_path}")
                if gif_path.exists():
                    logging.info(f"Found ligand-specific density GIF for {ligand_id}: {gif_path}")
                    return str(gif_path)
                else:
                    logging.warning(f"Ligand-specific GIF not found for {ligand_id}: expected {gif_path}")
                    logging.info(f"Available GIFs in {gif_dir}: {[g.name for g in all_gifs]}")
        
        # Fall back to first GIF found (for single ligand or if specific not found)
        fallback_gif = str(all_gifs[0])
        logging.info(f"Using fallback density GIF: {fallback_gif}")
        return fallback_gif


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

def parse_edstats_for_ligand_rm(edstats_file: str, ligand_name: str = "LIG") -> Union[float, str]:
    """
    Parse edstats.out file to extract Rm value.
    
    Args:
        edstats_file: Path to the edstats.out file
        ligand_name: Name of the ligand to search for (default: "LIG")
        
    Returns:
        Rm value as float, or 'NA' if not found
    """
    if not os.path.exists(edstats_file):
        logging.debug(f"edstats file not found: {edstats_file}")
        return 'NA'
    
    try:
        with open(edstats_file, 'r', encoding='utf-8') as f:
            for line in f:
                # Skip empty lines and header lines
                line = line.strip()
                if not line or line.startswith('//') or line.startswith('RT'):
                    continue
                
                # Split the line into fields
                fields = line.split()
                if len(fields) < 6:
                    continue
                
                # First field is residue type (RT)
                residue_type = fields[0]
                
                if residue_type == ligand_name:
                    # Rm is the 5th column (index 4) in the format:
                    # RT  CI  RN     BAm  NPm   Rm    RGm  SRGm   CCSm   CCPm ...
                    try:
                        rm_value = float(fields[5])
                        logging.debug(f"Found Rm value for {ligand_name}: {rm_value}")
                        return rm_value
                    except (ValueError, IndexError) as e:
                        logging.warning(f"Could not parse Rm value from line: {line}, error: {e}")
                        continue
        
        logging.warning(f"No {ligand_name} residue found")
        return 'NA'  # Consistent string return instead of None
        
    except (IOError, OSError) as e:
        logging.error(f"Error reading edstats file {edstats_file}: {e}")
        return 'NA'
    except Exception as e:
        logging.error(f"Unexpected error parsing edstats file: {e}")
        logging.debug(f"Full traceback: {traceback.format_exc()}")
        return 'NA'


def parse_clash_out_for_ligand(clash_file: str, ligand_id: str = "LIG") -> int:
    """
    Parse clash.out file to count ligand clash instances.
    
    Args:
        clash_file: Path to the clash.out file
        ligand_id: Ligand ID to search for (e.g., "LIG" or "A4000 LIG")
        
    Returns:
        Number of clash instances involving the ligand, or 0 if none found
    """
    if not os.path.exists(clash_file):
        logging.debug(f"clash.out file not found: {clash_file}")
        return 0
    
    try:
        # Parse ligand_id to extract chain and residue info
        # ligand_id format: "LIG A4000" or just "LIG"
        ligand_parts = ligand_id.split()
        
        logging.debug(f"Parsing clash.out for ligand_id='{ligand_id}', parts={ligand_parts}")
        
        if len(ligand_parts) == 2:
            restype, chain_resnum = ligand_parts
        else:
            restype = ligand_id
            chain_resnum = None
        
        logging.debug(f"Looking for restype='{restype}', chain_resnum='{chain_resnum}'")
        
        clash_count = 0
        
        with open(clash_file, 'r', encoding='utf-8') as f:
            in_clash_section = False
            line_count = 0
            for line in f:
                line = line.strip()
                
                # Look for the clash section start
                if 'Bad Clashes' in line:
                    in_clash_section = True
                    logging.debug(f"Found 'Bad Clashes' section")
                    continue
                
                # Stop at clashscore line
                if line.startswith('clashscore'):
                    logging.debug(f"Reached clashscore line, stopping")
                    break
                
                # Parse clash lines
                if in_clash_section and line and not line.startswith('Using'):
                    line_count += 1
                    # Clash line format: " B 163  HIS  NE2  B4000  LIG  H1C :0.671"
                    # Note: For ligands, chain+resnum is combined in one field (e.g., "B4000")
                    # Split by whitespace
                    parts = line.split()
                    if len(parts) < 8:
                        continue
                    
                    # Check if either position contains the ligand
                    # Position 1: parts[0-3] (chain, resnum, restype, atom) - protein format
                    # Position 2: parts[4-6] (chain_resnum_combined, restype, atom) - ligand format
                    try:
                        # Check position 1 (protein residue format)
                        pos1_restype = parts[2]
                        pos1_chain_resnum = parts[0] + parts[1]  # e.g., "B163"
                        
                        # Check position 2 (ligand format - chain+resnum already combined)
                        pos2_chain_resnum = parts[4]  # e.g., "B4000" - already combined!
                        pos2_restype = parts[5]
                        
                        # Match against our ligand
                        match = False
                        if chain_resnum:
                            # Match specific chain+resnum (e.g., "A4000 LIG")
                            if (pos1_restype == restype and pos1_chain_resnum == chain_resnum) or \
                               (pos2_restype == restype and pos2_chain_resnum == chain_resnum):
                                match = True
                        else:
                            # Match just residue type (e.g., "LIG")
                            if pos1_restype == restype or pos2_restype == restype:
                                match = True
                        
                        if match:
                            clash_count += 1
                            logging.debug(f"Found ligand clash: {line}")
                    except (IndexError, ValueError) as e:
                        logging.debug(f"Could not parse clash line: {line}, error: {e}")
                        continue
            
            logging.debug(f"Processed {line_count} clash lines in total")
        
        if clash_count > 0:
            logging.info(f"Found {clash_count} clash(es) for ligand {ligand_id}")
        else:
            logging.debug(f"No clashes found for ligand {ligand_id}")
        
        return clash_count
        
    except (IOError, OSError) as e:
        logging.error(f"Error reading clash.out file {clash_file}: {e}")
        return 0
    except Exception as e:
        logging.error(f"Unexpected error parsing clash.out file: {e}")
        logging.debug(f"Full traceback: {traceback.format_exc()}")
        return 0


def validate_structure_file(file_path: str, file_type: str) -> bool:
    """
    Validate that a structure file exists and is safe to use.
    
    Args:
        file_path: Path to the file to validate
        file_type: Type of file (e.g., 'pdb', 'mtz', 'map', 'cif')
        
    Returns:
        True if file is valid, False otherwise
    """
    if not file_path or file_path == 'NA':
        return False
    
    # Check file exists
    if not os.path.exists(file_path):
        return False
    
    # Check it's a file (not a directory)
    if not os.path.isfile(file_path):
        return False
    
    # Check file is readable
    if not os.access(file_path, os.R_OK):
        return False
    
    return True


def _determine_xce_export(composite_score: Any, r_factor: Any, rfree: Any) -> str:
    """
    Determine if a structure should be exported to XCE.
    
    Criteria:
    - Composite Quality Score >= 75
    - R <= 0.3 (30%)
    - Rfree <= 0.3 (30%)
    
    Args:
        composite_score: Composite Quality Score
        r_factor: R-factor
        rfree: Rfree value
        
    Returns:
        'True' or 'False' string
    """
    # Check composite score
    if composite_score in ['NA', None]:
        return 'False'
    
    try:
        score = float(composite_score)
        if score < 75:
            return 'False'
    except (ValueError, TypeError):
        return 'False'
    
    # Check R-factor
    if r_factor not in ['NA', None]:
        try:
            r_val = float(r_factor)
            if r_val > THRESHOLDS.r_max:
                return 'False'
        except (ValueError, TypeError):
            pass  # If can't parse, don't reject based on this
    
    # Check Rfree
    if rfree not in ['NA', None]:
        try:
            rfree_val = float(rfree)
            if rfree_val > THRESHOLDS.rfree_max:
                return 'False'
        except (ValueError, TypeError):
            pass  # If can't parse, don't reject based on this
    
    return 'True'


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
    high_resolution: str = 'NA',
    low_resolution: str = 'NA',
    ligand_rsr: Optional[float] = None,
    heavy_atoms: int = 0,
    ligand_clashes: int = 0
) -> Dict[str, Any]:
    """Build a result dictionary for a dataset."""
    # Construct ligand SVG diagram path from Input_dir and CompoundCode
    input_dir = info.get('Input_dir', 'NA')
    compound_code = info.get('CompoundCode', 'NA')
    ligand_png = get_ligand_png_path(input_dir, compound_code)
    
    # Calculate quality scores
    ligand_cc = safe_dict_get(ligand_stats, 'ligandcc')
    mogul_z_angle = safe_dict_get(ligand_stats, 'mogulzangl')
    mogul_z_bond = safe_dict_get(ligand_stats, 'mogulzbond')
    ligand_mean_b = safe_dict_get(ligand_stats, 'ligandbavg')
    global_mean_b = safe_dict_get(postrefinement_stats, 'MeanB')
    
    quality_scores = calculate_quality_scores(
        ligand_cc, 
        ligand_rsr, 
        mogul_z_bond, 
        mogul_z_angle,
        ligand_mean_b,
        global_mean_b,
        ligand_clashes
    )
    
    # Calculate clashes per atom
    if heavy_atoms > 0 and ligand_clashes >= 0:
        clashes_per_atom = ligand_clashes / heavy_atoms
    else:
        clashes_per_atom = 'NA'
    
    return {
        'Crystal Name': dataset,
        'Compound Code': compound_code,
        'Ligand Structure': ligand_png,
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
        'Ligand CC': ligand_cc,
        'Ligand RSR': ligand_rsr if ligand_rsr is not None else 'NA',
        'Ligand Fit Quality Score': quality_scores['Ligand Fit Quality Score'],
        'Ligand Fit PCA Percentile (Batch)': 'NA',  # Will be calculated after all results collected
        'Ligand Fit PCA Percentile (PDB)': 'NA',  # Will be calculated after all results collected
        'Mogul Z Angle': mogul_z_angle,
        'Mogul Z Bond': mogul_z_bond,
        'Ligand Geometry Quality Score': quality_scores['Ligand Geometry Quality Score'],
        'Ligand Geometry PCA Percentile (Batch)': 'NA',  # Will be calculated after all results collected
        'Ligand Geometry PCA Percentile (PDB)': 'NA',  # Will be calculated after all results collected
        'PCA Percentile Plot (Batch)': 'NA',  # Will be generated after batch percentiles calculated
        'PCA Percentile Plot (PDB)': 'NA',  # Will be generated after PDB percentiles calculated
        'Spider Plot': 'NA',  # Will be generated after all metrics collected
        'Composite Quality Score': quality_scores['Composite Quality Score'],
        'R': safe_round(safe_dict_get(postrefinement_stats, 'R')),
        'Rfree': safe_round(safe_dict_get(postrefinement_stats, 'Rfree')),
        'Ligand Mean B Factor': safe_dict_get(ligand_stats, 'ligandbavg'),
        'Ligand Occupancy': safe_dict_get(ligand_stats, 'ligandomin'),
        'Ligand Atoms': heavy_atoms if heavy_atoms > 0 else 'NA',
        'Ligand Clashes': ligand_clashes if ligand_clashes >= 0 else 'NA',
        'Ligand Clashes per Atom': clashes_per_atom if clashes_per_atom != 'NA' else 'NA',
        'C-beta Deviations': safe_dict_get(molprobity_stats, 'cbetadeviations'),
        'Rama Outlier Percent': safe_dict_get(molprobity_stats, 'ramaoutlierpercent'),
        'Rama Favored Percent': safe_dict_get(molprobity_stats, 'ramafavoredpercent'),
        'Poor Rotamer Percent': safe_dict_get(molprobity_stats, 'poorrotamerspercent'),
        'Clash Score': safe_dict_get(molprobity_stats, 'clashscore'),
        'Mol Probity Score': safe_dict_get(molprobity_stats, 'molprobityscore'),
        'RMS Bonds': safe_dict_get(molprobity_stats, 'rmsbonds'),
        'RMS Angles': safe_dict_get(molprobity_stats, 'rmsangles'),
        'Mean B Factor': safe_dict_get(postrefinement_stats, 'MeanB'),
        'High Resolution (Å)': high_resolution,
        'Low Resolution (Å)': low_resolution,
        'Comments': '',
        'Export to XCE': _determine_xce_export(quality_scores['Composite Quality Score'], 
                                               safe_dict_get(postrefinement_stats, 'R'),
                                               safe_dict_get(postrefinement_stats, 'Rfree'))
    }


def collect_results_from_json(json_data: dict, config: Optional['Config'] = None) -> List[Dict[str, Any]]:
    """Collect results from the JSON data for all datasets."""
    if config is None:
        config = Config()
    
    results = []
    total_datasets = len(json_data)
    logging.info(f"Processing {total_datasets} datasets")
    
    # Progress bar setup
    if HAS_TQDM:
        # Show progress bar
        print(f"Processing {total_datasets} datasets...")
        iterator = tqdm(
            json_data.items(),
            total=total_datasets,
            desc="Processing datasets",
            unit="dataset",
            ncols=80,
            bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]'
        )
    else:
        # No tqdm: show periodic updates
        print(f"Processing {total_datasets} datasets...")
        iterator = json_data.items()
    
    for i, (dataset, info) in enumerate(iterator, 1):
        logging.info(f"Processing dataset {i}/{total_datasets}: {dataset}")
        
        # Show progress updates without tqdm
        if not HAS_TQDM:
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
        
        # Extract low resolution from dataprocessing.inputdata.reslo
        low_resolution = safe_get(summary, ['dataprocessing', 'inputdata', 'reslo'], 'NA')
        if low_resolution != 'NA':
            try:
                low_resolution = f"{float(low_resolution):.2f}"
            except (ValueError, TypeError):
                low_resolution = 'NA'
        
        # Extract validation statistics - handle multiple ligands
        ligands_array = safe_get(summary, ['ligandfitting', 'ligands'], [])
        if not isinstance(ligands_array, list):
            ligands_array = [ligands_array] if ligands_array else []
        
        # Check if ligands were found
        if not ligands_array:
            logging.warning(f"No ligands found in summary for {dataset}")
            ligands_array = [{}]
        
        # Get the first ligand entry (usually there's only one at this level)
        first_ligand = ligands_array[0] if ligands_array else {}
        
        # Extract ligandstatistics array - this is where multiple ligands actually are
        ligandstatistics_array = safe_get(first_ligand, ['validationstatistics', 'ligandstatistics'], [])
        if not isinstance(ligandstatistics_array, list):
            ligandstatistics_array = [ligandstatistics_array] if ligandstatistics_array else []
        
        # If no ligand statistics found, create a dummy entry
        if not ligandstatistics_array:
            logging.debug(f"No ligand statistics found for {dataset}, using empty entry")
            ligandstatistics_array = [{}]
        
        num_ligands = len(ligandstatistics_array)
        if num_ligands > 1:
            logging.info(f"Dataset {dataset} contains {num_ligands} ligands - will create separate entries")
        
        # Extract key directory and compound information (shared across all ligands)
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
        
        # Get all file paths using helper classes (shared files)
        pdb_file = pipedream_paths.postrefine_file(REFINE_PDB_FILENAME)
        output_cif_file = pipedream_paths.rhofit_file(BEST_CIF_FILENAME)
        mtz_file = pipedream_paths.postrefine_file(REFINE_MTZ_FILENAME)
        # Note: electron_density_gif will be retrieved per-ligand below

        # Get input files
        input_files = input_paths.all_files()
        input_smi_path = input_files['smiles']
        input_pdb_path = input_files['pdb']
        input_cif_path = input_files['cif']
        
        # Fallback to JSON parsing if needed (for first ligand)
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
        
        # Look for existing edstats output file (don't regenerate)
        map_2fofc_file = 'NA'
        map_fofc_file = 'NA'
        edstats_output = 'NA'
        
        # Check for existing edstats.out file in postrefine directory
        # This aligns with pipedream_xchem.py and pipedream_post_process.py which save edstats.out in postrefine-{compound_code}/
        edstats_path = pipedream_paths.postrefine_file('edstats.out')
        if edstats_path != 'NA':
            edstats_output = edstats_path
            logging.debug(f"Using existing edstats output: {edstats_output}")
        else:
            logging.debug(f"No edstats output found for dataset: {dataset}")
        
        # Use best.pdb for chirality analysis
        chirality_pdb_file = pipedream_paths.rhofit_file(BEST_PDB_FILENAME)
        
        # Process each ligand separately
        for ligand_idx, ligand_stats in enumerate(ligandstatistics_array):
            # ligand_stats now contains the individual ligand statistics
            # Get molprobity and postrefinement stats from the parent ligand entry
            molprobity_stats = safe_get(first_ligand, ['validationstatistics', 'molprobity'])
            postrefinement_stats = safe_get(first_ligand, ['postrefinement', 1])
            
            # Get ligand ID if available
            ligand_id = safe_dict_get(ligand_stats, 'ligandid', 'NA')
            
            # Get ligand-specific electron density GIF
            electron_density_gif = pipedream_paths.find_electron_density_gif(ligand_id)
            
            # Parse ligand-specific RSR from edstats if available
            ligand_rsr = None
            if edstats_output != 'NA' and ligand_id != 'NA':
                # Extract just the residue type from ligand_id (e.g., "LIG A4000" -> "LIG")
                # ligand_id format is typically: "RESTYPE CHAIN RESNUM" or just "RESTYPE"
                residue_type = ligand_id.split()[0] if ' ' in ligand_id else ligand_id
                ligand_rsr = parse_edstats_for_ligand_rm(edstats_output, residue_type)
                if ligand_rsr is not None and ligand_rsr != 'NA':
                    logging.debug(f"Extracted Ligand RSR (Rm) for {dataset} ligand {ligand_id}: {ligand_rsr}")
            
            # Parse ligand-specific clashes from clash.out if available
            ligand_clashes = 0
            if pipedream_dir != 'NA' and compound_code != 'NA' and ligand_id != 'NA':
                # Construct path to clash.out: PipedreamDirectory/report-<compound_code>/molprobe/clash.out
                clash_file = os.path.join(pipedream_dir, f"report-{compound_code}", "molprobe", "clash.out")
                if os.path.exists(clash_file):
                    ligand_clashes = parse_clash_out_for_ligand(clash_file, ligand_id)
                    logging.debug(f"Extracted {ligand_clashes} clash(es) for {dataset} ligand {ligand_id}")
                else:
                    logging.debug(f"clash.out not found at: {clash_file}")
            
            # Perform chirality analysis for each ligand
            # Always prioritize using CIF input file when available
            # Unified chirality detection - prefers CIF, falls back to SMILES template
            heavy_atoms = 0  # Initialize atom count
            
            if input_smi_path != 'NA' and chirality_pdb_file != 'NA':
                logging.debug(f"Performing chirality analysis for {dataset} ligand {ligand_id}")
                
                # Unified function handles both CIF (preferred) and SMILES fallback approaches
                chirality_flip, updated_smiles, heavy_atoms = compare_chiral_centers_cif_pdb(
                    input_cif_path, input_pdb_path, output_cif_file, chirality_pdb_file, input_smi_path, dataset)
                
                # Determine method used from result prefix
                method = "CIF+PDB" if "[CIF]" in chirality_flip else "SMILES+PDB"
                
                # Save SMILES file
                smiles_file = _save_smiles_file(pipedream_dir, compound_code, 
                                                chirality_flip, updated_smiles, method, dataset)
            
            else:
                missing_info = []
                if input_smi_path == 'NA':
                    missing_info.append("SMILES file")
                if chirality_pdb_file == 'NA':
                    missing_info.append("output PDB file")
                chirality_flip = f"Missing: {', '.join(missing_info)}"
                smiles_file = "NA"
                logging.warning(f"Skipping chirality analysis for {dataset} ligand {ligand_id}: {chirality_flip}")
            
            dataset_name = dataset
            if num_ligands > 1:
                logging.debug(f"Creating entry for {dataset_name} with ligand {ligand_id}")
            
            # Build result for this ligand
            result = build_result(
                dataset=dataset_name,
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
                high_resolution=high_resolution,
                low_resolution=low_resolution,
                ligand_rsr=ligand_rsr,
                heavy_atoms=heavy_atoms,
                ligand_clashes=ligand_clashes
            )
            results.append(result)
    
    return results


def get_cell_style_and_value(col: str, cell_val: Any, row: pd.Series, df: pd.DataFrame = None) -> tuple[str, str]:
    """Get the appropriate styling and display value for a table cell."""
    style = ''
    
    # Check for NaN/nan/None and standardize to "NA"
    if pd.isna(cell_val) or cell_val is None or str(cell_val).lower() in ['nan', 'none', '']:
        display_val = 'NA'
    else:
        display_val = str(cell_val)
        
        # Try to format numbers to 2 decimal places if possible (except for integer columns)
        if col not in ['Ligand Atoms', 'Ligand Clashes', 'C-beta Deviations']:  # Integer columns - no decimals
            try:
                num_val = float(cell_val)
                display_val = f"{num_val:.2f}"
            except (ValueError, TypeError):
                pass
    
    # Special styling for Crystal Name - purple border if dataset has multiple ligands
    if col == 'Crystal Name' and df is not None:
        crystal_name = str(cell_val)
        # Count how many entries have this crystal name (same dataset)
        same_dataset_count = len(df[df['Crystal Name'] == crystal_name])
        if same_dataset_count > 1:
            style = 'font-weight: bold; border: 3px solid #800080;'  # Purple border for multi-ligand
    
    # Special styling for Compound Code - green/yellow/red based on quality ranking
    if col == 'Compound Code' and df is not None:
        compound_code = str(cell_val)
        if compound_code != 'NA':
            # Get all rows with this compound code
            compound_rows = df[df['Compound Code'] == compound_code]
            num_occurrences = len(compound_rows)
            
            if num_occurrences > 1:
                # Multiple entries for this compound - rank by Composite Quality Score
                try:
                    # Get composite quality scores for this compound
                    quality_scores = pd.to_numeric(compound_rows['Composite Quality Score'], errors='coerce')
                    max_score = quality_scores.max()
                    min_score = quality_scores.min()
                    current_score = pd.to_numeric(row['Composite Quality Score'], errors='coerce')
                    
                    if pd.notna(current_score) and pd.notna(max_score) and pd.notna(min_score):
                        # If this row has the highest score, highlight it green
                        if current_score == max_score:
                            base_color = '#B8E6B8'  # Green for best quality
                        # If this row has the lowest score, highlight it red
                        elif current_score == min_score:
                            base_color = '#FFB3B3'  # Red for worst quality
                        # Everything in between is yellow
                        else:
                            base_color = '#FFF4C4'  # Yellow for middle quality
                        
                        # Add purple border for duplicates
                        style = f'background-color: {base_color}; font-weight: bold; border: 3px solid #800080;'
                except (ValueError, TypeError, KeyError):
                    pass
            else:
                # Single occurrence - color it green
                style = 'background-color: #B8E6B8; font-weight: bold;'
    
    # Apply numeric highlighting rules (only if no style already set)
    if not style:
        try:
            val = float(cell_val)
            if col == 'Ligand CC':
                # Color-coded Ligand CC thresholds with consistent pastel colors
                if val >= 0.95:  # Very good fit
                    style = 'background-color: #B8E6B8; font-weight: bold;'  # Pastel green
                elif val >= 0.90:  # Good fit
                    style = 'background-color: #D4F1D4; font-weight: bold;'  # Lighter pastel green
                elif val >= 0.80:  # OK fit
                    style = 'background-color: #FFF4C4; font-weight: bold;'  # Pastel yellow
                else:  # Poor fit (< 0.80)
                    style = 'background-color: #FFB3B3; font-weight: bold;'  # Pastel red
            elif col == 'Ligand RSR':
                # Color-coded Ligand RSR thresholds (lower is better)
                if val < 0.2:  # Very good
                    style = 'background-color: #B8E6B8; font-weight: bold;'  # Pastel green
                elif val < 0.3:  # Good
                    style = 'background-color: #D4F1D4; font-weight: bold;'  # Lighter pastel green
                elif val < 0.4:  # OK
                    style = 'background-color: #FFF4C4; font-weight: bold;'  # Pastel yellow
                else:  # Poor (>= 0.4)
                    style = 'background-color: #FFB3B3; font-weight: bold;'  # Pastel red
            elif col in ['Ligand Fit Quality Score', 'Ligand Geometry Quality Score', 'Composite Quality Score',
                        'Ligand Fit PCA Percentile (Batch)', 'Ligand Geometry PCA Percentile (Batch)',
                        'Ligand Fit PCA Percentile (PDB)', 'Ligand Geometry PCA Percentile (PDB)']:
                # Quality scores are 0-100, higher is better
                if val >= 75:  # Excellent
                    style = 'background-color: #B8E6B8; font-weight: bold;'  # Pastel green
                elif val >= 50:  # Good/OK
                    style = 'background-color: #FFF4C4; font-weight: bold;'  # Pastel yellow
                else:  # Poor (< 50)
                    style = 'background-color: #FFB3B3; font-weight: bold;'  # Pastel red
            elif col in ['Mogul Z Angle', 'Mogul Z Bond']:
                # Color-coded Mogul Z score thresholds with consistent pastel colors
                if val > THRESHOLDS.mogul_z_poor:  # Z > 4: bad
                    style = 'background-color: #D4B8E6; font-weight: bold;'  # Pastel purple
                elif val > THRESHOLDS.mogul_z_ok:  # 2.5 < Z < 4: poor
                    style = 'background-color: #E6D4F1; font-weight: bold;'  # Lighter pastel purple
                elif val > THRESHOLDS.mogul_z_good:  # 1.5 < Z < 2.5: ok
                    style = 'background-color: #FFF4C4; font-weight: bold;'  # Pastel yellow
                else:  # Z < 1.5: good
                    style = 'background-color: #B8E6B8; font-weight: bold;'  # Pastel green
            elif col == 'Rfree' and val > THRESHOLDS.rfree_max:
                style = 'background-color: #FFB3B3; font-weight: bold;'  # Pastel red
            elif col == 'R' and val > THRESHOLDS.r_max:
                style = 'background-color: #FFB3B3; font-weight: bold;'  # Pastel red
            elif col == 'Ligand Mean B Factor':
                try:
                    mean_b = float(row['Mean B Factor'])
                    if mean_b != 0 and val / mean_b >= THRESHOLDS.b_factor_ratio:
                        style = 'background-color: #FFB3B3; font-weight: bold;'  # Pastel red
                except (ValueError, TypeError, ZeroDivisionError):
                    pass
            elif col == 'Rama Favored Percent' and val >= 98.0:
                style = 'background-color: #B8E6B8; font-weight: bold;'  # Pastel green
            elif col == 'Rama Outlier Percent' and val > THRESHOLDS.rama_outlier_max:
                style = 'background-color: #FFB3B3; font-weight: bold;'  # Pastel red
            elif col == 'Poor Rotamer Percent' and val > THRESHOLDS.poor_rotamer_max:
                style = 'background-color: #FFB3B3; font-weight: bold;'  # Pastel red
            elif col == 'Clash Score' and val > THRESHOLDS.clash_score_max:
                style = 'background-color: #FFB3B3; font-weight: bold;'  # Pastel red
            elif col == 'Mol Probity Score' and val > THRESHOLDS.molprobity_score_max:
                style = 'background-color: #FFB3B3; font-weight: bold;'  # Pastel red
            elif col == 'RMS Angles' and val > THRESHOLDS.rms_angles_max:
                style = 'background-color: #FFB3B3; font-weight: bold;'  # Pastel red
            elif col == 'RMS Bonds' and val > THRESHOLDS.rms_bonds_max:
                style = 'background-color: #FFB3B3; font-weight: bold;'  # Pastel red
        except (ValueError, TypeError):
            pass
    
    return style, display_val


def get_cell_html(col: str, cell_val: Any, row: pd.Series, df: pd.DataFrame = None) -> str:
    """Generate HTML for a specific table cell based on column type."""
    if col == 'Export to XCE':
        return f'<td><input type="checkbox" {"checked" if cell_val == "True" else ""}></td>'
    elif col == 'Comments':
        return f'<td><input type="text" value="{cell_val}"></td>'
    elif col == 'Ligand Structure':
        return f'<td><img src="file://{cell_val}" alt="Ligand Image" width="150"></td>' if cell_val and cell_val != 'NA' else '<td></td>'
    elif col == 'Ligand Density':
        return f'<td><img src="file://{cell_val}" alt="Ligand Density" width="180" style="max-width:180px;max-height:120px;"></td>' if cell_val and cell_val != 'NA' else '<td></td>'
    elif col == 'PCA Percentile Plot (Batch)':
        # Check if cell_val is a valid path (not the column name itself)
        if cell_val and cell_val != 'NA' and str(cell_val).endswith('.html'):
            # Return button styled like the other plot buttons
            return f'''<td><a href="#" class="btn btn-sm btn-outline-primary file-link" 
                data-file="file://{os.path.abspath(cell_val)}" data-ext=".html" 
                style="padding:5px 10px;text-decoration:none;">View</a></td>'''
        else:
            return '<td></td>'
    elif col == 'PCA Percentile Plot (PDB)':
        # Check if cell_val is a valid path (not the column name itself)
        if cell_val and cell_val != 'NA' and str(cell_val).endswith('.html'):
            # Return button styled like the other plot buttons
            return f'''<td><a href="#" class="btn btn-sm btn-outline-primary file-link" 
                data-file="file://{os.path.abspath(cell_val)}" data-ext=".html" 
                style="padding:5px 10px;text-decoration:none;">View</a></td>'''
        else:
            return '<td></td>'
    elif col == 'Spider Plot':
        # Check if cell_val is a valid path (not the column name itself)
        if cell_val and cell_val != 'NA' and str(cell_val).endswith('.html'):
            # Return button styled like the other plot buttons
            return f'''<td><a href="#" class="btn btn-sm btn-outline-primary file-link" 
                data-file="file://{os.path.abspath(cell_val)}" data-ext=".html" 
                style="padding:5px 10px;text-decoration:none;">View</a></td>'''
        else:
            return '<td></td>'
    elif col == 'Chiral inversion':
        # Check if stereochemistry matches (starts with 'Stereochemistry matches input')
        is_match = cell_val and (cell_val.startswith('Stereochemistry matches input') or cell_val == 'Not checked')
        return f'<td>{cell_val}</td>' if is_match else f'<td style="color:red;font-weight:bold">{cell_val}</td>'
    elif col in [
        'Pipedream Directory', 'Buster Report HTML', 'Ligand Report HTML', 'Pipedream Summary',
        'PDB File', 'MTZ File', '2Fo-Fc Map File', 'Fo-Fc Map File', 'Output SMILES']:
        return get_file_link_html(cell_val)
    else:
        # Regular cell with potential highlighting
        style, display_val = get_cell_style_and_value(col, cell_val, row, df)
        return f'<td style="{style}">{display_val}</td>'


def save_results_to_html(results: List[Dict[str, Any]], output_file: str, open_browser: bool = True,
                         quality_plot: str = 'NA', ligand_pca_plot: str = 'NA', 
                         geometry_pca_plot: str = 'NA', pdb_quality_plot: str = 'NA',
                         pca_ref: Optional[PCAReference] = None) -> None:
    """Save results to an interactive HTML file.
    
    Note: PCA scores should already be calculated before calling this function.
    """
    logging.info(f"Generating HTML output with {len(results)} results")
    print(f"Generating HTML report...")


    # Implementation

    df = pd.DataFrame(results)
    
    # Ensure all columns from COLUMN_ORDER exist in the DataFrame
    for col in COLUMN_ORDER:
        if col not in df.columns:
            df[col] = 'NA'
            logging.debug(f"Added missing column to DataFrame: {col}")

    # Ensure 'Composite Quality Score' is numeric for sorting
    df['Composite Quality Score'] = pd.to_numeric(df['Composite Quality Score'], errors='coerce')
    df = df.sort_values(by='Composite Quality Score', ascending=False)
    logging.debug(f"Sorted results by Composite Quality Score")

    column_order = COLUMN_ORDER

    # Convert DataFrame to list of dicts for template
    results_data = df.to_dict('records')
    
    quality_score_index = column_order.index('Composite Quality Score') if 'Composite Quality Score' in column_order else 0

    # Helper functions for template
    def get_cell_content(col, cell_val, row_dict):
        """Get cell content for template - returns dict with type and data."""
        # Convert dict back to Series for existing helper functions
        row_series = pd.Series(row_dict)
        
        if col == 'Export to XCE':
            return {'type': 'checkbox', 'checked': cell_val == "True"}
        elif col == 'Comments':
            return {'type': 'input', 'value': cell_val}
        elif col == 'Ligand Structure':
            return {'type': 'image', 'src': cell_val, 'alt': 'Ligand Image', 'width': '150'} if cell_val and cell_val != 'NA' else {'type': 'empty'}
        elif col == 'Ligand Density':
            return {'type': 'image', 'src': cell_val, 'alt': 'Ligand Density', 'width': '180', 'style': 'max-width:180px;max-height:120px;'} if cell_val and cell_val != 'NA' else {'type': 'empty'}
        elif col in ['PCA Percentile Plot (Batch)', 'PCA Percentile Plot (PDB)', 'Spider Plot']:
            if cell_val and cell_val != 'NA' and str(cell_val).endswith('.html'):
                return {'type': 'file_button', 'path': os.path.abspath(cell_val), 'ext': '.html'}
            return {'type': 'empty'}
        elif col == 'Chiral inversion':
            # Check if stereochemistry matches (starts with 'Stereochemistry matches input')
            is_match = cell_val and (cell_val.startswith('Stereochemistry matches input') or cell_val == 'Not checked')
            return {'type': 'text', 'value': cell_val, 'style': 'color:red;font-weight:bold' if not is_match else ''}
        elif col in ['Pipedream Directory', 'Buster Report HTML', 'Ligand Report HTML', 'Pipedream Summary',
                     'PDB File', 'MTZ File', '2Fo-Fc Map File', 'Fo-Fc Map File', 'Output SMILES']:
            if not cell_val or cell_val == 'NA':
                return {'type': 'empty'}
            file_ext = os.path.splitext(cell_val)[1].lower()
            return {'type': 'file_button', 'path': cell_val, 'ext': file_ext}
        else:
            style, display_val = get_cell_style_and_value(col, cell_val, row_series, df)
            return {'type': 'text', 'value': display_val, 'style': style}

    # Prepare template context
    template_context = {
        'original_json_filename': os.path.basename(output_file).replace('.html', '.json'),
        'quality_score_index': quality_score_index,
        'column_order': column_order,
        'column_tooltips': COLUMN_TOOLTIPS,
        'results': results_data,
        'get_cell_content': get_cell_content,
        'os': os,  # For os.path functions
        'quality_plot_exists': quality_plot != 'NA' and os.path.exists(quality_plot),
        'quality_plot_path': os.path.abspath(quality_plot) if quality_plot != 'NA' and os.path.exists(quality_plot) else '',
        'pdb_quality_plot_exists': pdb_quality_plot != 'NA' and os.path.exists(pdb_quality_plot),
        'pdb_quality_plot_path': os.path.abspath(pdb_quality_plot) if pdb_quality_plot != 'NA' and os.path.exists(pdb_quality_plot) else '',
        'ligand_pca_plot_exists': ligand_pca_plot != 'NA' and os.path.exists(ligand_pca_plot),
        'ligand_pca_plot_path': os.path.abspath(ligand_pca_plot) if ligand_pca_plot != 'NA' and os.path.exists(ligand_pca_plot) else '',
        'geometry_pca_plot_exists': geometry_pca_plot != 'NA' and os.path.exists(geometry_pca_plot),
        'geometry_pca_plot_path': os.path.abspath(geometry_pca_plot) if geometry_pca_plot != 'NA' and os.path.exists(geometry_pca_plot) else '',
    }

    # Load and render template
    template_dir = os.path.join(os.path.dirname(__file__), 'templates')
    env = Environment(
        loader=FileSystemLoader(template_dir),
        autoescape=select_autoescape(['html', 'xml'])
    )
    template = env.get_template('results_table.html')
    html_content = template.render(**template_context)

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
        dest='input',
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
        '--no-plots',
        action='store_true',
        help='Don\'t automatically open plot PNG files'
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
    if not os.path.isfile(args.input):
        parser.error(f"Input JSON file not found: {args.input}")
    
    # Validate input file extension
    if not args.input.lower().endswith('.json'):
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
        output_dir = os.path.dirname(args.input)
    
    # Determine output file names
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    if args.output_name:
        base_name = args.output_name
    else:
        base_name = f'Pipedream_results_{timestamp}'
    
    output_file_json = os.path.join(output_dir, f'{base_name}.json')
    output_file_html = os.path.join(output_dir, f'{base_name}.html')
    
    log_file = os.path.join(output_dir, f'{base_name}.log')
    
    # Configure logging
    logger = logging.getLogger()
    logger.setLevel(getattr(logging, args.log_level))
    
    # File handler - captures everything based on log level
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(getattr(logging, args.log_level))
    file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(file_formatter)
    
    # Console handler - respects verbose flag
    console_handler = logging.StreamHandler()
    if args.verbose:
        console_handler.setLevel(getattr(logging, args.log_level))
    else:
        console_handler.setLevel(logging.WARNING)  # Only warnings and errors when not verbose
    console_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console_handler.setFormatter(console_formatter)
    
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    
    # Create config without file path
    config = Config()
    
    # Try to load PDB PCA reference data
    pca_ref = PCAReference()
    for ref_path in config.pca_reference_paths:
        if os.path.exists(str(ref_path)):
            if pca_ref.load_reference(str(ref_path)):
                print(f"Loaded PDB archive reference data: {ref_path}")
                break
    
    if not pca_ref.loaded:
        print("Warning: PDB archive reference data not found - PDB percentile scores will be skipped")
        print("Place 'PCA_PDB_references.csv' in script directory or Data/ subdirectory to enable")
    
    logging.info("Starting Pipedream results collation")
    print("Starting Pipedream results collation...")
    
    # Show note if tqdm is not available
    if not HAS_TQDM:
        print("Note: Install tqdm for progress bar support (pip install tqdm)")
    
    logging.info(f"Reading JSON file: {args.input}")
    json_data = read_json(args.input)
    
    logging.info("Validating JSON structure")
    validate_json_structure(json_data)
    
    logging.info("Collecting results from datasets")
    results = collect_results_from_json(json_data, config)
    
    # Calculate PCA scores BEFORE saving JSON or generating HTML
    logging.info("Calculating PCA percentile scores (batch)")
    print("Calculating PCA percentile scores...")
    results = calculate_pca_scores_batch(results, pca_ref)
    
    # Generate individual batch percentile plots for each dataset
    logging.info("Generating individual batch percentile plots...")
    print("Generating individual batch percentile plots...")
    for result in results:
        plot_path = generate_individual_batch_percentile_plot(result, output_dir)
        result['PCA Percentile Plot (Batch)'] = plot_path
        logging.info(f"Set PCA Percentile Plot (Batch) for {result.get('Crystal Name')}: {plot_path}")
    
    logging.info("Calculating PCA percentile scores (PDB)")
    results = calculate_pca_scores_pdb(results, pca_ref)
    
    # Generate individual PDB percentile plots for each dataset
    logging.info("Generating individual PDB percentile plots...")
    print("Generating individual PDB percentile plots...")
    for result in results:
        plot_path = generate_individual_pdb_percentile_plot(result, output_dir)
        result['PCA Percentile Plot (PDB)'] = plot_path
        logging.info(f"Set PCA Percentile Plot (PDB) for {result.get('Crystal Name')}: {plot_path}")
    
    # Generate spider plots for each dataset
    logging.info("Generating spider plots...")
    print("Generating spider plots...")
    for result in results:
        spider_path = generate_spider_plot(result, output_dir)
        result['Spider Plot'] = spider_path
        logging.info(f"Set Spider Plot for {result.get('Crystal Name')}: {spider_path}")
    
    # Generate 2D PDB percentile plot (after PDB scores are calculated)
    logging.info("Generating 2D PDB percentile plot")
    print("Generating 2D PDB percentile plot...")
    pdb_quality_plot = generate_pdb_percentile_quality_plot(results, output_dir, base_name)
    
    # Generate PCA plots
    logging.info("Generating PCA plots")
    print("Generating PCA plots...")
    ligand_pca_plot, geometry_pca_plot = generate_pca_plots(results, output_dir, base_name)
    
    # Generate 2D ligand quality plot (batch comparison)
    logging.info("Generating 2D ligand quality plot (batch)")
    print("Generating 2D ligand quality plot (batch)...")
    quality_plot = generate_ligand_quality_plot(results, output_dir, base_name)
    
    # Generate outputs based on format selection
    if args.format in ['json', 'both']:
        logging.info(f"Saving results to JSON: {output_file_json}")
        print(f"Saving JSON results...")
        save_results_to_json(results, output_file_json)
    
    if args.format in ['html', 'both']:
        logging.info(f"Generating HTML report: {output_file_html}")
        save_results_to_html(results, output_file_html, open_browser=not args.no_browser,
                           quality_plot=quality_plot, ligand_pca_plot=ligand_pca_plot,
                           geometry_pca_plot=geometry_pca_plot, pdb_quality_plot=pdb_quality_plot,
                           pca_ref=pca_ref)
    
    logging.info("Result collection completed.")
    print("Pipedream results collation completed!")
    
    # Print summary of what was generated
    print(f"\nOutput summary:")
    if args.format in ['json', 'both']:
        print(f"  JSON: {output_file_json}")
    if args.format in ['html', 'both']:
        print(f"  HTML: {output_file_html}")
    if ligand_pca_plot != 'NA':
        print(f"  PCA (Ligand): {ligand_pca_plot}")
    if geometry_pca_plot != 'NA':
        print(f"  PCA (Geometry): {geometry_pca_plot}")
    print(f"  Log:  {log_file}")


if __name__ == "__main__":
    main()
