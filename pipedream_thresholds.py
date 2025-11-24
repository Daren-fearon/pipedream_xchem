"""
Pipedream Validation Thresholds Module

This module defines the ValidationThresholds class and a global THRESHOLDS
instance that contains all quality assessment thresholds used throughout
the Pipedream analysis pipeline.

Author: DFearon
Date: November 2025
"""


class ValidationThresholds:
    """Validation thresholds for crystallographic quality assessment."""
    ligand_cc: float = 0.8
    ligand_rsr_max: float = 0.3
    rfree_max: float = 0.3
    r_max: float = 0.3
    mogul_z_max: float = 2.0  # Kept for backward compatibility
    mogul_z_good: float = 1.5
    mogul_z_ok: float = 2.5
    mogul_z_poor: float = 4.0
    # Above 4.0 is bad (purple)
    rama_outlier_max: float = 0.5
    poor_rotamer_max: float = 2.0
    clash_score_max: float = 20.0
    molprobity_score_max: float = 2.0
    b_factor_ratio: float = 1.5
    rms_angles_max: float = 3.0
    rms_bonds_max: float = 0.03


# Global thresholds instance
THRESHOLDS = ValidationThresholds()
