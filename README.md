# Pipedream for XChem

Last edited time: July 31, 2025 1:54 PM

Last edited by: Daren Fearon

These python scripts automate the submission of Pipedream refinement jobs as an array on the Diamond cluster, collation of Pipedream outputs and re-integration of data back into the XChem environment. They support dataset selection from a soakDB.sqlite database file and prepare jobs based on user-defined parameters in the YAML parameters file.

This guide explains how to set up your environment, run the Pipedream for XChem pipeline, and troubleshoot common issues using the provided Python scripts and configuration files.

---

## Environment Setup 🔧

Before running any scripts, ensure the following environment is correctly configured:

**1. Micromamba Setup**
If not already initialized:

    `micromamba init`

**2. Activate XChem Environment**

    `micromamba activate xchem`

**3. Required Python Packages**
Ensure the following packages are installed in your environment (if using the xchem env at DLS then they should be installed):
    - rdkit
    - pandas
    - pyyaml
    - sqlite3
    - paramiko
    - gemmi

Install missing packages using:

   `micromamba install -c conda-forge rdkit pandas pyyaml gemmi`

---

## Files Overview 📁 

| File                          | Purpose                                                 |
|-------------------------------|---------------------------------------------------------|
| pipedream_xchem.py            | Submits SLURM jobs for Pipedream refinement             |
| collate_pipedream_results.py  | Aggregates results and performs chirality analysis      |
| export_pipedream_to_xce.py    | Exports results to XCE format and updates the database  |
| pipedream_parameters.yaml     | Configuration file for refinement parameters            |
| datasets.csv                  | List of specific datasets to process                    |
| soakDBDataFile.sqlite         | SQLite database with dataset metadata                   |
|-------------------------------|---------------------------------------------------------|

---

## YAML Configuration ⚙️

Create a YAML file (e.g., `pipedream_parameters.yaml`) with the following structure (or copy the template file from /dls/science/groups/i04-1/software/pipedream/pipedream_parameters.yaml):

## Parameter Descriptions

```yaml

Mode: "specific_datasets" # 'pending_analysis' or 'specific_datasets' - the former will parse your database file for all datasets with RefinementOutcome "1 - Pending Analysis", the latter will use a specified list of datasets provided in the csv file specified below

Processing_directory: /dls/labxchem/data/proposal/visit/processing
Output_directory: /dls/labxchem/data/proposal/visit/processing/Pipedream/Pipedream_output_123 # Optional - defaults to Processing_directory/analysis/Pipedream/Pipedream_<timestamp> if not set
Database_path: /dls/labxchem/data/proposal/visit/processing/database/soakDBDataFile.sqlite
Dataset_csv_path: /dls/labxchem/data/proposal/visit/processing/Pipedream/datasets.csv # Only required if Mode is 'specific_datasets'

Remove_crystallisation_components: true  # Optional - removes DMS, EDO, GOL, SO4, PO4, PEG from input PDBs if true (can skip if not modelled in site of interest in MR model)
Refinement_parameters: #For more information see https://www.globalphasing.com/buster/manual/pipedream/manual/index.html#_details_of_command_line_arguments
  keepwater: true #DO NOT remove waters that are present in the input model (default is to remove them)
  mrefine: "TLSbasic" #"TLSbasic" turns on TLS refinement and autoncs. Leave blank for no TLS.
  remediate: true #Run SideAide to refit side chains
  sidechainrebuild: true #Allow SideAide to rebuild stubbed sidechains
  runpepflip: true #Run pepflip to check for and correct peptide bond flips
  rhocommands:
    - -xclusters # Produces ligand fits for the <n> best possible binding sites. Leave blank for default and fit to <NCS> best sites.
    - -nochirals # Ignore CHIRAL restraints in fitting/output. Chiral centres can then invert as needed.
```

- **Mode**:
    - `"pending_analysis"`: Selects datasets from SoakDB with `RefinementOutcome = '1 - Analysis Pending'`.
    - `"specific_datasets"`: Uses a CSV file listing specific `CrystalName` entries.
- **Processing_directory**: Root directory, most likely your XChem visit, e.g. `/dls/labxchem/data/lb32633/lb32633-6/processing`
- **Output_directory:** Optional - if not set, defaults to `Processing_directory/analysis/Pipedream/Pipedream_<timestamp>`
- **Database_path**: Path to the SQLite database containing dataset metadata, e.g. `/dls/labxchem/data/lb32633/lb32633-6/processing/database/soakDBDataFile.sqlite`
- **Dataset_csv_path**: Path to a CSV file with a `CrystalName` column (used only in `specific_datasets` mode). Example file structure:
    
    `CrystalName
    Protein-x100
    Protein-x101
    Protein-x102`
    
- **Remove_crystallisation_components**: Optional - removes DMS, EDO, GOL, SO4, PO4, PEG from input PDBs if true (can skip if not modelled in site of interest in MR model)
- **Refinement_parameters**: Options passed to the `pipedream` command. For more information see [here](https://www.globalphasing.com/buster/manual/pipedream/manual/index.html#_details_of_command_line_arguments).
    - `keepwater: true` - DO NOT remove waters that are present in the input model (default is to remove them)
    - `mrefine: "TLSbasic"` - "TLSbasic" turns on TLS refinement and autoncs. Leave blank for no TLS.
    - `remediate: true` - Run SideAide to refit side chains
    - `sidechainrebuild: true` - Allow SideAide to rebuild stubbed sidechains
    - `runpepflip: true` - Run pepflip to check for and correct peptide bond flips
    - Rhocommands:
    - `"-xclusters <n>"` - produces ligand fits for the <n> best possible binding sites. Leave blank (i.e. “-xclusters “ for default (fit to <ncs> best sites)
    - `-nochirals` - Ignore CHIRAL restraints in fitting/output. Chiral centres can then invert as needed.

---

## Usage 🚀

**Running the Pipeline**

Python scripts can be found on DLS file system at `/dls/science/groups/i04-1/software/pipedream_xchem/`

---

### Step 1: Submit Refinement Jobs
`python pipedream_xchem.py --parameters /path/to/pipedream_parameters.yaml`

This script submits SLURM jobs for Pipedream analysis with ligand restraints generated from SMILES using Grade2. It processes crystallographic datasets from an SQLite database and prepares them for automated refinement.

It will:
1. Read the YAML parameters file
2. Query your sqlite database file for datasets
3. Prepare input folders and job scripts
4. Submit jobs to the SLURM cluster
5. Output data in specified directory or `<processing_directory>/analysis/<Pipedream_yyyymmdd_hhmmss>`

📝 **Logging**

A log file named `pipedream.log` will be created in the output directory for each run, capturing job submission details and any errors encountered.

📈 **Outputs**

Details of the pipedream output file structure can be found [here](https://www.globalphasing.com/buster/manual/pipedream/manual/index.html#_location_of_pipedream_output).

---

### Step 2: Collate Results
`python collate_pipedream_results.py --input /path/to/<Pipedream_yyyymmdd_hhmmss>_output.json`

This script collates and analyzes results from completed Pipedream crystallographic refinement runs. It processes multiple datasets, performs chirality analysis by comparing input and output ligand structures, extracts validation statistics, and generates comprehensive interactive HTML and JSON reports for quality assessment.

Features:
- Chirality inversion detection comparing input SMILES/PDB with output PDB structures
- Ligand validation statistics extraction (correlation coefficients, B-factors, etc.)
- MolProbity and structural quality metrics analysis
- SMILES extraction from refined ligand structures
- Interactive HTML reports with file links and data filtering
- Map file generation from MTZ files using gemmi

<img width="1910" height="973" alt="image" src="https://github.com/user-attachments/assets/50260180-7989-45d6-8730-3ab59e0bd599" />

This HTML file contains key validation statistics which can be used to assess and select datasets for re-integration back into the XChem/XCE environment and annotate with comments. By default ligands with a correlation coefficient >0.8 are selected for export.

*Should you wish to amend the datasets for export and/or add comments, ensure you click `Export Table to JSON` and use this output file for the next step.You can also load previously edited JSON files using `Load JSON and Update Table`*

📈 **Outputs**
1. `Pipedream_results_yyyymmdd_hhmmss.html` 
2. `Pipedream_results_yyyymmdd_hhmmss.json`
3. `Pipedream_results_yyyymmdd_hhmmss.log`

Details of the pipedream output file structure can be found [here](https://www.globalphasing.com/buster/manual/pipedream/manual/index.html#_location_of_pipedream_output).

📝 **Logging**

A log file named `Pipedream_results_yyyymmdd_hhmmss.log` will be created in the output directory for each run, capturing actions and any errors encountered.

---

### Step 3: Export to XCE
`python export_pipedream_to_xce.py --input /path/to/results.json --parameters /path/to/pipedream_parameters.yaml`

This script exports completed Pipedream refinement results to XChem Explorer (XCE) format and updates the associated SQLite database. It processes validated results from the collation script, integrates them into the XCE file structure, and handles special cases like chirality inversions by updating restraint files and molecular diagrams.

Features:
- Exports Pipedream results to XCE directory structure
- Updates SQLite database with refinement statistics and quality metrics
- Applies traffic-light scoring for validation metrics
- Handles chirality inversions by copying updated restraint files (CIF/PDB)
- Generates molecular structure diagrams (PNG) from SMILES using RDKit
- Copies BUSTER output files for downstream analysis in XCE
- Creates symlinks for PDB, MTZ, and map files in XCE format

📝 **Logging**

A log file named `export_yyyymmdd_hhmmss.log` will be created in the parent directory for each run, capturing actions and any errors encountered.

---

## Troubleshooting 🛠

- Ensure all required file paths in the database are valid and accessible.
- If a dataset is missing required files, it will be skipped and logged.
- Make sure your Wilson SSH credentials and SLURM token permissions are correctly configured.

❌ Error: “Missing required file paths”
Cause: One or more datasets lack MTZ or PDB paths in the database.
Fix:
    - Check soakDBDataFile.sqlite for missing entries.
    - Use --dry-run to identify problematic datasets:
        python pipedream_xchem.py --parameters pipedream_parameters.yaml --dry-run

❌ Error: “RDKit not available – PNG generation will be disabled”
Cause: RDKit is not installed in the environment.
Fix:
    micromamba install -c conda-forge rdkit

❌ Error: “Could not find CIF file for compound”
Cause: CIF file missing in expected location.
Fix:
    - Check input_files/ or rhofit-<CompoundCode>/ directories.
    - Ensure CIF files are named correctly (<CompoundCode>.cif).

❌ Error: “SLURM job not submitted”
Cause: SSH or SLURM token misconfiguration.
Fix:
    - Ensure SSH access to wilson.diamond.ac.uk is working.
    - Check SLURM token permissions.

## Contact 📬

For questions or issues, please contact xchem@diamond.ac.uk or daren.fearon@diamond.ac.uk

---
