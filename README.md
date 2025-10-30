# 🧬 Pipedream for XChem

**Last edited:** October 30 2025  
**Last edited by:** Daren Fearon  

---

These Python scripts automate the submission of **Pipedream refinement jobs** as an array on the Diamond cluster, the **collation of Pipedream outputs**, and **re-integration of data back into the XChem environment**.

They support dataset selection from a `soakDB.sqlite` database and prepare jobs based on user-defined parameters in a YAML configuration file.

This guide explains how to set up your environment, run the pipeline, and troubleshoot common issues using the provided scripts and configuration files.

---

## ⚙️ Environment Setup

Before running any scripts, ensure your environment is configured correctly.

### 1. Initialize Micromamba

If not already initialized:

```bash
micromamba-init
```

### 2. Activate the XChem Environment

```bash
micromamba activate xchem
```

### 3. Required Python Packages

Ensure the following packages are installed in your environment (the `xchem` environment at DLS already includes them):

- rdkit  
- pandas  
- pyyaml  
- sqlite3  
- paramiko  
- gemmi  

Install missing packages using:

```bash
micromamba install -c conda-forge rdkit pandas pyyaml gemmi
```

---

## 📁 Files Overview

| File | Purpose |
|------|----------|
| `pipedream_xchem.py` | Submits SLURM jobs for Pipedream refinement |
| `collate_pipedream_results.py` | Aggregates results and performs chirality analysis |
| `export_pipedream_to_xce.py` | Exports results to XCE format and updates the database |
| `pipedream_parameters.yaml` | Configuration file for refinement parameters |
| `datasets.csv` | List of specific datasets to process |
| `soakDBDataFile.sqlite` | SQLite database with dataset metadata |

---

## 🧾 YAML Configuration

Create a YAML file (e.g., `pipedream_parameters.yaml`) with the following structure or copy the template from  
`/dls/science/groups/i04-1/software/pipedream_xchem/pipedream_parameters.yaml`.

```yaml
Mode: "specific_datasets"
Processing_directory: /dls/labxchem/data/proposal/visit/processing
Output_directory: /dls/labxchem/data/proposal/visit/processing/Pipedream/Pipedream_output_123
Database_path: /dls/labxchem/data/proposal/visit/processing/database/soakDBDataFile.sqlite
Dataset_csv_path: /dls/labxchem/data/proposal/visit/processing/Pipedream/datasets.csv

Remove_crystallisation_components: true
Refinement_parameters:
  keepwater: true
  WaterUpdatePkmaps: true
  TLS: "TLSbasic"
  remediate: true
  sidechainrebuild: true
  runpepflip: true
  rhocommands:
    - -xclusters
    - -nochirals
```

See the [Global Phasing Pipedream Manual](https://www.globalphasing.com/buster/manual/pipedream/manual/index.html#_details_of_command_line_arguments) for full parameter details.

---

## 🚀 Usage

Python scripts are available at:  
`/dls/science/groups/i04-1/software/pipedream_xchem/`

### Step 1: Submit Refinement Jobs

```bash
python pipedream_xchem.py --parameters /path/to/pipedream_parameters.yaml
```

This script:
1. Reads the YAML parameters file  
2. Queries the SQLite database  
3. Prepares input folders and job scripts  
4. Submits jobs to the SLURM cluster  
5. Outputs data to the specified directory or `<processing_directory>/analysis/Pipedream_<timestamp>`

#### 📊 Output Structure

<details>
<summary><b>🗂️ Top-Level Structure</b></summary>

```
Pipedream_YYYYMMDD_HHMMSS/
├── pipedream.log
├── datasets_metadata.csv
├── Pipedream_YYYYMMDD_HHMMSS_output.json
├── array_logs/
└── CrystalName/
```

| Path | Description |
|------|--------------|
| `pipedream.log` | Main log file for the entire pipeline execution. |
| `datasets_metadata.csv` | Metadata summary of processed datasets. |
| `Pipedream_YYYYMMDD_HHMMSS_output.json` | Consolidated output summary containing file paths and results. |
| `array_logs/` | Contains SLURM job array logs and index mappings. |
| `CrystalName/` | Directory for each processed crystal, including inputs and outputs. |

</details>

<details>
<summary><b>📜 array_logs/</b></summary>

```
array_logs/
├── slurm_array_index_map_*.csv
└── pipedream_array_*.out
```

| File | Description |
|------|--------------|
| `slurm_array_index_map_*.csv` | Maps SLURM array indices to dataset identifiers. |
| `pipedream_array_*.out` | Individual SLURM job logs for each array task. |

</details>

<details>
<summary><b>💎 CrystalName/</b></summary>

```
CrystalName/
├── input_files/
│   ├── *.mtz
│   ├── *.pdb
│   └── CompoundCode.{cif,pdb,smiles}
├── Pipedream_YYYYMMDD_HHMMSS/    ← created by pipeline
└── CrystalName_slurm_*.out
```

| Path | Description |
|------|--------------|
| `input_files/` | Contains structure and ligand input files. |
| `*.mtz`, `*.pdb` | Crystallographic data and model files. |
| `CompoundCode.{cif,pdb,smiles}` | Ligand definition files in multiple formats. |
| `Pipedream_YYYYMMDD_HHMMSS/` | Pipeline-generated subfolder containing processed outputs. |
| `CrystalName_slurm_*.out` | SLURM job log specific to this crystal. |

</details>

> 💡 **Tip:** The timestamped directory (`Pipedream_YYYYMMDD_HHMMSS`) uniquely identifies each Pipedream run.

---

### Step 2: Collate Results

```bash
python collate_pipedream_results.py --input /path/to/<Pipedream_yyyymmdd_hhmmss>_output.json
```

This script collates and analyzes results from completed Pipedream runs, performing chirality checks, validation metric extraction, and HTML/JSON report generation.

Outputs:
1. `Pipedream_results_yyyymmdd_hhmmss.html`
2. `Pipedream_results_yyyymmdd_hhmmss.json`
3. `Pipedream_results_yyyymmdd_hhmmss.log`

---

### Step 3: Export to XCE

```bash
python export_pipedream_to_xce.py --input /path/to/results.json --parameters /path/to/pipedream_parameters.yaml
```

This script exports validated Pipedream refinement results to XChem Explorer (XCE) format and updates the SQLite database.

---

## 🛠 Troubleshooting

| Error | Cause | Fix |
|-------|--------|-----|
| `Missing required file paths` | One or more datasets lack MTZ or PDB paths in the database. | Check `soakDBDataFile.sqlite` for missing entries or use `--dry-run`. |
| `Could not find CIF file for compound` | CIF file missing in expected location. | Check `input_files/` or `rhofit-<CompoundCode>/` directories. |
| `SLURM job not submitted` | SSH or SLURM token misconfiguration. | Ensure SSH access to `wilson.diamond.ac.uk` and valid SLURM permissions. |

---

## 📬 Contact

For questions or issues, please contact:  
📧 [xchem@diamond.ac.uk](mailto:xchem@diamond.ac.uk) or [daren.fearon@diamond.ac.uk](mailto:daren.fearon@diamond.ac.uk)
