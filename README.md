# Pipedream for XChem

**Last edited:** 26th November 2025  
**Last edited by:** Daren Fearon  

---

These Python scripts automate the submission of **Pipedream jobs** as SLURM arrays on the Diamond cluster, the **collation of Pipedream outputs**, and **re-integration of data back into the XChem environment**.

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
- tqdm (optional, for progress bars)

Install missing packages using:

```bash
micromamba install -c conda-forge rdkit pandas pyyaml gemmi paramiko tqdm
```
---

## 🔑 SSH Key Setup (Optional & first-time only)

The scripts connect to `wilson.diamond.ac.uk` to submit SLURM jobs. You can set up SSH key authentication before running the pipeline or the script should default to requesting a password.

### 1. Generate SSH Key Pair

If you don't already have an SSH key:

```bash
ssh-keygen -t rsa -b 4096
```

- Press **Enter** to accept the default location (`~/.ssh/id_rsa`)
- Enter a **passphrase** (minimum 5 characters - **required by Diamond IT policy**)
- Confirm the passphrase

### 2. Copy Public Key to Wilson

```bash
touch ~/.ssh/authorized_keys && chmod 600 ~/.ssh/authorized_keys && cat ~/.ssh/id_rsa.pub >> ~/.ssh/authorized_keys
```

### 3. Test Your Connection

```bash
ssh wilson.diamond.ac.uk
```

✅ If successful, you should connect without entering a password (only your SSH key passphrase if you set one).

### How the Script Uses SSH Keys

The pipeline automatically:
1. **Tries key-based authentication first** (looks for keys in `~/.ssh/`)
2. **Falls back to password** if key authentication fails
3. Uses `ssh-agent` if available to avoid repeated passphrase prompts

> ⚠️ **Security Note:** Your **private key** (`~/.ssh/id_rsa`) should **never be shared**. Only the **public key** (`~/.ssh/id_rsa.pub`) is copied to remote servers.

---

## 📁 Files Overview

| File | Purpose |
|------|----------|
| `pipedream_xchem.py` | Submits SLURM jobs for Pipedream refinement (includes automatic map generation and edstats) |
| `pipedream_post_process.py` | Standalone script for reprocessing old data (backwards compatibility) |
| `collate_pipedream_results.py` | Aggregates results and performs chirality analysis |
| `export_pipedream.py` | Exports results to XCE format and updates the database |
| `pipedream_parameters.yaml` | Configuration file for refinement parameters |
| `datasets.csv` | List of specific datasets to process |
| `soakDBDataFile.sqlite` | SQLite database with dataset metadata |

---

## 🧾 YAML Configuration

Create a YAML file (e.g., `pipedream_parameters.yaml`) with the following structure or copy the template from  
`/dls/science/groups/i04-1/software/pipedream_xchem/pipedream_parameters.yaml`.

```yaml
Mode: "specific_datasets"  # 'pending_analysis' or 'specific_datasets'
Processing_directory: /dls/labxchem/data/proposal/visit/processing
Output_directory: /dls/labxchem/data/proposal/visit/processing/Pipedream/Pipedream_output_123  # Optional
Database_path: /dls/labxchem/data/proposal/visit/processing/database/soakDBDataFile.sqlite
Dataset_csv_path: /dls/labxchem/data/proposal/visit/processing/Pipedream/datasets.csv  # Required for 'specific_datasets' mode

Cluster_partition: "cs05r"  # Optional: cs05r or cs04r (default: cs05r)
Job_priority: "low"          # Optional: normal, low, high (default: normal)

Remove_crystallisation_components: true  # Optional
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

### Parameter Descriptions

- **Mode**: `'pending_analysis'` (queries database for RefinementOutcome "1 - Pending Analysis") or `'specific_datasets'` (uses CSV list)
- **Processing_directory**: Base directory for XChem processing
- **Output_directory**: Optional - defaults to `<processing_directory>/analysis/Pipedream/Pipedream_<timestamp>` if not set
- **Database_path**: Path to SQLite database with dataset metadata
- **Dataset_csv_path**: **Required for 'specific_datasets' mode** - CSV file with `CrystalName` column
- **Cluster_partition**: Optional - SLURM partition to use (`cs05r` or `cs04r`, default: `cs05r`)
- **Job_priority**: Optional - Job priority using SLURM nice value (`normal`, `low`, `high`, default: `normal`). `low` sets nice=1000 (runs after other jobs), `high` sets nice=-100 (higher priority)
- **Remove_crystallisation_components**: Optional - removes DMS, EDO, GOL, SO4, PO4, PEG from input PDBs
- **Refinement_parameters**: Options passed to Pipedream ([see documentation](https://www.globalphasing.com/buster/manual/pipedream/manual/index.html#_details_of_command_line_arguments))

---

## 🚀 Usage

Python scripts are available at:  
`/dls/science/groups/i04-1/software/pipedream_xchem/`

### Step 1: Submit Refinement Jobs

```bash
python /dls/science/groups/i04-1/software/pipedream_xchem/pipedream_xchem.py --parameters /path/to/pipedream_parameters.yaml [OPTIONS]
```

#### Command-Line Options

**Required:**
- `--parameters`, `-p` - Path to YAML parameters file

**Output Options:**
- `--output-json PATH` - Custom path for JSON metadata output (default: auto-generated with timestamp)
- `--output-csv PATH` - Custom path for CSV datasets output (default: `datasets_metadata.csv` in output directory)

**Execution Options:**
- `--dry-run` - Prepare files but don't submit SLURM jobs (useful for testing)

**Logging Options:**
- `--verbose`, `-v` - Enable verbose console output
- `--log-level {DEBUG,INFO,WARNING,ERROR}` - Set logging level (default: INFO)
- `--version` - Show version number

#### Example Commands

```bash
# Standard run
python /dls/science/groups/i04-1/software/pipedream_xchem/pipedream_xchem.py --parameters params.yaml

# Test run with detailed output
python /dls/science/groups/i04-1/software/pipedream_xchem/pipedream_xchem.py --parameters params.yaml --dry-run --verbose

# Custom output locations
python /dls/science/groups/i04-1/software/pipedream_xchem/pipedream_xchem.py --parameters params.yaml --output-json /path/to/output.json --output-csv /path/to/datasets.csv
```

#### What It Does

1. Reads the YAML parameters file  
2. Queries the SQLite database for datasets
3. Validates SMILES strings and file paths
4. Prepares input folders with MTZ, PDB, and SMILES files
5. Generates SLURM array job scripts
6. Submits jobs to the cluster (unless `--dry-run`)
7. Outputs metadata to JSON and CSV files

#### 📊 Output Structure

<details>
<summary><b>🗂️ Top-Level Structure</b></summary>

```
Pipedream_YYYYMMDD_HHMMSS/
├── pipedream.log
├── datasets_metadata.csv
├── Pipedream_YYYYMMDD_HHMMSS_output.json
├── pipedream_array.sh  (or pipedream_array_chunk_*.sh for >1000 datasets)
├── array_logs/
└── CrystalName/
```

| Path | Description |
|------|--------------|
| `pipedream.log` | Main log file for the entire pipeline execution |
| `datasets_metadata.csv` | Metadata summary of processed datasets |
| `Pipedream_YYYYMMDD_HHMMSS_output.json` | Consolidated output containing file paths and commands |
| `pipedream_array.sh` | SLURM job script (or multiple chunks if >1000 datasets) |
| `array_logs/` | Contains SLURM job logs and index mappings |
| `CrystalName/` | Directory for each processed crystal |

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
| `slurm_array_index_map_*.csv` | Maps SLURM array indices to dataset identifiers |
| `pipedream_array_*.out` | Individual SLURM job logs for each array task |

</details>

<details>
<summary><b>💎 CrystalName/</b></summary>

```
CrystalName/
├── input_files/
│   ├── *.mtz
│   ├── *.pdb
│   └── CompoundCode.{cif,pdb,smiles}
├── Pipedream_YYYYMMDD_HHMMSS/    ← created by pipeline on cluster
└── CrystalName_slurm_*.out
```

| Path | Description |
|------|--------------|
| `input_files/` | Contains crystallographic data and ligand files |
| `*.mtz`, `*.pdb` | Crystallographic data and model files |
| `CompoundCode.smiles` | SMILES string for ligand (restraints generated on cluster) |
| `Pipedream_YYYYMMDD_HHMMSS/` | Pipeline-generated output directory (created during job execution) |
| `CrystalName_slurm_*.out` | SLURM job log specific to this crystal |

</details>

> 💡 **Note:** 
> - Restraint files (`.cif`, `.pdb`) are generated **on the cluster** using grade2 from the `.smiles` file
> - For datasets > 1000, jobs are automatically split into chunks of 1000 datasets each
> - The timestamped directory uniquely identifies each Pipedream run

---

### Step 2: Collate Results

```bash
python /dls/science/groups/i04-1/software/pipedream_xchem/collate_pipedream_results.py --input /path/to/<Pipedream_yyyymmdd_hhmmss>_output.json [OPTIONS]
```

#### Command-Line Options

**Required:**
- `--input`, `--json`, `-i` - Path to Pipedream output JSON file

**Output Options:**
- `--output-dir DIR` - Output directory for reports (default: same as input JSON)
- `--output-name NAME` - Base name for output files (default: `Pipedream_results_<timestamp>`)
- `--format {json,html,both}` - Output format(s) (default: both)
- `--no-browser` - Don't automatically open HTML report
- `--no-plots` - Don't automatically open plot PNG files

**Logging Options:**
- `--verbose`, `-v` - Enable verbose console output
- `--log-level {DEBUG,INFO,WARNING,ERROR}` - Set logging level (default: DEBUG)
- `--version` - Show version number

#### What It Does

This script:
- Collates results from completed Pipedream runs
- Performs chirality inversion detection
- Extracts ligand validation statistics and quality metrics
- Creates interactive HTML reports with filtering
- Outputs comprehensive JSON metadata for export

> **Note:** Maps and edstats are generated automatically during the pipeline, making collation fast.

#### Outputs

1. `Pipedream_results_yyyymmdd_hhmmss.html` - Interactive HTML report
2. `Pipedream_results_yyyymmdd_hhmmss.json` - JSON metadata for export script
3. `Pipedream_results_yyyymmdd_hhmmss.log` - Detailed log file

---

### Step 3: Export to XCE

```bash
python /dls/science/groups/i04-1/software/pipedream_xchem/export_pipedream.py --input /path/to/results.json --parameters /path/to/pipedream_parameters.yaml [OPTIONS]
```

#### Command-Line Options

**Required:**
- `--input`, `--json`, `-i` - Path to collation results JSON file
- `--parameters`, `-p` - Path to YAML parameters file

**Logging Options:**
- `--verbose`, `-v` - Enable verbose console output
- `--log-level {DEBUG,INFO,WARNING,ERROR}` - Set logging level (default: INFO)
- `--version` - Show version number

#### What It Does

This script:
- Exports validated Pipedream results to XChem Explorer (XCE) format
- Updates SQLite database with refinement statistics and quality metrics
- Applies traffic-light scoring for validation metrics
- Handles chirality inversions by copying updated restraint files (CIF/PDB)
- Generates molecular structure diagrams (PNG) from SMILES using RDKit
- Creates proper symlinks for PDB, MTZ, and map files in XCE format
- Copies BUSTER output files for downstream analysis

---

## 🛠 Troubleshooting

| Error | Cause | Fix |
|-------|--------|-----|
| `Missing required file paths` | Datasets lack MTZ/PDB paths | Check database or use `--dry-run` |
| `Could not find CIF file` | Restraints missing | Check `input_files/` or `rhofit-*/` directories |
| `SLURM job not submitted` | SSH authentication issue | Configure SSH keys to `wilson.diamond.ac.uk` |
| `RDKit not available` | RDKit not installed | `micromamba install -c conda-forge rdkit` |
| `grade2 failed` | Invalid SMILES or CSD unavailable | Verify SMILES in database, check CSD installation |
| `No valid SMILES available` | Missing SMILES in database | Update database with valid SMILES |
| `Chirality inversion detected` | Stereochemistry changed | Review output; restraints auto-copied |
| `Array job chunked` | > 1000 datasets | Normal - auto-split into chunks of 1000 |
| `Maps or edstats missing` | Post-processing failed | Check SLURM logs; reprocess with `pipedream_post_process.py --input <output.json>` |

### Additional Notes

**SSH Authentication:**
- Scripts use SSH key-based authentication (`~/.ssh/id_rsa`)
- **Not** password-based authentication
- Ensure key is set up for `wilson.diamond.ac.uk`

**Job Chunking:**
- For > 1000 datasets, jobs are automatically split into chunks
- Each chunk processes up to 1000 datasets
- This is a SLURM limitation, not an error

**YAML Requirements:**
- For `specific_datasets` mode, `Dataset_csv_path` is **required**
- CSV must have a `CrystalName` column header

**Restraint Generation:**
- Happens **on the cluster** during job execution
- Uses grade2 with SMILES from database
- Requires valid SMILES strings (checked before submission)

---

## 📬 Contact

For questions or issues, please contact:  
📧 [xchem@diamond.ac.uk](mailto:xchem@diamond.ac.uk)

---

## 📝 Version Information

- **pipedream_xchem.py**: v1.0.2
- **collate_pipedream_results.py**: v1.0.2
- **export_pipedream.py**: v1.0.2
- **pipedream_post_process.py**: v1.0.0

Check version: `python <script>.py --version`
