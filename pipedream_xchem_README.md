This python script automates the submission of Pipedream refinement jobs to the Diamond cluster as an array. It supports dataset selection from a SQLite database and prepares job scripts based on user-defined parameters in a YAML configuration file.

# **⚙️ YAML Configuration**

Create a YAML file (e.g., `pipedream_parameters.yaml`) with the following structure (or copy the template file from /dls/science/groups/i04-1/software/pipedream/pipedream_parameters.yaml):

### **Parameter Descriptions**

```yaml

Mode: "specific_datasets" # 'pending_analysis' or 'specific_datasets' - the former will parse your database file for all datasets with RefinementOutcome "1 - Pending Analysis", the later will use a specified list of datasets provided in the csv file specified below

Processing_directory: /dls/labxchem/data/proposal/visit/processing
Output_directory: /dls/labxchem/data/proposal/visit/processing/Pipedream/Pipedream_output_123 # Optional - defaults to Processing_directory/analysis/Pipedream/Pipedream_<timestamp> if not set
Database_path: /dls/labxchem/data/proposal/visit/processing/database/soakDBDataFile.sqlite
Dataset_csv_path: /dls/labxchem/data/proposal/visit/processing/Pipedream/datasets.csv # Only required if Mode is 'specific_datasets'

Copy_files_to_working_dir: true
Remove_crystallisation_components: true  # Optional - removes DMS, EDO, GOL, SO4, PO4, PEG from input PDBs if true (can skip if not modelled in site of interest in MR model)
Refinement_parameters: #For more information see https://www.globalphasing.com/buster/manual/pipedream/manual/index.html#_details_of_command_line_arguments
  keepwater: true #DO NOT remove waters that are present in the input model (default is to remove them)
  mrefine: "TLSbasic" #"TLSbasic" turns on TLS refinement and autoncs. Leave blank for no TLS.
  remediate: true #Run SideAide to refit side chains
  sidechainrebuild: true #Allow SideAide to rebuild stubbed sidechains
  runpepflip: true #Run pepflip to check for and correct peptide bond flips
  rhocommands:
	  - -xclusters 1 # Produces ligand fits for the <n> best possible binding sites. Leave blank for default and fit to <NCS> best sites.
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
    
- **Copy_files_to_working_dir**: If `true`, input files are copied into job-specific folders. **Required.**
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

# **🚀 Usage**

To run the script with the YAML file as an argument:
`micromamba-init`

`micromamba activate xchem`

`python /dls/science/groups/i04-1/software/pipedream/pipedream_xchem.py --yaml pipedream_parameters.yaml`

The script will:

1. Read the YAML configuration.
2. Query the database for datasets.
3. Prepare input folders and job scripts.
4. Submit jobs to the SLURM cluster via REST API.
5. Output data in specified directory or `<processing_directory>/analysis/<Pipedream_yyyymmdd_hhmmss>`

---

# **📝 Logging**

A log file named `pipedream.log` will be created in the output directory for each run, capturing job submission details and any errors encountered.

---

# **📈 Outputs**

Details of the pipedream output file structure can be found [here](https://www.globalphasing.com/buster/manual/pipedream/manual/index.html#_location_of_pipedream_output).

The files you are most likely looking for are:

Refine.pdb: 
`<Output_directory>/<CrystalName>/Pipedream_<timestamp>/postrefine-<CompoundCode>/refine.pdb`

Refine.mtz:
`<Output_directory>/<CrystalName>/Pipedream_<timestamp>/postrefine-<CompoundCode>/refine.mtz`

Buster-report index.html: `Output_directory>/<CrystalName>/Pipedream_<timestamp>/report-<CompoundCode>/index.html.`

Summary statistics for all datasets can be aggregated into a readable csv/html following these [instructions](https://www.notion.so/Pipedream-results-collator-21694b98aff680e6aee2c81d94b02639?pvs=21).

# **🛠 Troubleshooting**

- Ensure all required file paths in the database are valid and accessible.
- If a dataset is missing required files, it will be skipped and logged.
- Make sure your Wilson SSH credentials and SLURM token permissions are correctly configured.# pipedream_xchem
Python scripts for running pipedream on XChem data and collating results
