# Pipedream XChem - Quick Reference

**Complete pipeline workflow for automated Pipedream refinement at Diamond Light Source**

---

## 📋 Standard Workflow

```bash
# 1. Setup environment
micromamba activate xchem

# 2. Submit refinement jobs
python /dls/science/groups/i04-1/software/pipedream_xchem/pipedream_xchem.py \
  --parameters pipedream_parameters.yaml

# 3. Wait for jobs to complete, then collate results
python /dls/science/groups/i04-1/software/pipedream_xchem/collate_pipedream_results.py \
  --input Pipedream_YYYYMMDD_HHMMSS_output.json

# 4. Export to XChem Explorer
python /dls/science/groups/i04-1/software/pipedream_xchem/export_pipedream.py \
  --input Pipedream_results_YYYYMMDD_HHMMSS.json \
  --parameters pipedream_parameters.yaml
```

---

## ⚙️ YAML Configuration

**Minimal example:**
```yaml
Mode: "specific_datasets"
Processing_directory: /dls/labxchem/data/proposal/visit/processing
Database_path: /path/to/soakDBDataFile.sqlite
Dataset_csv_path: /path/to/datasets.csv

Cluster_partition: "cs05r"  # Optional: cs05r or cs04r
Job_priority: "low"          # Optional: normal, low, high

Refinement_parameters:
  keepwater: true
  TLS: "TLSbasic"
  remediate: true
```

**Key parameters:**
- `Mode`: `"pending_analysis"` or `"specific_datasets"`
- `Dataset_csv_path`: Required for `"specific_datasets"` mode
- `Cluster_partition`: Optional cluster selection (`cs05r` or `cs04r`)
- `Job_priority`: Optional priority (`normal`, `low`, `high`) - `low` uses nice=1000 for background jobs
- `Remove_crystallisation_components`: Optional, removes DMS/EDO/GOL/SO4/PO4/PEG
- `Output_directory`: Optional, defaults to timestamped folder

---

## 🚀 Common Commands

### Test run (no submission)
```bash
python pipedream_xchem.py --parameters params.yaml --dry-run --verbose
```

### Monitor SLURM jobs
```bash
squeue -u $USER
tail -f array_logs/pipedream_array_*.out
```

### Check post-processing
```bash
# Verify maps and edstats were generated
ls -lh CrystalName/Pipedream_*/postrefine-*/refine_*.map
cat CrystalName/Pipedream_*/postrefine-*/edstats.out
```

### Reprocess old data
```bash
# For outputs before automatic post-processing
python pipedream_post_process.py --json old_Pipedream_output.json
```

---

## 📁 Output Structure

```
Pipedream_YYYYMMDD_HHMMSS/
├── pipedream.log
├── datasets_metadata.csv
├── Pipedream_YYYYMMDD_HHMMSS_output.json
├── pipedream_array.sh
├── array_logs/
│   ├── slurm_array_index_map_*.csv
│   └── pipedream_array_*.out
└── CrystalName/
    ├── input_files/
    │   ├── *.mtz, *.pdb
    │   └── CompoundCode.smiles
    ├── Pipedream_YYYYMMDD_HHMMSS/
    │   └── postrefine-CompoundCode/
    │       ├── refine.mtz, refine.pdb
    │       ├── refine_2fofc.map  ← auto-generated
    │       ├── refine_fofc.map   ← auto-generated
    │       └── edstats.out       ← auto-generated
    └── CrystalName_slurm_*.out
```

---

## 🔍 Checking Status

### Job submission
```bash
# Check SLURM queue
squeue -u $USER

# View array job logs
tail -f array_logs/pipedream_array_1_*.out

# Check specific dataset
tail -f CrystalName/CrystalName_slurm_*.out
```

### Results verification
```bash
# Count completed datasets
grep -c "successfully" pipedream.log

# Check for errors
grep -i "error\|failed" pipedream.log

# Verify maps generated
find . -name "refine_2fofc.map" | wc -l
```

---

## 🐛 Common Issues

| Problem | Solution |
|---------|----------|
| `Missing required file paths` | Check database has MTZ/PDB paths |
| `SLURM job not submitted` | Verify SSH key setup to wilson.diamond.ac.uk |
| `grade2 failed` | Check SMILES validity in database |
| `No valid SMILES` | Update database (avoid 'none', 'null', 'nan') |
| `Chirality inversion` | Review and confirm; restraints auto-updated |
| `Maps missing` | Reprocess: `pipedream_post_process.py --json <file>` |

---

## 📊 Key Files

| File | Purpose |
|------|---------|
| `pipedream_xchem.py` | Submit jobs with auto post-processing |
| `collate_pipedream_results.py` | Aggregate results, HTML report |
| `export_pipedream.py` | Export to XCE, update database |
| `pipedream_post_process.py` | Reprocess old data (backwards compat) |
| `pipedream_parameters.yaml` | Configuration file |
| `datasets.csv` | Dataset list for specific_datasets mode |

---

## 💡 Tips

- **Use `--dry-run`** to test before submitting
- **Check logs** for errors before collation
- **Maps are auto-generated** during pipeline (87% faster collation)
- **Jobs auto-chunk** at 1000 datasets (normal behavior)
- **SSH keys required** for job submission to wilson.diamond.ac.uk

---

## 📞 Support

- **Full documentation:** See `README.md`
- **Email:** xchem@diamond.ac.uk
- **Versions:** Run `python <script>.py --version`

---

*Quick Reference v2.0 | Updated: November 24, 2025*
