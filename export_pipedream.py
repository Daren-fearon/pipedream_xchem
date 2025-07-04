import os
import json
import shutil
import sqlite3
import yaml
import getpass
from datetime import datetime

# Helper function to safely create symlinks
def safe_symlink(src, dst):
    try:
        if os.path.islink(dst) or os.path.exists(dst):
            os.remove(dst)
        os.symlink(src, dst)
    except Exception as e:
        print(f"Error creating symlink from {src} to {dst}: {e}")

# Helper function to check if a column exists in a table
def column_exists(cursor, table, column):
    cursor.execute(f"PRAGMA table_info({table})")
    columns = [info[1] for info in cursor.fetchall()]
    return column in columns

# Traffic light function with corrected comparison
def traffic_light(value, green, orange=None):
    try:
        val = float(value)
        if orange:
            if val < green:
                return "green"
            elif val < orange:
                return "orange"
            else:
                return "red"
        else:
            return "green" if val < green else "red"
    except:
        return "NA"

# Load the JSON results file
json_filename = "Pipedream_results_20250630_131310.json"
with open(json_filename, "r") as f:
    results = json.load(f)

# Load the YAML parameters file
with open("pipedream_parameters.yaml", "r") as f:
    yaml_params = yaml.safe_load(f)

db_path = yaml_params["Database_path"]

# Connect to the SQLite database
conn = sqlite3.connect(db_path)
cursor = conn.cursor()

# Check required columns exist
required_columns = [
    "RefinementResolution", "RefinementRcryst", "RefinementRfree", "RefinementOutcome",
    "RefinementOutcomePerson", "RefinementOutcomeDate", "RefinementPDB_latest",
    "RefinementMTZ_latest", "RefinementMMCIFmodel_latest", "RefinementMMCIFreflections_latest",
    "RefinementLigandCC", "RefinementBoundConformation", "RefinementMolProbityScore",
    "RefinementRamachandranOutliers", "RefinementRamachandranFavored",
    "RefinementStatus", "RefinementBusterReportHTML", "RefinementRefiner",
    "RefinementDate", "LastUpdated", "LastUpdated_by", "CrystalName"
]

missing_columns = [col for col in required_columns if not column_exists(cursor, "mainTable", col)]
if missing_columns:
    print(f"Warning: Missing columns in database table: {missing_columns}")

# Get current user and timestamp
user = getpass.getuser()
now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# Process each dataset marked for export
for entry in results:
    if entry.get("Export to XCE") != "True":
        continue

    try:
        crystal = entry["Crystal Name"]
        pipedream_dir = entry["Pipedream Directory"]
        target_dir = os.path.join(os.path.dirname(pipedream_dir), "..", "..", "model_building", crystal)
        os.makedirs(target_dir, exist_ok=True)

        # Copy pipedream directory
        dest_pipedream = os.path.join(target_dir, os.path.basename(pipedream_dir))
        if not os.path.exists(dest_pipedream):
            shutil.copytree(pipedream_dir, dest_pipedream)

        # Create symlinks
        pdb_file = entry["PDB File"]
        mtz_file = entry["MTZ File"]
        bound_pdb = entry["PDB File"]

        safe_symlink(pdb_file, os.path.join(target_dir, "refine.pdb"))
        safe_symlink(mtz_file, os.path.join(target_dir, "refine.mtz"))
        safe_symlink(bound_pdb, os.path.join(target_dir, "refine.split.bound-state.pdb"))

        # Load pipedream_summary.json
        summary_path = entry["Pipedream Summary"]
        with open(summary_path, "r") as f:
            summary = json.load(f)

        ligand = summary.get("ligandfitting", {}).get("ligands", [{}])[0]
        molprobity = ligand.get("validationstatistics", {}).get("molprobity", {})
        postref = summary.get("pipedream_outputs", {}).get("ligandfitting", {}).get("ligands", [{}])[0].get("postrefinement", [])

        mmcif_model = next((x["filename"] for x in postref if x.get("type") == "model" and x.get("format") == "mmCIF"), "NA")
        mmcif_reflections = next((x["filename"] for x in postref if x.get("type") == "map" and x.get("format") == "mmCIF"), "NA")

        rama_out = molprobity.get("ramaoutlierpercent", "NA")
        rama_fav = molprobity.get("ramafavoredpercent", "NA")
        molprob = molprobity.get("molprobityscore", "NA")

        r = entry["R"]
        rfree = entry["Rfree"]
        ligandcc = entry["Ligand CC"]
        report = entry["Ligand Report HTML"]

        update_query = f"""
        UPDATE mainTable SET
            RefinementResolution = 'NA',
            RefinementResolutionTL = 'NA',
            RefinementOutcome = '3 - In Refinement',
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
            RefinementBoundConformation = ?,
            RefinementMolProbityScore = ?,
            RefinementMolProbityScoreTL = ?,
            RefinementRamachandranOutliers = ?,
            RefinementRamachandranOutliersTL = ?,
            RefinementRamachandranFavored = ?,
            RefinementRamachandranFavoredTL = ?,
            RefinementStatus = 'finished',
            RefinementBusterReportHTML = ?,
            RefinementRefiner = ?,
            RefinementDate = ?,
            LastUpdated = ?,
            LastUpdated_by = ?
        WHERE CrystalName = ?
        """

        cursor.execute(update_query, (
            r,
            traffic_light(r, 0.3),
            rfree,
            traffic_light(rfree, 0.3),
            user,
            now,
            pdb_file,
            mtz_file,
            mmcif_model,
            mmcif_reflections,
            ligandcc,
            bound_pdb,
            molprob,
            traffic_light(molprob, 2),
            rama_out,
            traffic_light(rama_out, 1, 5),
            rama_fav,
            traffic_light(rama_fav, 99, 95),
            report,
            user,
            now,
            now,
            user,
            crystal
        ))

    except Exception as e:
        print(f"Error processing crystal {entry.get('Crystal Name', 'UNKNOWN')}: {e}")

conn.commit()
conn.close()

print("Processing complete.")

