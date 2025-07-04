import os
import json
import pandas as pd
import logging
import argparse
from datetime import datetime
import webbrowser
import subprocess
from collections import OrderedDict

def read_json(json_file):
    try:
        with open(json_file, 'r') as file:
            data = json.load(file)
            logging.info("JSON file read successfully.")
            return data
    except Exception as e:
        logging.error(f"Error reading JSON file: {e}")
        raise

def validate_json_structure(json_data):
    required_keys = ['ExpectedSummary', 'LigandPNG']
    for dataset, info in json_data.items():
        for key in required_keys:
            if key not in info:
                logging.warning(f"Missing required key '{key}' in dataset '{dataset}'")

def safe_get(d, keys, default="NA"):
    for i, key in enumerate(keys):
        if isinstance(d, dict):
            d = d.get(key, {})
        elif isinstance(d, list) and isinstance(key, int) and len(d) > key:
            d = d[key]
        else:
            return default if i == len(keys) - 1 else {}
    return d if d != {} else default

def safe_round(value, digits=3):
    try:
        return round(float(value), digits)
    except (TypeError, ValueError):
        return 'NA'

def convert_mtz_to_map(mtz_file_path, map_type='2fofc'):
    map_file_path = mtz_file_path.replace('.mtz', f'_{map_type}.map')
    if os.path.exists(map_file_path):
        return map_file_path
    if map_type == '2fofc':
        subprocess.run(['gemmi', 'sf2map', mtz_file_path, map_file_path])
    elif map_type == 'fofc':
        subprocess.run(['gemmi', 'sf2map', mtz_file_path, '-d', map_file_path])
    return map_file_path

def collect_results_from_json(json_data):
    results = []
    for dataset, info in json_data.items():
        summary_file = info.get('ExpectedSummary')
        if not os.path.exists(summary_file):
            logging.warning(f"Summary file not found for {dataset}: {summary_file}")
            result = create_result(dataset, info, summary_file)
            results.append(result)
            continue

        try:
            with open(summary_file, 'r') as f:
                summary = json.load(f)
        except json.JSONDecodeError as e:
            logging.error(f"Error decoding JSON in {summary_file} for dataset {dataset}: {e}")
            result = create_result(dataset, info, summary_file)
            results.append(result)
            continue

        ligand = safe_get(summary, ['ligandfitting', 'ligands', 0])
        ligand_stats = safe_get(ligand, ['validationstatistics', 'ligandstatistics', 0])
        molprobity_stats = safe_get(ligand, ['validationstatistics', 'molprobity'])
        postrefinement_stats = safe_get(ligand, ['postrefinement', 1])

        # Extract PDB and MTZ file paths
        pdb_file = 'NA'
        mtz_file = 'NA'
        postrefinement = safe_get(summary, ['pipedream_outputs', 'ligandfitting', 'ligands', 0, 'postrefinement'], [])
        if isinstance(postrefinement, list):
            for entry in postrefinement:
                if entry.get('description') == 'final':
                    if entry.get('type') == 'model' and entry.get('format') == 'PDB':
                        pdb_file = os.path.join(info.get('PipedreamDirectory', 'NA'), entry.get('relative_path', ''), entry.get('filename', ''))
                    elif entry.get('type') == 'map' and entry.get('format') == 'MTZ':
                        mtz_file = os.path.join(info.get('PipedreamDirectory', 'NA'), entry.get('relative_path', ''), entry.get('filename', ''))

        # Generate map files if MTZ file is valid
        map_2fofc_file = convert_mtz_to_map(mtz_file, '2fofc') if mtz_file != 'NA' else 'NA'
        map_fofc_file = convert_mtz_to_map(mtz_file, 'fofc') if mtz_file != 'NA' else 'NA'

        result = {
            'Crystal Name': dataset,
            'Compound Code': info.get('CompoundCode', 'NA'),
            'Ligand Structure': info.get('LigandPNG', 'NA'),
            'Pipedream Directory': info.get('PipedreamDirectory', 'NA'),
            'Buster Report HTML': info.get('ReportHTML', 'NA'),
            'Ligand Report HTML': info.get('LigandReportHTML', 'NA'),
            'Pipedream Summary': summary_file,
            'PDB File': pdb_file,
            'MTZ File': mtz_file,
            '2Fo-Fc Map File': map_2fofc_file,
            'Fo-Fc Map File': map_fofc_file,
            'Ligand ID': ligand_stats.get('ligandid', 'NA') if isinstance(ligand_stats, dict) else 'NA',
            'Ligand CC': ligand_stats.get('ligandcc', 'NA') if isinstance(ligand_stats, dict) else 'NA',
            'R': safe_round(postrefinement_stats.get('R')) if isinstance(postrefinement_stats, dict) else 'NA',
            'Rfree': safe_round(postrefinement_stats.get('Rfree')) if isinstance(postrefinement_stats, dict) else 'NA',
            'Ligand avg B factor': ligand_stats.get('ligandbavg', 'NA') if isinstance(ligand_stats, dict) else 'NA',
            'Ligand occupancy': ligand_stats.get('ligandomin', 'NA') if isinstance(ligand_stats, dict) else 'NA',
            'Mogul Z angle': ligand_stats.get('mogulzangl', 'NA') if isinstance(ligand_stats, dict) else 'NA',
            'Mogul Z bond': ligand_stats.get('mogulzbond', 'NA') if isinstance(ligand_stats, dict) else 'NA',
            'c beta deviations': molprobity_stats.get('cbetadeviations', 'NA') if isinstance(molprobity_stats, dict) else 'NA',
            'Rama outlier percent': molprobity_stats.get('ramaoutlierpercent', 'NA') if isinstance(molprobity_stats, dict) else 'NA',
            'Rama favored percent': molprobity_stats.get('ramafavoredpercent', 'NA') if isinstance(molprobity_stats, dict) else 'NA',
            'Poor rotamer percent': molprobity_stats.get('poorrotamerspercent', 'NA') if isinstance(molprobity_stats, dict) else 'NA',
            'Clash score': molprobity_stats.get('clashscore', 'NA') if isinstance(molprobity_stats, dict) else 'NA',
            'Mol probity score': molprobity_stats.get('molprobityscore', 'NA') if isinstance(molprobity_stats, dict) else 'NA',
            'RMS bonds': molprobity_stats.get('rmsbonds', 'NA') if isinstance(molprobity_stats, dict) else 'NA',
            'RMS angles': molprobity_stats.get('rmsangles', 'NA') if isinstance(molprobity_stats, dict) else 'NA',
            'Mean B factor': postrefinement_stats.get('MeanB', 'NA') if isinstance(postrefinement_stats, dict) else 'NA',
            'Comments': '',
            'Export to XCE': 'True' if (isinstance(ligand_stats, dict) and ligand_stats.get('ligandcc', 0) > 0.8) else 'False'
        }
        results.append(result)
    return results

def create_result(dataset, info, summary_file):
    return {
        'Crystal Name': dataset,
        'Compound Code': info.get('CompoundCode', 'NA'),
        'Ligand Structure': info.get('LigandPNG', 'NA'),
        'Pipedream Directory': info.get('PipedreamDirectory', 'NA'),
        'Buster Report HTML': info.get('ReportHTML', 'NA'),
        'Ligand Report HTML': info.get('LigandReportHTML', 'NA'),
        'Pipedream Summary': summary_file,
        'PDB File': 'NA',
        'MTZ File': 'NA',
        '2Fo-Fc Map File': 'NA',
        'Fo-Fc Map File': 'NA',
        'Ligand ID': 'NA',
        'Ligand CC': 'NA',
        'R': 'NA',
        'Rfree': 'NA',
        'Ligand avg B factor': 'NA',
        'Ligand occupancy': 'NA',
        'Mogul Z angle': 'NA',
        'Mogul Z bond': 'NA',
        'c beta deviations': 'NA',
        'Rama outlier percent': 'NA',
        'Rama favored percent': 'NA',
        'Poor rotamer percent': 'NA',
        'Clash score': 'NA',
        'Mol probity score': 'NA',
        'RMS bonds': 'NA',
        'RMS angles': 'NA',
        'Mean B factor': 'NA',
        'Comments': '',
        'Export to XCE': 'False'
    }

def save_results_to_csv(results, output_file):
    # Define the desired column order
    column_order = [
        'Export to XCE',
        'Comments',
        'Crystal Name',
        'Compound Code',
        'Ligand Structure',
        'Pipedream Directory',
        'Buster Report HTML',
        'Ligand Report HTML',
        'Ligand CC',
        'Ligand occupancy',
        'Ligand avg B factor',
        'Mean B factor',
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
        'Fo-Fc Map File'
    ]
    df = pd.DataFrame(results)
    df = df.reindex(columns=column_order)
    df.to_csv(output_file, index=False)
    logging.info(f"Results saved to {output_file}")

def save_results_to_html(results, output_file):
    df = pd.DataFrame(results)

    # Ensure 'Ligand CC' is numeric for sorting, non-numeric values become NaN
    df['Ligand CC'] = pd.to_numeric(df['Ligand CC'], errors='coerce')
    # Auto sort by Ligand CC descending (NaNs will be at the end)
    df = df.sort_values(by='Ligand CC', ascending=False)

    # Define the desired column order (ensure it's always defined)
    column_order = [
        'Export to XCE',
        'Comments',
        'Crystal Name',
        'Compound Code',
        'Ligand Structure',
        'Pipedream Directory',
        'Buster Report HTML',
        'Ligand Report HTML',
        'Ligand CC',
        'Ligand occupancy',
        'Ligand avg B factor',
        'Mean B factor',
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
        'Fo-Fc Map File'
    ]

    # Build table rows with actual cell values in the new order
    table_rows = []
    for _, row in df.iterrows():
        row_html = '<tr>'
        for col in column_order:
            cell_val = row[col]
            style = ''
            display_val = cell_val
            # Try to format numbers to 2 decimal places if possible
            try:
                num_val = float(cell_val)
                display_val = f"{num_val:.2f}"
            except (ValueError, TypeError):
                pass
            # Highlighting rules
            if col == 'Export to XCE':
                row_html += f'<td><input type="checkbox" {"checked" if cell_val == "True" else ""}></td>'
            elif col == 'Comments':
                row_html += f'<td><input type="text" value="{cell_val}"></td>'
            elif col == 'Ligand Structure':
                row_html += f'<td><img src="file://{cell_val}" alt="Ligand Image" width="150"></td>' if cell_val and cell_val != 'NA' else '<td></td>'
            elif col in [
                'Pipedream Directory', 'Buster Report HTML', 'Ligand Report HTML', 'Pipedream Summary',
                'PDB File', 'MTZ File', '2Fo-Fc Map File', 'Fo-Fc Map File']:
                if cell_val and cell_val != 'NA':
                    file_ext = os.path.splitext(cell_val)[1].lower()
                    # Use a data attribute for the file path
                    row_html += f'<td><a href="#" class="file-link" data-file="file://{cell_val}" data-ext="{file_ext}">{os.path.basename(cell_val)}</a></td>'
                else:
                    row_html += '<td></td>'
            else:
                # Numeric highlighting rules
                try:
                    val = float(cell_val)
                    if col == 'Rfree' and val > 0.3:
                        style = 'background-color: #ffcccc; font-weight: bold;'
                    elif col == 'R' and val > 0.3:
                        style = 'background-color: #ffcccc; font-weight: bold;'
                    elif col == 'Mogul Z angle' and val > 2:
                        style = 'background-color: #ffcccc; font-weight: bold;'
                    elif col == 'Mogul Z bond' and val > 2:
                        style = 'background-color: #ffcccc; font-weight: bold;'
                    elif col == 'Ligand avg B factor':
                        try:
                            mean_b = float(row['Mean B factor'])
                            if mean_b != 0 and val / mean_b >= 1.5:
                                style = 'background-color: #ffcccc; font-weight: bold;'
                        except (ValueError, TypeError, ZeroDivisionError):
                            pass
                    elif col == 'Rama outlier percent' and val > 0.5:
                        style = 'background-color: #ffcccc; font-weight: bold;'
                    elif col == 'Poor rotamer percent' and val > 2:
                        style = 'background-color: #ffcccc; font-weight: bold;'
                    elif col == 'Clash score' and val > 20:
                        style = 'background-color: #ffcccc; font-weight: bold;'
                    elif col == 'Mol probity score' and val > 2:
                        style = 'background-color: #ffcccc; font-weight: bold;'
                    elif col == 'c beta deviations' and val > 0:
                        style = 'background-color: #ffcccc; font-weight: bold;'
                    elif col == 'RMS angles' and val > 3:
                        style = 'background-color: #ffcccc; font-weight: bold;'
                    elif col == 'RMS bonds' and val > 0.03:
                        style = 'background-color: #ffcccc; font-weight: bold;'
                except (ValueError, TypeError):
                    pass
                row_html += f'<td style="{style}">{display_val}</td>'
        row_html += '</tr>'
        table_rows.append(row_html)

    # Render table header HTML outside JS template
    table_header_html = ''.join([f'<th>{col}</th>' for col in column_order])

    html_content = f"""
    <!DOCTYPE html>
    <html lang=\"en\">
    <head>
        <meta charset=\"UTF-8\">
        <title>Pipedream XChem Results</title>
        <link rel=\"stylesheet\" href=\"https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css\">
        <link rel=\"stylesheet\" href=\"https://cdn.datatables.net/1.13.4/css/dataTables.bootstrap5.min.css\">
        <script src=\"https://code.jquery.com/jquery-3.6.0.min.js\"></script>
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
        </style>
        <script>
            // Embed the original JSON filename for export
            const originalJsonFilename = "{os.path.basename(output_file).replace('.html', '.json')}";
            $(document).ready(function() {{
                $('#resultsTable').DataTable({{
                    responsive: true,
                    pageLength: 25,
                    order: [[8, 'desc']], // Ligand CC column index (0-based)
                    columnDefs: [
                        {{ targets: '_all', type: 'num' }}
                    ],
                    scrollX: true,
                    scrollY: '60vh',
                    scrollCollapse: true,
                    fixedHeader: true,
                    dom: 'Bflrt<"bottom-scrollbar"ip>',
                }});
                // Only sync scroll for bottom
                function syncScroll() {{
                    var scrollBody = $('.dataTables_scrollBody')[0];
                    var scrollFoot = $('.dataTables_scrollFoot')[0];
                    if(scrollBody && scrollFoot) {{
                        scrollFoot.scrollLeft = scrollBody.scrollLeft;
                    }}
                }}
                $('.dataTables_scrollBody').on('scroll', syncScroll);

                // Modal popup for file links
                $(document).on('click', '.file-link', function(e) {{
                    e.preventDefault();
                    var fileUrl = $(this).data('file');
                    var ext = $(this).data('ext');
                    var modalBody = $('#fileModal .modal-body');
                    modalBody.empty();
                    if(ext === '.html') {{
                        modalBody.append(`<iframe src="${{fileUrl}}" style="width:100%;height:70vh;border:none;"></iframe>`);
                    }} else if(ext === '.pdb' || ext === '.mtz' || ext === '.map') {{
                        modalBody.append(`<p>File: <a href="${{fileUrl}}" download>Download</a></p>`);
                    }} else {{
                        modalBody.append(`<p>File: <a href="${{fileUrl}}" target="_blank">Open</a></p>`);
                    }}
                    $('#fileModal').modal('show');
                }});
            }});
        </script>
        <!-- Bootstrap Modal -->
        <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">
        <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js"></script>
    </head>
    <body>
        <div class="container">
            <div class="mb-3 d-flex flex-wrap gap-2 align-items-center">
                <input type="file" id="jsonFileInput" accept="application/json" style="display:none" onchange="loadJsonAndUpdateTable(event)">
                <button class="btn btn-warning" onclick="document.getElementById('jsonFileInput').click()">Load JSON and Update Table</button>
                <button class="btn btn-primary" onclick="exportToCSV()">Export Selected to CSV</button>
                <button class="btn btn-success" onclick="exportToJSON()">Export Table to JSON</button>
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
              <div class="modal-body">
              </div>
            </div>
          </div>
        </div>
    </body>
    </html>
    """

    with open(output_file, 'w') as f:
        f.write(html_content)

    logging.info(f"Results saved to {output_file}")
    webbrowser.open(f'file://{os.path.abspath(output_file)}')

def save_results_to_json(results, output_file):
    # Use the same column order as for CSV/HTML
    column_order = [
        'Export to XCE',
        'Comments',
        'Crystal Name',
        'Compound Code',
        'Ligand Structure',
        'Pipedream Directory',
        'Buster Report HTML',
        'Ligand Report HTML',
        'Ligand CC',
        'Ligand occupancy',
        'Ligand avg B factor',
        'Mean B factor',
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
        'Fo-Fc Map File'
    ]
    ordered_results = [OrderedDict((col, row.get(col, '')) for col in column_order) for row in results]
    with open(output_file, 'w') as f:
        json.dump(ordered_results, f, indent=2)
    logging.info(f"Results saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Collect results from the latest Pipedream_N directories.')
    parser.add_argument('--json', required=True, help='Path to the JSON file with dataset metadata.')
    parser.add_argument('--output_file', help='Optional: specify output file name (with .json extension) for results. Defaults to timestamped file in same dir.')
    args = parser.parse_args()

    output_dir = os.path.dirname(args.json)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    if args.output_file:
        output_file_json = args.output_file
        output_file_html = args.output_file.replace('.json', '.html')
    else:
        output_file_json = os.path.join(output_dir, f'Pipedream_results_{timestamp}.json')
        output_file_html = os.path.join(output_dir, f'Pipedream_results_{timestamp}.html')
    log_file = os.path.join(output_dir, 'collate_pipedream_results.log')

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )

    json_data = read_json(args.json)
    validate_json_structure(json_data)
    results = collect_results_from_json(json_data)
    save_results_to_json(results, output_file_json)
    save_results_to_html(results, output_file_html)
    logging.info("Result collection completed.")

if __name__ == "__main__":
    main()
