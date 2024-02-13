"""Code for generating an HTML report file with a table containing links to
Proksee projects and images for each sample. A single genome viewer is
positioned to the right of the table.

The output directory will be structured as in the following example:

    output_directory/
        cgview-js_code/
            ...
        html_report_code/
            style.css
            table-functions.js
            viewer-functions.js
            utilities.js
        data/
            genome_name_1.js
            genome_name_2.js
            ...
        report.html
"""

import os
from typing import Any
from typing import Dict


def generate_report_html(output_dir: str, genome_info: Dict[str, Any]) -> None:
    """
    Generates an HTML report file with a table containing links to Proksee
    projects and images for each sample. A single genome viewer is positioned
    to the right of the table.
    """
    assert os.path.isdir(output_dir), f"Output directory does not exist: {output_dir}"
    output_file = os.path.join(output_dir, "report.html")

    with open(output_file, "w") as file:
        # Write the DOCTYPE, html, head sections with CSS, and JavaScript imports
        file.write(
            """<!DOCTYPE html>
<html lang="en">
  <head>

    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">

    <!-- Import CGView.js code -->
    <link href="./cgview-js_code/docs/styles/bootstrap.min.css" rel="stylesheet">
    <script src='./cgview-js_code/docs/scripts/marked.min.js'></script>
    <script src='./cgview-js_code/docs/scripts/general.js'></script>
    <script src="./cgview-js_code/docs/scripts/d3.min.js"></script>
    <script src='./cgview-js_code/docs/dist/cgview.min.js'></script>
    <link rel="stylesheet" href="./cgview-js_code/docs/dist/cgview.css" />
    <link rel="stylesheet" href="./cgview-js_code/docs/styles/controls.css" />
    <link rel="stylesheet" href="./html_report_code/style.css">

    <title>Proksee Batch</title>

  </head>
<body>
    <header>
        <h1>Proksee Batch</h1>
    </header>

    <main>
        <div class="container">
        <div class="side-table">
            <div class="scrollable-table">
            <table id="sortable-table">
                <tr>
                    <th onclick="sortTable(0)">Sample ID</th>
                    <th onclick="sortTable(1)">Total size (bp)</th>
                    <th onclick="sortTable(2)">Contigs</th>
                    <th onclick="sortTable(3)">GC content</th>
                    <th>Action</th>
                </tr>
        """
        )

        # Create table rows based on the js and svg files
        for genome_code_name, info in genome_info.items():
            genome_name = info["Name"]
            # description = info["Description"]
            total_size = str(info["Total size"])
            number_of_contigs = str(info["Number of contigs"])
            gc_content = str(info["GC content"])

            file.write(
                f"""
                    <tr id='data/{genome_code_name}.js'>
                        <td>{genome_name}</td>
                        <td>{total_size}</td>
                        <td>{number_of_contigs}</td>
                        <td>{gc_content}</td>
                        <td class="generate-link" onclick="generateProkseeLink(this, 'data/{genome_code_name}')">Generate Proksee Project</td>
                    </tr>
            """
            )

        # Close the table.
        file.write(
            """
                </table>
            </div>
        </div>
        """
        )

        # Add the genome viewer.
        file.write(
            """
        <!-- Genome Viewer and Controls Container -->
        <div class="viewer-container">
            <!-- Genome Viewer -->
            <div id='my-viewer'></div>
            <!-- Controls -->
            <div class='cgv-controls'>
            <div class='cgv-btn' id='btn-reset' title='Reset Map'></div>
            <div class='cgv-btn' id='btn-zoom-in' title='Zoom In'></div>
            <div class='cgv-btn' id='btn-zoom-out' title='Zoom Out'></div>
            <div class='cgv-btn' id='btn-move-left' title='Move Left/Counterclockwise'></div>
            <div class='cgv-btn' id='btn-move-right' title='Move Right/Clockwise'></div>
            <div class='cgv-btn' id='btn-toggle-format' title='Toggle Linear/Circular Format'></div>
            <div class='cgv-btn' id='btn-invert-colors' title='Invert Map Colors'></div>
            <div class='cgv-btn' id='btn-download' title='Download Map PNG'></div>
            <div class='cgv-btn' id='btn-toggle-labels' title='Toggle Labels'></div>
            <div class='cgv-btn' id='btn-toggle-legend' title='Toggle Legend'></div>
            </div>
        </div>
        """
        )

        # Load more scripts and close the main content.
        file.write(
            """
        <script src="./cgview-js_code/docs/scripts/bootstrap.min.js"></script>
        <script src="./cgview-js_code/docs/scripts/controls.js"></script>
        <script src='./html_report_code/table-functions.js'></script>
        <script src='./html_report_code/viewer-functions.js'></script>
        <script src='./html_report_code/utilities.js'></script>

    </div>
  </main>
  </body>
</html>"""
        )
