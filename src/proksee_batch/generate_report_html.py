import os
from typing import Any
from typing import Dict
from typing import List


#  <!doctype html>
#  <html lang="en">
#    <head>
#      <!-- Required meta tags -->
#      <meta charset="utf-8">
#      <meta name="viewport" content="width=device-width, initial-scale=1">
#
#      <!-- Bootstrap CSS -->
#      <link href="./styles/bootstrap.min.css" rel="stylesheet">
#
#      <link rel="stylesheet" href="./styles/prism.css" />
#      <link rel="stylesheet" href="./styles/tables.css" />
#      <link rel="stylesheet" href="./styles/general.css" />
#      <script src='./scripts/marked.min.js'></script>
#      <script src='./scripts/general.js'></script>
#
#      <!-- D3 -->
#      <script src="./scripts/d3.min.js"></script>
#      <!-- CGView -->
#      <script src='./dist/cgview.min.js'></script>
#      <link rel="stylesheet" href="./dist/cgview.css" />
#      <link rel="stylesheet" href="./styles/controls.css" />
#
#      <title>CGView.js - Examples</title>
#    </head>
#    <body>
#
#  <body>
#    <main>
#
#  <div id='my-viewer'></div>
#  <div class='cgv-controls'>
#    <div class='cgv-btn' id='btn-reset' title='Reset Map'></div>
#    <div class='cgv-btn' id='btn-zoom-in' title='Zoom In'></div>
#    <div class='cgv-btn' id='btn-zoom-out' title='Zoom Out'></div>
#    <div class='cgv-btn' id='btn-move-left' title='Move Left/Counterclockwise'></div>
#    <div class='cgv-btn' id='btn-move-right' title='Move Right/Clockwise'></div>
#    <div class='cgv-btn' id='btn-toggle-format' title='Toggle Linear/Circular Format'></div>
#    <div class='cgv-btn' id='btn-invert-colors' title='Invert Map Colors'></div>
#    <div class='cgv-btn' id='btn-download' title='Download Map PNG'></div>
#    <div class='cgv-btn' id='btn-toggle-labels' title='Toggle Labels'></div>
#  </div>
#
#  </main>
#
#      <script src="./scripts/prism.js"></script>
#
#      <!-- Bootstrap JavaScript -->
#      <script src="./scripts/bootstrap.min.js"></script>
#
#      <script>
#        function createViewerAndLoadJS(path) {
#          // Create Viewer in default div: #my-viewer
#          const cgv = new CGV.Viewer('#my-viewer', {height: 500});
#
#          // Auto resize viewer
#          autoResizeMyViewer();
#
#          // Add viewer as global variable 'cgv'
#          window.cgv = cgv;
#
#          // Load the JavaScript file
#          var script = document.createElement('script');
#          script.src = path;
#          script.onload = function() {
#            // Use the 'json' variable directly
#            cgv.io.loadJSON(window.json);
#            cgv.draw();
#          };
#          document.head.appendChild(script);
#        }
#
#        createViewerAndLoadJS('./data/js/NZ_CP028842.js');
#      </script>
#      <script src='./scripts/controls.js'></script>
#
#    </body>
#  </html>


def generate_report_html(
    output_dir: str,
    genome_info: Dict[str, Any],
) -> None:
    """
    Generates an HTML report file with a table containing links to Proksee
    projects and images for each sample.
    """
    assert os.path.isdir(output_dir), f"Output directory does not exist: {output_dir}"
    output_file = os.path.join(output_dir, "report.html")
    with open(output_file, "w") as file:
        # Write the DOCTYPE, html, and head sections with CSS
        file.write(
            """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">

    <!-- Import CGView.js code -->
    <link href="./cgview-js_code/docs/styles/bootstrap.min.css" rel="stylesheet">
    <link rel="stylesheet" href="./cgview-js_code/docs/styles/prism.css" />
    <link rel="stylesheet" href="./cgview-js_code/docs/styles/tables.css" />
    <link rel="stylesheet" href="./cgview-js_code/docs/styles/general.css" />
    <script src='./cgview-js_code/docs/scripts/marked.min.js'></script>
    <script src='./cgview-js_code/docs/scripts/general.js'></script>
    <script src="./cgview-js_code/docs/scripts/d3.min.js"></script>
    <script src='./cgview-js_code/docs/dist/cgview.min.js'></script>
    <link rel="stylesheet" href="./cgview-js_code/docs/dist/cgview.css" />
    <link rel="stylesheet" href="./cgview-js_code/docs/styles/controls.css" />
    <script src="./cgview-js_code/docs/scripts/prism.js"></script>
    <script src="./cgview-js_code/docs/scripts/bootstrap.min.js"></script>

    <title>Results</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            line-height: 1.6;
            color: #333;
            background-color: #f4f4f4;
        }
        h1 {
            font-family: Georgia, serif;
        }
        table {
            width: 100%;
            border-collapse: collapse;
        }
        th, td {
            border: 1px solid #ddd;
            padding: 10px;
            text-align: left;
        }
        td {
            padding: 10px;
            vertical-align: middle;
        }
        th {
            background-color: #eaeaea;
        }
        th {
            background-color: #eaeaea;
            cursor: pointer; /* Adds the pointer cursor on hover */
        }
        th:hover {
            color: #0275d8; /* Changes text color on hover */
            font-weight: bold; /* Makes the font bold on hover */
        }
        tr:nth-child(even) {
            background-color: #f9f9f9;
        }
        .dropdown-arrow {
            cursor: pointer;
            user-select: none;
        }
        .viewer-row {
            display: none;
        }
        .image-cell img {
            max-width: 100%;
            max-height: 600px;
            object-fit: contain;
        }
        @media screen and (max-width: 600px) {
            table, thead, tbody, th, td, tr {
                display: block;
            }
        }
        .generate-link {
            color: blue;
            cursor: pointer;
            text-decoration: underline;
        }
        .generated-link {
            color: green;
            text-decoration: none;
        }
    </style>
</head>
"""
        )

        # Write the body section with a dynamic table
        file.write(
            """
    <body>
        <h1>Proksee Batch results</h1>
        <table id="resultsTable">
            <thead>
                <tr>
                    <th onclick="sortTable(0)"></th>
                    <th onclick="sortTable(1)">Sample ID</th>
                    <th onclick="sortTable(2)">Description</th>
                    <th onclick="sortTable(3)">Total size (bp)</th>
                    <th onclick="sortTable(4)">Contigs</th>
                    <th onclick="sortTable(5)">GC content</th>
                    <th>Action</th>
                </tr>
            </thead>
            <tbody>
                """
        )

        ## Create table rows based on the js and svg files
        for genome_code_name in genome_info.keys():
            genome_name = genome_info[genome_code_name]["Name"]
            description = genome_info[genome_code_name]["Description"]
            total_size = str(genome_info[genome_code_name]["Total size"])
            number_of_contigs = str(genome_info[genome_code_name]["Number of contigs"])
            gc_content = str(genome_info[genome_code_name]["GC content"])

            file.write(
                f"""
        <tr>
            <td class="dropdown-arrow">&#9660;</td>
            <td>{genome_name}</td>
            <td>{description}</td>
            <td>{total_size}</td>
            <td>{number_of_contigs}</td>
            <td>{gc_content}</td>
            <td class="generate-link" onclick="generateProkseeLink(this, 'data/{genome_code_name}')">Generate Proksee Project</td>
        </tr>
        <tr class="viewer-row">
            <td colspan="7">
                <div id="{genome_code_name}" class="genome-viewer"></div>
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
                </div>
            </td>
        </tr>


"""
            )

        # Write the script section
        file.write(
            """
    </table>

    <script>
        function sortTable(column) {
            var table, rows, switching, i, x, y, shouldSwitch;
            table = document.getElementById("resultsTable");
            switching = true;

            while (switching) {
                switching = false;
                rows = table.getElementsByTagName("TR");

                // Skip the header row and iterate by 2 since each data row is followed by an image row
                for (i = 1; i < rows.length - 3; i += 2) {
                    shouldSwitch = false;

                    x = rows[i].getElementsByTagName("TD")[column];
                    y = rows[i + 2].getElementsByTagName("TD")[column]; // Compare with the next data row

                    if (column < rows[i].cells.length) { // Ensure we are not comparing the last column
                        if (x.innerHTML.toLowerCase() > y.innerHTML.toLowerCase()) {
                            shouldSwitch = true;
                            break;
                        }
                    }
                }

                if (shouldSwitch) {
                    // Move the data row and its following image row together
                    rows[i].parentNode.insertBefore(rows[i + 2], rows[i]);
                    rows[i].parentNode.insertBefore(rows[i + 3], rows[i + 1]); // Move the image row
                    switching = true;
                }
            }
        }

        function loadScript(scriptUrl, callback) {
            const existingScript = document.querySelector(`script[src="${scriptUrl}"]`);
            if (existingScript) {
                existingScript.remove();
            }
            const script = document.createElement('script');
            script.src = scriptUrl;
            script.type = 'text/javascript';
            document.body.appendChild(script);
            script.onload = function() {
                console.log(`Script loaded: ${scriptUrl}`);
                if (typeof callback === "function") {
                    callback();
                }
            };
            script.onerror = function() {
                console.error(`Error loading script: ${scriptUrl}`);
                alert(`Failed to load data for script: ${scriptUrl}`);
            };
        }

        function generateProkseeLink(element, sampleId) {
            loadScript(`${sampleId}.js`, function() {
                if (typeof window.jsonData === 'undefined') {
                    console.error('No data found for sample ID:', sampleId);
                    alert('Failed to generate Proksee project: No data available for this sample.');
                    return;
                }
                const data = { origin: 'proksee-batch', data: JSON.stringify(window.jsonData) };
                const url = 'https://proksee.ca/api/v1/projects.json';
                fetch(url, {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify(data)
                })
                .then(response => response.json())
                .then(data => {
                    if (data?.status === 'success' && data?.url) {
                        const link = document.createElement('a');
                        link.textContent = 'Link to Proksee Project';
                        link.href = data.url;
                        link.className = 'generated-link';
                        link.target = '_blank';
                        element.parentNode.replaceChild(link, element);
                    } else {
                        console.error(`Failed to create Proksee project: ${data?.error}`);
                    }
                })
                .catch(error => {
                    console.error('Error:', error);
                });
            });
        }


        function createViewerAndLoadJS(viewer_id, path) {
            // Create Viewer
            const cgv = new CGV.Viewer('#' + viewer_id, {height: 500});

            // Auto resize viewer
            autoResizeMyViewer();

            // Add viewer as global variable 'cgv'
            window.cgv = cgv;

            // Load the JavaScript file
            var script = document.createElement('script');
            script.src = path;
            script.onload = function() {
                // Use the 'jsonData' variable directly
                cgv.io.loadJSON(window.jsonData);
                cgv.draw();
            };
            document.head.appendChild(script);
        }

        document.querySelectorAll('.dropdown-arrow').forEach(arrow => {
            arrow.addEventListener('click', function() {
                let nextRow = this.parentNode.nextElementSibling;
                if (nextRow && nextRow.classList.contains('viewer-row')) {
                    nextRow.style.display = nextRow.style.display === 'table-row' ? 'none' : 'table-row';
                    this.innerHTML = nextRow.style.display === 'table-row' ? '&#9650;' : '&#9660;';
                    // Load the viewer if it is not already loaded
                    if (nextRow.style.display === 'table-row') {
                        const viewerId = nextRow.querySelector('.genome-viewer').id;
                        const viewer = document.getElementById(viewerId);
                        if (viewer && viewer.children.length === 0) {
                            createViewerAndLoadJS(viewerId, `./data/${viewerId}.js`);
                        }
                    }
                }
            });
        });


    </script>
</body>
</html>
"""
        )
