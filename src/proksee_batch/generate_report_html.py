"""Code for generating an HTML report file with a table containing links to
Proksee projects and images for each sample. A single genome viewer is
positioned to the right of the table.

The output directory will be structured as in the following example:

    output_directory/
        cgview-js_code/
            ...
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

    <title>Proksee Batch</title>

    <style>
      header {
        background-color: #f8f9fa;
        padding: 10px;
        text-align: left;
      }

      header h1 {
        margin: 0;
        color: #333;
      }

      body, html {
        font-family: 'Helvetica', 'Arial', sans-serif;
        margin: 0;
        padding: 0;
        box-sizing: border-box;
      }

      .container {
        display: flex;
        margin: 0;
        width: 100%;
      }

      .side-table {
        margin: 0;
        flex-shrink: 0;
        width: 60%;
        max-width: none;
      }

      .scrollable-table {
        height: 90%;
        overflow-y: auto;
        margin-top: 10px;
      }

      table tr:nth-child(even) {
        background-color: #f2f2f2;
      }

      #sortable-table tr td:last-child {
          display: flex;
          justify-content: center;
          align-items: center;
          height: 100%;
      }

      .viewer-container {
        display: flex;
        flex-direction: column;
        align-items: center;
        width: 40%;
      }

      #my-viewer {
        width: 100%;
      }

      .cgv-controls {
        width: 100%;
        display: flex;
        justify-content: center;
      }

      .generate-link {
          color: blue;
          text-decoration: underline;
          cursor: pointer;
      }

      .generated-link {
        color: blue;
        text-decoration: underline;
        cursor: pointer;
      }

      .btn-toggle-legend,
      #btn-toggle-legend {
        content: url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 30 30'%3E%3Ccircle cx='5' cy='7' r='2' fill='black'/%3E%3Cpath d='M10 7 h15' stroke='black' stroke-width='2'/%3E%3Ccircle cx='5' cy='15' r='2' fill='black'/%3E%3Cpath d='M10 15 h15' stroke='black' stroke-width='2'/%3E%3Ccircle cx='5' cy='23' r='2' fill='black'/%3E%3Cpath d='M10 23 h15' stroke='black' stroke-width='2'/%3E%3C/svg%3E");
      }

    </style>
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

        # Close the main content.
        file.write(
            """
    </div>
  </main>
        """
        )

        # Add script content.
        file.write(
            """
    <script src="./cgview-js_code/docs/scripts/bootstrap.min.js"></script>
    <script src="./cgview-js_code/docs/scripts/controls.js"></script>
    <script>
        function autoResizeMyViewer() {
            const mainPadding = 5;
            function myResize() {
                const myViewer = document.querySelector('#my-viewer');
                const main = document.getElementsByTagName('main')[0];
                const mainWidth = main.clientWidth * 0.50 - mainPadding;

                const mainHeight = main.clientHeight - 50;

                const minHeight = 500; // Set a minimum height if necessary
                const height = Math.max(mainHeight, minHeight);

                cgv.resize(mainWidth, height);
            }

            window.onresize = myResize;
            window.onload = function () {
                setTimeout( () => {
                    myResize();
                }, 100);
            }
        }

        function createViewer() {
            // Create Viewer in default div: #my-viewer
            const cgv = new CGV.Viewer('#my-viewer', {height: 500});
            autoResizeMyViewer();
            window.cgv = cgv;
        }

        function loadDataFile(filePath) {
            var script = document.createElement('script');
            script.src = filePath;
            script.onload = function() {
            cgv.io.loadJSON(window.json);
            cgv.draw();
            };
            document.head.appendChild(script);
        }

        document.addEventListener('DOMContentLoaded', (event) => {
            createViewer();

            // Add click event listener to each table row
            const rows = document.querySelectorAll('.side-table table tr');
            rows.forEach(row => {
            row.addEventListener('click', () => {
                const filePath = row.id;
                loadDataFile(filePath);
            });
            });

            // Check if there is at least one row and load the first genome by default
            if (rows.length > 0) {
            const firstRowFilePath = rows[1].id;
            loadDataFile(firstRowFilePath);
            }
        });

        function sortTable(column) {
            var table, rows, switching, i, x, y, shouldSwitch, dir, switchcount = 0;
            table = document.getElementById("sortable-table");
            switching = true;
            // Set the sorting direction to ascending:
            dir = "asc";
            // Make a loop that will continue until no switching has been done:
            while (switching) {
            // Start by saying: no switching is done:
            switching = false;
            rows = table.rows;
            // Loop through all table rows (except the first, which contains table headers):
            for (i = 1; i < (rows.length - 1); i++) {
                // Start by saying there should be no switching:
                shouldSwitch = false;
                // Get the two elements you want to compare, one from current row and one from the next:
                x = rows[i].getElementsByTagName("TD")[column];
                y = rows[i + 1].getElementsByTagName("TD")[column];
                // Check if the two rows should switch place, based on the direction, asc or desc:
                if (dir == "asc") {
                if (x.innerHTML.toLowerCase() > y.innerHTML.toLowerCase()) {
                    // If so, mark as a switch and break the loop:
                    shouldSwitch= true;
                    break;
                }
                } else if (dir == "desc") {
                if (x.innerHTML.toLowerCase() < y.innerHTML.toLowerCase()) {
                    // If so, mark as a switch and break the loop:
                    shouldSwitch = true;
                    break;
                }
                }
            }
            if (shouldSwitch) {
                // If a switch has been marked, make the switch and mark that a switch has been done:
                rows[i].parentNode.insertBefore(rows[i + 1], rows[i]);
                switching = true;
                // Each time a switch is done, increase this count by 1:
                switchcount ++;
            } else {
                // If no switching has been done AND the direction is "asc", set the direction to "desc" and run the while loop again.
                if (switchcount == 0 && dir == "asc") {
                dir = "desc";
                switching = true;
                }
            }
            }
        }

        // Code for generating Proksee projects.

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
                if (typeof window.json === 'undefined') {
                    console.error('No data found for sample ID:', sampleId);
                    alert('Failed to generate Proksee project: No data available for this sample.');
                    return;
                }

                const data = { origin: 'proksee-batch', data: JSON.stringify(window.json) };
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
                        link.textContent = 'Go to Proksee Project';
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

        function adjustTableHeight() {
            const headerHeight = document.querySelector('header').offsetHeight;
            const availableHeight = window.innerHeight - headerHeight;
            const scrollableTable = document.querySelector('.scrollable-table');
            scrollableTable.style.height = (availableHeight * 0.9) + 'px';
        }

        // Event listener for window resize
        window.addEventListener('resize', adjustTableHeight);

        // Initial adjustment on page load
        window.onload = function() {
            adjustTableHeight();
        };

        // Toggle Labels (not included in the controls.js file)
        onClick('btn-toggle-legend', () => {
        cgv.legend.update({visible: !cgv.legend.visible});
        cgv.draw();
        });

    </script>
</body>
</html>"""
        )
