"""Main module for the Proksee Batch tool."""
import json
import os
import shutil
import sys
from typing import Any
from typing import Dict
from typing import Optional

import click

from .genbank_to_cgview_json import genbank_to_cgview_json
from .generate_proksee_link import generate_proksee_link
from .generate_report_html import generate_report_html
from .get_stats_from_genbank import get_stats_from_genbank
from .merge_cgview_json_with_template import merge_cgview_json_with_template
from .scrape_proksee_image import scrape_proksee_image


@click.command()
@click.option(
    "--genomes",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help="Path to the directory containing input genomes in GenBank format.",
)
@click.option(
    "--output",
    required=True,
    type=click.Path(file_okay=False, dir_okay=True),
    help="Path to the output directory.",
)
@click.option(
    "--template",
    type=click.Path(exists=True),
    help="Path to the Proksee configuration file in JSON format.",
)
def main(genomes: str, output: str, template: Optional[str]) -> None:
    """Proksee Batch: A tool for batch processing of genomes using Proksee
    (https://proksee.ca/)."""
    # Process genomes directory.
    # Check that there is at least one .gbk file in the directory.
    if not any(genome.endswith(".gbk") for genome in os.listdir(genomes)):
        print("No GenBank files found in the genomes directory.", file=sys.stderr)
        sys.exit(1)

    # Process output directory
    # If it exists, delete it and create a new one.
    if os.path.exists(output):
        shutil.rmtree(output)
    os.mkdir(output)

    # Process template file
    if template:
        try:
            with open(template) as template_file:
                configuration = json.load(template_file)
                # Process the configuration as needed
                print("Configuration loaded successfully.")
        except json.JSONDecodeError as e:
            print(f"Error reading the template file: {e}", file=sys.stderr)
            sys.exit(1)

    # Set template to the path to a default Proksee configuration file if no template was provided.
    else:
        template = os.path.join(
            os.path.dirname(__file__), "data", "default_proksee_template.json"
        )

    # Define path to a temporary output directory within the output directory.
    temp_output = os.path.join(output, "temp")
    os.mkdir(temp_output)

    # Initiate a dictionary to store file paths for each genome.
    genome_files: Dict[str, Any] = {}

    # Iterate over the .gbk files in the genomes directory, and process each one to create a .json file in the output directory.
    for genome in os.listdir(genomes):
        if genome.endswith(".gbk"):
            print(f"Processing {genome}...")

            # Get basic stats from the GenBank file.
            (
                genbank_description,
                genbank_total_size,
                genbank_number_of_contigs,
                genbank_gc_content,
            ) = get_stats_from_genbank(os.path.join(genomes, genome))

            genome_id = genome.rsplit(".", 1)[0]
            genome_files[genome_id] = {
                "JSON_file": "",
                "Proksee_image_file": "",
                "Description": genbank_description,
                "Total size": genbank_total_size,
                "Number of contigs": genbank_number_of_contigs,
                "GC content": genbank_gc_content,
            }

            # Convert the GenBank file to a basic cgview map in JSON format.
            basic_json_file = os.path.join(temp_output, genome.replace(".gbk", ".json"))
            genbank_to_cgview_json(os.path.join(genomes, genome), basic_json_file)

            # Merge the basic cgview map with the template Proksee configuration file.
            merged_json_file = os.path.join(
                temp_output, genome.replace(".gbk", ".merged.json")
            )
            merge_cgview_json_with_template(basic_json_file, template, merged_json_file)

            # Convert the merged JSON file to .js file by wrapping it in a variable assignment.
            js_file = os.path.join(temp_output, genome.replace(".gbk", ".js"))
            genome_files[genome_id]["JSON_file"] = js_file
            with open(js_file, "w") as file:
                # Get the JSON data from the merged JSON file as a string.
                json_data = None
                with open(merged_json_file) as json_file:
                    json_data = json_file.read()
                # Write the variable assignment.
                file.write(f"var jsonData = {json_data};")

            # Generate a Proksee link using the merged JSON file.
            proksee_project_link_file = os.path.join(
                temp_output, genome.replace(".gbk", ".proksee_link.txt")
            )
            generate_proksee_link(merged_json_file, proksee_project_link_file)

            # Scrape proksee image.
            proksee_image_file = os.path.join(
                temp_output, genome.replace(".gbk", ".svg")
            )
            genome_files[genome_id]["Proksee_image_file"] = proksee_image_file
            scrape_proksee_image(proksee_project_link_file, proksee_image_file)

    # Make "data" and "images" subdirectories in the report directory.
    data_dir = os.path.join(output, "data")
    os.mkdir(data_dir)
    images_dir = os.path.join(output, "images")
    os.mkdir(images_dir)

    # Copy the .js files to the data directory, and the .svg files to the images
    # directory, and generate lists of both sets of file paths.
    js_files = []
    image_files = []
    for genome_id in genome_files:
        js_files.append(genome_files[genome_id]["JSON_file"])
        image_files.append(genome_files[genome_id]["Proksee_image_file"])
        shutil.copy(genome_files[genome_id]["JSON_file"], data_dir)
        shutil.copy(genome_files[genome_id]["Proksee_image_file"], images_dir)

    # Generate the HTML report file.
    report_file = os.path.join(output, "report.html")
    generate_report_html(js_files, image_files, genome_files, report_file)

    # Delete the temporary output directory.
    shutil.rmtree(temp_output)


if __name__ == "__main__":
    main()
