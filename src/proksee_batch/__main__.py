"""Main module for the Proksee Batch tool."""
import json
import os
import shutil
import sys
from importlib import resources
from importlib.metadata import version
from typing import Any
from typing import Dict
from typing import Optional

import click
import toml

from .download_example_genbank_files import download_example_genbank_files
from .genbank_to_cgview_json import genbank_to_cgview_json
from .generate_proksee_link import generate_proksee_link
from .generate_report_html import generate_report_html
from .get_stats_from_genbank import get_stats_from_genbank
from .merge_cgview_json_with_template import merge_cgview_json_with_template
from .parse_additional_features import add_bed_features_and_tracks
from .parse_additional_features import add_blast_features_and_tracks
from .scrape_proksee_image import scrape_proksee_image


@click.command()
@click.version_option(version=version("proksee-batch"), prog_name="Proksee Batch")
@click.option(
    "--genomes",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help="Path to the directory containing input genomes in GenBank format.",
)
@click.option(
    "--blast-results",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help="Path to a directory containing BLAST results in tabular (outfmt 6) format.",
)
@click.option(
    "--bed-features",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help="Path to a directory containing BED files with additional genomic features.",
)
@click.option(
    "--output",
    type=click.Path(file_okay=False, dir_okay=True),
    help="Path to the output directory.",
)
@click.option(
    "--template",
    type=click.Path(exists=True),
    help="Path to the Proksee configuration file in JSON format.",
)
@click.option(
    "--download-example",
    type=click.Path(file_okay=False, dir_okay=True),
    help="Download example GenBank files to the specified directory, and exit.",
)
def main(
    genomes: Optional[str],
    blast_results: Optional[str],
    bed_features: Optional[str],
    output: Optional[str],
    template: Optional[str],
    download_example: Optional[str],
) -> None:
    """Proksee Batch: A tool for visualizing batches of genomes via https://www.proksee.ca."""

    # Check if the download example GenBank files option is used
    if download_example:
        # Download logic
        download_example_genbank_files(download_example)
        print(f"Example GenBank files downloaded to {download_example}")
        sys.exit(0)

    # Check if genomes directory is provided
    if not genomes:
        print("Error: Missing option '--genomes'.", file=sys.stderr)
        sys.exit(1)

    # Check if output directory is provided
    if not output:
        print("Error: Missing option '--output'.", file=sys.stderr)
        sys.exit(1)

    # Process genomes directory.
    # Check that there is at least one .gbk file in the directory.
    if not any(
        genome.endswith(".gbk") or genome.endswith(".gbff") or genome.endswith(".gb")
        for genome in os.listdir(genomes)
    ):
        print("No GenBank files found in the genomes directory.", file=sys.stderr)
        sys.exit(1)

    # Process BLAST results directory (if provided).
    if blast_results:
        # Check that there is at least one .txt or .tsv file in the directory.
        if not any(
            blast_result.endswith(".txt") or blast_result.endswith(".tsv")
            for blast_result in os.listdir(blast_results)
        ):
            print(
                "No BLAST result files (.txt or .tsv) found in the BLAST results directory.",
                file=sys.stderr,
            )
            sys.exit(1)
        # Check that the basename of each BLAST result file matches the basename of a GenBank file.
        blast_result_genome_ids = [
            blast_result.rsplit(".", 1)[0]
            for blast_result in os.listdir(blast_results)
            if blast_result.endswith(".txt") or blast_result.endswith(".tsv")
        ]
        genome_ids = [genome.rsplit(".", 1)[0] for genome in os.listdir(genomes)]
        for blast_result_genome_id in blast_result_genome_ids:
            corresponding_genomes = []
            for genome_id in genome_ids:
                if blast_result_genome_id.startswith(genome_id):
                    corresponding_genomes.append(genome_id)
            if len(corresponding_genomes) == 0:
                print(
                    f"No GenBank file found for {blast_result_genome_id}. BLAST result file names must start with the genome name.",
                    file=sys.stderr,
                )
                sys.exit(1)
            elif len(corresponding_genomes) > 1:
                print(
                    f"BLAST result file {blast_result_genome_id} is ambiguously assigned to multiple genomes: {corresponding_genomes}. BLAST result file names must start with a unique genome name.",
                    file=sys.stderr,
                )
                sys.exit(1)

    # Process BED features directory (if provided).
    if bed_features:
        # Check that there is at least one .bed file in the directory.
        if not any(
            bed_feature.endswith(".bed") for bed_feature in os.listdir(bed_features)
        ):
            print("No BED files found in the BED features directory.", file=sys.stderr)
            sys.exit(1)
        # Check that the basename of each BED file matches the basename of a GenBank file.
        bed_feature_genome_ids = [
            bed_feature.rsplit(".", 1)[0]
            for bed_feature in os.listdir(bed_features)
            if bed_feature.endswith(".bed")
        ]
        genome_ids = [genome.rsplit(".", 1)[0] for genome in os.listdir(genomes)]
        for bed_feature_genome_id in bed_feature_genome_ids:
            corresponding_genomes = []
            for genome_id in genome_ids:
                if bed_feature_genome_id.startswith(genome_id):
                    corresponding_genomes.append(genome_id)
            if len(corresponding_genomes) == 0:
                print(
                    f"No GenBank file found for {bed_feature_genome_id}. BED file names must start with the genome name.",
                    file=sys.stderr,
                )
                sys.exit(1)
            elif len(corresponding_genomes) > 1:
                print(
                    f"BED file {bed_feature_genome_id} is ambiguously assigned to multiple genomes: {corresponding_genomes}. BED file names must start with a unique genome name.",
                    file=sys.stderr,
                )
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

    else:
        with resources.path(
            "proksee_batch.data", "default_proksee_template.json"
        ) as template_path:
            template = str(template_path)
        assert os.path.exists(template)

    # Define path to a temporary output directory within the output directory.
    temp_output = os.path.join(output, "temp")
    os.mkdir(temp_output)

    # Initiate a dictionary to store file paths for each genome.
    genome_files: Dict[str, Any] = {}

    # Iterate over the .gbk files in the genomes directory, and process each one to create a .json file in the output directory.
    for genome in os.listdir(genomes):
        if (
            genome.endswith(".gbk")
            or genome.endswith(".gbff")
            or genome.endswith(".gb")
        ):
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
            basic_json_file = os.path.join(
                temp_output, genome.rsplit(".", 1)[0] + ".json"
            )
            genbank_to_cgview_json(os.path.join(genomes, genome), basic_json_file)

            # Find any BLAST result files for this genome.
            blast_files = []
            if blast_results:
                for blast_result in os.listdir(blast_results):
                    if blast_result.endswith(".txt") or blast_result.endswith(".tsv"):
                        if blast_result.startswith(genome_id):
                            blast_files.append(
                                os.path.join(blast_results, blast_result)
                            )

            # Parse the BLAST result files (if any) to create additional features and tracks.
            basic_json_file_with_blast_features = os.path.join(
                temp_output, genome.rsplit(".", 1)[0] + ".with_blast_features.json"
            )
            if blast_files:
                add_blast_features_and_tracks(
                    blast_files, basic_json_file, basic_json_file_with_blast_features
                )
            else:
                shutil.copy(basic_json_file, basic_json_file_with_blast_features)

            # Find any BED files for this genome.
            bed_files = []
            if bed_features:
                for bed_feature in os.listdir(bed_features):
                    if bed_feature.endswith(".bed"):
                        if bed_feature.startswith(genome_id):
                            bed_files.append(os.path.join(bed_features, bed_feature))

            # Parse the BED files (if any) to create additional features and tracks.
            basic_json_file_with_bed_features = os.path.join(
                temp_output, genome.rsplit(".", 1)[0] + ".with_bed_features.json"
            )
            if bed_files:
                add_bed_features_and_tracks(
                    bed_files,
                    basic_json_file_with_blast_features,
                    basic_json_file_with_bed_features,
                )
            else:
                shutil.copy(
                    basic_json_file_with_blast_features,
                    basic_json_file_with_bed_features,
                )

            # Merge the basic cgview map with the template Proksee configuration file.
            merged_json_file = os.path.join(
                temp_output, genome.rsplit(".", 1)[0] + ".merged.json"
            )
            merge_cgview_json_with_template(
                basic_json_file_with_bed_features, template, merged_json_file
            )

            # Convert the merged JSON file to .js file by wrapping it in a variable assignment.
            js_file = os.path.join(temp_output, genome.rsplit(".", 1)[0] + ".js")
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
                temp_output, genome.rsplit(".", 1)[0] + ".proksee_link.txt"
            )
            generate_proksee_link(merged_json_file, proksee_project_link_file)

            # Scrape proksee image.
            proksee_image_file = os.path.join(
                temp_output, genome.rsplit(".", 1)[0] + ".svg"
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
