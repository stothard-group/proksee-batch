"""Main module for the Proksee Batch tool."""
import json
import os
import shutil
import sys
from importlib import resources
from importlib.metadata import version
from typing import Any
from typing import Dict
from typing import List
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

    # Validate contents of the genomes directory.
    validate_directory_contents(
        genomes, [".gbk", ".gbff", ".gb"], "GenBank files", "--genomes"
    )
    genomes_path = ""
    if genomes:
        genomes_path = os.path.abspath(genomes)

    # Get a list of IDs and a list of filenames for the genomes.
    genome_filenames = [
        genome
        for genome in os.listdir(genomes_path)
        if genome.endswith(".gbk") or genome.endswith(".gbff") or genome.endswith(".gb")
    ]
    genome_ids = [genome.rsplit(".", 1)[0] for genome in genome_filenames]

    output_path = ""
    if not output:
        handle_error_exit("Missing option '--output'.")
    else:
        output_path = os.path.abspath(output)
        if os.path.exists(output_path):
            shutil.rmtree(output_path)
        os.mkdir(output_path)

    if blast_results:
        validate_directory_contents(
            blast_results, [".txt", ".tsv"], "BLAST result files", "--blast-results"
        )
        check_filenames_in_directory(
            blast_results, [".txt", ".tsv"], genome_ids, "BLAST result"
        )
    if bed_features:
        validate_directory_contents(
            bed_features, [".bed"], "BED files", "--bed-features"
        )
        check_filenames_in_directory(bed_features, [".bed"], genome_ids, "BED")

    # Process template file
    if template:
        try:
            with open(template) as template_file:
                # Load the template file as a JSON object.
                configuration = json.load(template_file)
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
    temp_output = os.path.join(output_path, "temp")
    os.mkdir(temp_output)

    # Initiate a dictionary to store file paths for each genome.
    genome_files: Dict[str, Any] = {}

    # Iterate over the .gbk files in the genomes directory, and process each one to create a .json file in the output directory.
    for genome in genome_filenames:
        print(f"Processing {genome}...")

        # Get basic stats from the GenBank file.
        (
            genbank_description,
            genbank_total_size,
            genbank_number_of_contigs,
            genbank_gc_content,
        ) = get_stats_from_genbank(os.path.join(genomes_path, genome))

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
        basic_json_file = os.path.join(temp_output, genome.rsplit(".", 1)[0] + ".json")
        genbank_to_cgview_json(os.path.join(genomes_path, genome), basic_json_file)

        # Find any BLAST result files for this genome.
        blast_files = []
        if blast_results:
            for blast_result in os.listdir(blast_results):
                if blast_result.endswith(".txt") or blast_result.endswith(".tsv"):
                    if blast_result.startswith(genome_id):
                        blast_files.append(os.path.join(blast_results, blast_result))

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
    data_dir = os.path.join(output_path, "data")
    os.mkdir(data_dir)
    images_dir = os.path.join(output_path, "images")
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
    report_file = os.path.join(output_path, "report.html")
    generate_report_html(js_files, image_files, genome_files, report_file)

    # Delete the temporary output directory.
    shutil.rmtree(temp_output)


def validate_directory_contents(
    directory: Optional[str],
    extensions: List[str],
    file_description: str,
    option_name: str,
) -> None:
    """
    Validates if the provided directory contains files with the specified extensions.

    Args:
        directory (Optional[str]): The path to the directory to be checked.
        extensions (List[str]): A list of acceptable file extensions.
        file_description (str): A description of the expected files, for error messages.
        option_name (str): The name of the command-line option corresponding to the directory.

    Raises:
        SystemExit: If the directory is not provided, doesn't exist, or doesn't contain files with the required extensions.
    """
    if not directory:
        handle_error_exit(f"Missing option '{option_name}'.")

    else:
        if not os.path.exists(directory):
            handle_error_exit(
                f"The directory specified in '{option_name}' does not exist: {directory}"
            )

        if not any(file.endswith(tuple(extensions)) for file in os.listdir(directory)):
            handle_error_exit(
                f"No {file_description} found in the directory specified in '{option_name}': {directory}"
            )


def handle_error_exit(error_message: str, exit_code: int = 1) -> None:
    """
    Handles errors by printing a message to sys.stderr and exiting the program.

    Args:
        error_message (str): The error message to be printed.
        exit_code (int): The exit code to be used for sys.exit. Defaults to 1.

    Exits:
        SystemExit: Exits the program with the provided exit code.
    """
    print(error_message, file=sys.stderr)
    sys.exit(exit_code)


def check_filenames_in_directory(
    directory: str, file_extensions: List[str], genome_ids: List[str], file_type: str
) -> None:
    """
    Checks that the basename of each file in the directory matches the basename of a GenBank file.

    Args:
        directory (str): The path to the directory to be checked.
        file_extensions (List[str]): A list of acceptable file extensions.
        genome_ids (List[str]): A list of genome IDs.
        file_type (str): A description of the type of file being checked, for error messages.

    Raises:
        SystemExit: If a file in the directory does not have a basename that matches the basename of a GenBank file.
    """
    # Check that the basename of each file matches the basename of a GenBank file.
    file_genome_ids = [
        file.rsplit(".", 1)[0]
        for file in os.listdir(directory)
        if any(file.endswith(ext) for ext in file_extensions)
    ]
    for file_genome_id in file_genome_ids:
        corresponding_genomes = []
        for genome_id in genome_ids:
            if file_genome_id.startswith(genome_id):
                corresponding_genomes.append(genome_id)
        if len(corresponding_genomes) == 0:
            print(
                f"No GenBank file found for {file_genome_id}. {file_type} file names must start with the genome name.",
                file=sys.stderr,
            )
            sys.exit(1)
        elif len(corresponding_genomes) > 1:
            print(
                f"{file_type} file {file_genome_id} is ambiguously assigned to multiple genomes: {corresponding_genomes}. {file_type} file names must start with a unique genome name.",
                file=sys.stderr,
            )
            sys.exit(1)


if __name__ == "__main__":
    main()
