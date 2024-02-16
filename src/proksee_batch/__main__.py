"""Main module for the Proksee Batch tool.

This file defines the proksee-batch command-line interface and implements the
high-level logic of the tool.

"""
import datetime
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

from .download_example_input import download_example_input

# from .generate_proksee_link import generate_proksee_link
from .generate_report_html import generate_report_html
from .get_stats_from_seq_file import get_stats_from_seq_file
from .merge_cgview_json_with_template import merge_cgview_json_with_template
from .parse_additional_features import add_bed_features_and_tracks
from .parse_additional_features import add_blast_features_and_tracks
from .parse_additional_features import add_gff_features_and_tracks
from .parse_additional_features import add_vcf_features_and_tracks
from .seq_file_to_cgview_json import fasta_to_cgview_json
from .seq_file_to_cgview_json import genbank_to_cgview_json

# from .scrape_proksee_image import scrape_proksee_image
from .validate_input_data import get_data_files
from .validate_input_data import validate_input_directory_contents


@click.command()
@click.version_option()
@click.option(
    "--input",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help="Path to the directory containing input genomes in GenBank format.",
)
@click.option(
    "--output",
    type=click.Path(file_okay=False, dir_okay=True),
    help="Path to the output directory.",
)
@click.option(
    "--download-example-data",
    type=click.Path(file_okay=False, dir_okay=True),
    help="Download example input files to the specified directory, and exit.",
)
def main(
    input: Optional[str],
    output: Optional[str],
    download_example_data: Optional[str],
) -> None:
    """
    Proksee Batch: A tool for visualizing batches of genomes via https://www.proksee.ca.
    """
    # Get current date and time.
    run_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Check if the download example GenBank files option is used
    if download_example_data:
        # Download logic
        download_example_input(download_example_data)
        print(f"Example data files downloaded to {download_example_data}")
        sys.exit(0)

    # Check if the input directory is provided.
    input_dir_path = ""
    if not input:
        handle_error_exit("Missing option '--input'.")
    else:
        assert input is not None
        if not os.path.exists(input):
            handle_error_exit(f"The input directory does not exist: {input}")
        input_dir_path = os.path.abspath(input)

    # Validate the contents of the input directory.
    validate_input_directory_contents(input_dir_path)

    # Make the output directory if it doesn't exist.
    output_path = ""
    if not output:
        handle_error_exit("Missing option '--output'.")
    else:
        output_path = os.path.abspath(output)
        if os.path.exists(output_path):
            shutil.rmtree(output_path)
        os.mkdir(output_path)
        os.mkdir(os.path.join(output_path, "data"))
        # os.mkdir(os.path.join(output_path, "images"))

    # Initiate a dictionary to store file paths for each genome.
    genome_info: Dict[str, Dict[str, Any]] = {}

    # Iterate over subdirectories in the input directory.
    genome_num = 0
    for genome_dir in os.listdir(input):
        # Skip any files in the input directory.
        if os.path.isfile(os.path.join(input_dir_path, genome_dir)):
            continue

        genome_num += 1
        genome_code_name = "genome_" + str(genome_num)

        # Define path to a temporary output directory within the output directory.
        temp_output = os.path.join(output_path, genome_code_name + "_temp")
        os.mkdir(temp_output)

        # Get the data file paths
        genbank_paths = get_data_files(
            os.path.join(input_dir_path, genome_dir), "genbank"
        )
        fasta_paths = get_data_files(os.path.join(input_dir_path, genome_dir), "fasta")
        json_paths = get_data_files(os.path.join(input_dir_path, genome_dir), "json")
        blast_paths = get_data_files(os.path.join(input_dir_path, genome_dir), "blast")
        bed_paths = get_data_files(os.path.join(input_dir_path, genome_dir), "bed")
        vcf_paths = get_data_files(os.path.join(input_dir_path, genome_dir), "vcf")
        gff_paths = get_data_files(os.path.join(input_dir_path, genome_dir), "gff")

        # Generate a dict with file types as keys and lists of file paths as values.
        file_names_dict = {
            "genbank": [os.path.basename(x) for x in genbank_paths],
            "fasta": [os.path.basename(x) for x in fasta_paths],
            "json": [os.path.basename(x) for x in json_paths],
            "blast": [os.path.basename(x) for x in blast_paths],
            "bed": [os.path.basename(x) for x in bed_paths],
            "vcf": [os.path.basename(x) for x in vcf_paths],
            "gff": [os.path.basename(x) for x in gff_paths],
        }

        # If the genome directory contains one or more GenBank files, remove the
        # list of FASTA files from the dict (because they won't be used).
        if genbank_paths:
            del file_names_dict["fasta"]

        # Define path to the basic cgview map in JSON format.
        basic_json_file = os.path.join(temp_output, genome_code_name + ".json")

        if genbank_paths:
            genbank_path = genbank_paths[0]

            # Get basic stats from the GenBank file.
            (
                genbank_description,
                genbank_total_size,
                genbank_number_of_contigs,
                genbank_gc_content,
            ) = get_stats_from_seq_file(genbank_path, "genbank")

            genome_info[genome_code_name] = {
                "Name": genome_dir,
                "Description": genbank_description,
                "Total size": genbank_total_size,
                "Number of contigs": genbank_number_of_contigs,
                "GC content": genbank_gc_content,
                "Files": file_names_dict,
            }

            # Convert the GenBank file to a basic cgview map in JSON format.
            genbank_to_cgview_json(genome_dir, genbank_path, basic_json_file)

        elif fasta_paths:
            fasta_path = fasta_paths[0]

            # Get basic stats from the FASTA file.
            (
                fasta_description,
                fasta_total_size,
                fasta_number_of_contigs,
                fasta_gc_content,
            ) = get_stats_from_seq_file(fasta_path, "fasta")

            genome_info[genome_code_name] = {
                "Name": genome_dir,
                "Description": fasta_description,
                "Total size": fasta_total_size,
                "Number of contigs": fasta_number_of_contigs,
                "GC content": fasta_gc_content,
                "Files": file_names_dict,
            }

            # Convert the FASTA file to a basic cgview map in JSON format.
            fasta_to_cgview_json(genome_dir, fasta_path, basic_json_file)

        else:
            handle_error_exit(
                f"Could not find a GenBank or FASTA file in the {genome_dir} directory."
            )

        # Parse the BLAST result files (if any) to create additional features and tracks.
        basic_json_file_with_blast_features = os.path.join(
            temp_output, genome_code_name + ".with_blast_features.json"
        )
        if blast_paths:
            add_blast_features_and_tracks(
                blast_paths, basic_json_file, basic_json_file_with_blast_features
            )
        else:
            shutil.copy(basic_json_file, basic_json_file_with_blast_features)

        # Parse the BED files (if any) to create additional features and tracks.
        basic_json_file_with_bed_features = os.path.join(
            temp_output, genome_code_name + ".with_bed_features.json"
        )
        if bed_paths:
            add_bed_features_and_tracks(
                bed_paths,
                basic_json_file_with_blast_features,
                basic_json_file_with_bed_features,
            )
        else:
            shutil.copy(
                basic_json_file_with_blast_features,
                basic_json_file_with_bed_features,
            )

        # Parse the VCF files (if any) to create additional features and tracks.
        basic_json_file_with_vcf_features = os.path.join(
            temp_output, genome_code_name + ".with_vcf_features.json"
        )
        if vcf_paths:
            add_vcf_features_and_tracks(
                vcf_paths,
                basic_json_file_with_bed_features,
                basic_json_file_with_vcf_features,
            )
        else:
            shutil.copy(
                basic_json_file_with_bed_features,
                basic_json_file_with_vcf_features,
            )

        # Parse the GFF files (if any) to create additional features and tracks.
        basic_json_file_with_gff_features = os.path.join(
            temp_output, genome_code_name + ".with_gff_features.json"
        )
        if gff_paths:
            add_gff_features_and_tracks(
                gff_paths,
                basic_json_file_with_vcf_features,
                basic_json_file_with_gff_features,
            )
        else:
            shutil.copy(
                basic_json_file_with_vcf_features,
                basic_json_file_with_gff_features,
            )

        # Determine path to the template Proksee configuration file.
        template_path = ""
        if json_paths:
            try:
                with open(json_paths[0]) as template_file:
                    # Load the template file as a JSON object.
                    configuration = json.load(template_file)
                    template_path = json_paths[0]
            except json.JSONDecodeError as e:
                print(f"Error reading the template file: {e}", file=sys.stderr)
                sys.exit(1)
        else:
            with resources.path(
                "proksee_batch.data", "default_proksee_template.json"
            ) as template_json_path:
                template_path = str(template_json_path)
        assert os.path.exists(template_path)

        # Merge the basic cgview map with the template Proksee configuration file.
        merged_json_file = os.path.join(temp_output, genome_code_name + ".merged.json")
        merge_cgview_json_with_template(
            basic_json_file_with_gff_features, template_path, merged_json_file
        )

        # Add GC content and GC skew tracks to the merged JSON.
        merged_json_file_with_gc_tracks = os.path.join(
            temp_output, genome_code_name + ".merged_with_gc_tracks.json"
        )
        with open(merged_json_file) as merged_json_fh, open(
            merged_json_file_with_gc_tracks, "w"
        ) as merged_json_file_with_gc_tracks_fh:
            merged_json = json.load(merged_json_fh)
            merged_json["cgview"]["tracks"].extend(
                [
                    {
                        "name": "CG Content",
                        "thicknessRatio": 2,
                        "position": "inside",
                        "dataType": "plot",
                        "dataMethod": "sequence",
                        "dataKeys": "gc-content",
                    },
                    {
                        "name": "CG Skew",
                        "thicknessRatio": 2,
                        "position": "inside",
                        "dataType": "plot",
                        "dataMethod": "sequence",
                        "dataKeys": "gc-skew",
                    },
                ]
            )
            json.dump(merged_json, merged_json_file_with_gc_tracks_fh)

        # Convert the merged JSON file to .js file by wrapping it in a variable assignment.
        js_file = os.path.join(output_path, "data", genome_code_name + ".js")
        with open(js_file, "w") as file:
            # Get the JSON data from the merged JSON file as a string.
            json_data = None
            with open(merged_json_file_with_gc_tracks) as json_file:
                json_data = json_file.read()
            # Write the variable assignment.
            file.write(f"json = {json_data};")

        # Generate a .js file with the contents of the genome_info dictionary.
        generate_js_data(output_path, genome_info, run_date, input_dir_path)

        ## Generate a Proksee link using the merged JSON file.
        # proksee_project_link_file = os.path.join(
        #    temp_output, genome_code_name + ".proksee_link.txt"
        # )
        # generate_proksee_link(merged_json_file, proksee_project_link_file)

        ## Scrape proksee image.
        # proksee_image_file = os.path.join(
        #    os.path.join(output_path, "images", genome_code_name) + ".svg"
        # )
        # scrape_proksee_image(proksee_project_link_file, proksee_image_file)

        # Delete the temporary output directory.
        shutil.rmtree(temp_output)

    # Copy the directory with CGView.js code from the package data to the output directory.
    with resources.path("proksee_batch.data", "cgview-js_code") as cgview_js_path:
        shutil.copytree(cgview_js_path, os.path.join(output_path, "cgview-js_code"))

    # Copy the directory with supporting CSS and JavaScript code for the HTML report from the package data to the output directory.
    with resources.path(
        "proksee_batch.data", "html_report_code"
    ) as report_support_files_path:
        shutil.copytree(
            report_support_files_path, os.path.join(output_path, "html_report_code")
        )

    ## Generate the HTML report file.
    # generate_report_html(output_path, genome_info)

    # Copy the HTML file with the report from the package data to the output directory.
    with resources.path("proksee_batch.data", "report.html") as report_html_path:
        shutil.copy(report_html_path, os.path.join(output_path, "report.html"))


def generate_js_data(
    output_dir: str, genome_info: Dict[str, Any], run_date: str, input_dir: str
) -> None:
    """
    Generates a JavaScript file with genome information wrapped in a variable assignment.
    """
    assert os.path.isdir(output_dir), f"Output directory does not exist: {output_dir}"
    output_file = os.path.join(output_dir, "data", "table_data.js")

    # Construct a dict to contain all the info.
    all_info = {
        "run_date": run_date,
        "input_dir": input_dir,
        "genomes": genome_info,
    }

    # Prepare the data to be wrapped in a JavaScript variable
    js_content = "const tableData = " + json.dumps(all_info, indent=4) + ";"

    # Save the data to a .js file
    with open(output_file, "w") as file:
        file.write(js_content)


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


if __name__ == "__main__":
    main()
