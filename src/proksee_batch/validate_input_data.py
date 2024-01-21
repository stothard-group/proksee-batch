"""

The input directory must be structured as in the following example:

    input_directory/
        genome_name_1/
            genbank/
                genome1.gbk
            blast/
                abc.txt
                def.tsv
            bed/
                ghi.bed
                jkl.bed
            json/
                template1.json
            vcf/
                mno.vcf
                pqr.vcf
            gff/
                stu.gff
                vwx.gff3
        genome_name_2/
            genbank/
                genome2.gbff
            blast/
                yza.txt
                bcd.tsv
            bed/
                efg.bed
                hij.bed
            json/
                template2.json
            vcf/
                klm.vcf
                nop.vcf
            gff/
                qrs.gff
                tuv.gff3
        ...

The genbank directory must a single GenBank file with the extension .gbk, .gbff,
or .gb. This is the genome that will be visualized.  The blast, bed, vcf, and
gff directories are optional. They contain files with additional genomic
features.  The json directory is also optional. It contains a custom Proksee
project JSON file that will be used as a template for the visualization.
"""
import json
import os
import sys
from typing import List

from Bio import SeqIO


def validate_input_directory_contents(input: str) -> None:
    """
    Validates if the provided input directory contains the required subdirectories and files.

    Args:
        input (str): The path to the input directory to be checked.

    Raises:
        SystemExit: If the input directory does not contain the required subdirectories.
    """
    # Iterate over subdirectories in the input directory (each subdirectory
    # should contain the data for one genome).
    for genome_dir in os.listdir(input):
        required_subdirectories = ["genbank"]
        for subdirectory in required_subdirectories:
            if not os.path.isdir(os.path.join(input, genome_dir, subdirectory)):
                handle_error_exit(
                    f"The input directory does not contain the required subdirectory: {subdirectory}"
                )

        # Check that there is exactly one GenBank file in the genbank directory.
        genbank_files = [
            file
            for file in os.listdir(os.path.join(input, genome_dir, "genbank"))
            if file.endswith(".gbk") or file.endswith(".gbff") or file.endswith(".gb")
        ]
        if len(genbank_files) == 0:
            handle_error_exit(
                "The input directory does not contain a GenBank file in the 'genbank' subdirectory."
            )
        elif len(genbank_files) > 1:
            handle_error_exit(
                "The input directory contains more than one GenBank file in the 'genbank' subdirectory."
            )

        # Parse the GenBank file with BioPython to check if the file format is valid.
        try:
            SeqIO.parse(
                os.path.join(input, genome_dir, "genbank", genbank_files[0]), "genbank"
            ).__next__()  # type: ignore
        except Exception as e:
            handle_error_exit(
                f"Error parsing the GenBank file: {e}. Please check that the file is in GenBank format."
            )

        # Check that all other subdirectories have valid names.
        optional_subdirectories = ["blast", "bed", "json", "vcf", "gff"]
        for subdirectory in os.listdir(os.path.join(input, genome_dir)):
            if (
                subdirectory not in required_subdirectories
                and subdirectory not in optional_subdirectories
            ):
                handle_error_exit(
                    f"The input directory contains an invalid subdirectory: {subdirectory}"
                )

        # Check the contents of the optional subdirectories, if they exist.
        for subdirectory in optional_subdirectories:
            subdirectory_path = os.path.join(input, genome_dir, subdirectory)
            if os.path.isdir(subdirectory_path):
                if subdirectory == "blast":
                    blast_result_files = []
                    for file in os.listdir(subdirectory_path):
                        if file.endswith(".txt") or file.endswith(".tsv"):
                            blast_result_files.append(file)
                            # Validate the BLAST result file. It should be in tabular format (outfmt 6).
                            with open(os.path.join(subdirectory_path, file)) as f:
                                for line in f:
                                    if not line.startswith("#"):
                                        line_split = line.split("\t")
                                        if len(line_split) != 12:
                                            handle_error_exit(
                                                f"The BLAST result file {file} does not have the correct number of columns. Please check that the file is in tabular format (outfmt 6)."
                                            )
                                        try:
                                            float(line_split[2])
                                            float(line_split[3])
                                            float(line_split[10])
                                            float(line_split[11])
                                        except ValueError:
                                            handle_error_exit(
                                                f"The BLAST result file {file} does not have the correct format. Please check that the file is in tabular format (outfmt 6)."
                                            )
                    if len(blast_result_files) == 0:
                        handle_error_exit(
                            f"The input directory {subdirectory_path} does not contain any BLAST result files. Valid file extensions are .txt and .tsv."
                        )
                elif subdirectory == "bed":
                    bed_files = []
                    for file in os.listdir(subdirectory_path):
                        if file.endswith(".bed"):
                            bed_files.append(file)
                            # Validate the BED file. It should be in BED format.
                            with open(os.path.join(subdirectory_path, file)) as f:
                                for line in f:
                                    line_split = line.split("\t")
                                    if len(line_split) < 3:
                                        handle_error_exit(
                                            f"The BED file {file} does not have the correct number of columns. Please check that the file is in BED format."
                                        )
                                    try:
                                        int(line_split[1])
                                        int(line_split[2])
                                    except ValueError:
                                        handle_error_exit(
                                            f"The BED file {file} does not have the correct format. Please check that the file is in BED format."
                                        )
                    if len(bed_files) == 0:
                        handle_error_exit(
                            f"The input directory {subdirectory_path} does not contain any BED files. Valid file extension is .bed."
                        )
                elif subdirectory == "json":
                    json_files = []
                    for file in os.listdir(subdirectory_path):
                        if file.endswith(".json"):
                            json_files.append(file)
                            # Validate the JSON file. It should be in JSON format and contain the required keys.
                            with open(os.path.join(subdirectory_path, file)) as f:
                                try:
                                    json_dict = json.load(f)
                                except Exception as e:
                                    handle_error_exit(
                                        f"Error parsing the JSON file {file}: {e}. Please check that the file is in JSON format."
                                    )
                                required_keys = [
                                    "cgview",
                                ]
                                for key in required_keys:
                                    if key not in json_dict:
                                        handle_error_exit(
                                            f"The JSON file {file} does not contain the required key: {key}"
                                        )
                    if len(json_files) == 0:
                        handle_error_exit(
                            f"The input directory {subdirectory_path} does not contain any JSON files. Valid file extension is .json."
                        )
                elif subdirectory == "vcf":
                    vcf_files = []
                    for file in os.listdir(subdirectory_path):
                        if file.endswith(".vcf"):
                            vcf_files.append(file)
                            # Validate the VCF file. It should be in VCF format.
                            with open(os.path.join(subdirectory_path, file)) as f:
                                for line in f:
                                    if not line.startswith("#"):
                                        line_split = line.split("\t")
                                        if len(line_split) != 8:
                                            handle_error_exit(
                                                f"The VCF file {file} does not have the correct number of columns. Please check that the file is in VCF format."
                                            )
                                        try:
                                            int(line_split[1])
                                            float(line_split[5])
                                        except ValueError:
                                            handle_error_exit(
                                                f"The VCF file {file} does not have the correct format. Please check that the file is in VCF format."
                                            )
                    if len(vcf_files) == 0:
                        handle_error_exit(
                            f"The input directory {subdirectory_path} does not contain any VCF files. Valid file extension is .vcf."
                        )
                elif subdirectory == "gff":
                    gff_files = []
                    for file in os.listdir(subdirectory_path):
                        if file.endswith(".gff") or file.endswith(".gff3"):
                            gff_files.append(file)
                            # Validate the GFF file. It should be in GFF format.
                            with open(os.path.join(subdirectory_path, file)) as f:
                                for line in f:
                                    if not line.startswith("#"):
                                        line_split = line.split("\t")
                                        if len(line_split) != 9:
                                            handle_error_exit(
                                                f"The GFF file {file} does not have the correct number of columns. Please check that the file is in GFF format."
                                            )
                                        try:
                                            int(line_split[3])
                                            int(line_split[4])
                                            float(line_split[5])
                                        except ValueError:
                                            handle_error_exit(
                                                f"The GFF file {file} does not have the correct format. Please check that the file is in GFF format."
                                            )
                    if len(gff_files) == 0:
                        handle_error_exit(
                            f"The input directory {subdirectory_path} does not contain any GFF files. Valid file extensions are .gff and .gff3."
                        )


def get_data_files(input_subdir: str, data_type: str) -> List[str]:
    """
    Returns the paths to the data files of the specified type in the provided subdirectory.

    Args:
        input_subdir (str): The path to the subdirectory containing the data files.
        data_type (str): The type of the data files to be returned. Valid values are "genbank", "blast", "bed", "json", "vcf", and "gff".

    Returns:
        list: The paths to the data files.
    """
    data_dir = os.path.join(input_subdir, data_type)
    if os.path.isdir(data_dir):
        if data_type == "genbank":
            genbank_paths = [
                os.path.join(data_dir, file)
                for file in os.listdir(data_dir)
                if file.endswith(".gbk")
                or file.endswith(".gbff")
                or file.endswith(".gb")
            ]
            assert (
                len(genbank_paths) == 1
            ), f"The input directory {data_dir} does not contain exactly one GenBank file. Valid file extensions are .gbk, .gbff, and .gb."
            return genbank_paths
        elif data_type == "blast":
            return [
                os.path.join(data_dir, file)
                for file in os.listdir(data_dir)
                if file.endswith(".txt") or file.endswith(".tsv")
            ]
        elif data_type == "bed":
            return [
                os.path.join(data_dir, file)
                for file in os.listdir(data_dir)
                if file.endswith(".bed")
            ]
        elif data_type == "json":
            return [
                os.path.join(data_dir, file)
                for file in os.listdir(data_dir)
                if file.endswith(".json")
            ]
        elif data_type == "vcf":
            return [
                os.path.join(data_dir, file)
                for file in os.listdir(data_dir)
                if file.endswith(".vcf")
            ]
        elif data_type == "gff":
            return [
                os.path.join(data_dir, file)
                for file in os.listdir(data_dir)
                if file.endswith(".gff") or file.endswith(".gff3")
            ]
        else:
            handle_error_exit(f"Invalid data type: {data_type}")
    else:
        if data_type == "genbank":
            handle_error_exit(f"The expected directory {data_dir} does not exist.")
        else:
            return []

    # This should never be reached.
    return []


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
