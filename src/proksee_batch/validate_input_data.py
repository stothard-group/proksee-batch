"""

The input directory must be structured as in the following example:

    input_directory/
        genome_name_1/
            genbank/
                genome1.gbk
            fasta/
                genome1.fna
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
            fasta/
                genome2.fa
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

The genbank directory must contain a single GenBank file with the extension
.gbk, .gbff, or .gb. This is the genome that will be visualized. If the genbank
directory is not present, then proksee-batch will use a file from the fasta
directory instead (otherwise the fasta directory is ignored).  The blast, bed,
vcf, and gff directories are optional. They contain files with additional
genomic features.  The json directory is also optional. It contains a custom
Proksee project JSON file that will be used as a template for the visualization.
"""
import json
import os
import sys
from typing import Dict
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
    # Check that the input directory exists.
    if not os.path.isdir(input):
        handle_error_exit(f"The input directory {input} does not exist.")

    # Iterate over subdirectories in the input directory (each subdirectory
    # should contain the data for one genome).
    for genome_dir in os.listdir(input):
        # Check that the subdirectory contains a genbank subdirectory, a fasta subdirectory, or both.
        required_subdirectories = ["genbank", "fasta"]
        if not any(
            os.path.isdir(os.path.join(input, genome_dir, subdirectory))
            for subdirectory in required_subdirectories
        ):
            handle_error_exit(
                f"The input directory {input} does not contain any of the required subdirectories: {required_subdirectories}"
            )

        # Check that all subdirectories have valid names.
        valid_subdirectories = [
            "genbank",
            "fasta",
            "blast",
            "bed",
            "json",
            "vcf",
            "gff",
        ]
        for subdirectory in os.listdir(os.path.join(input, genome_dir)):
            if subdirectory not in valid_subdirectories:
                handle_error_exit(
                    f"The input directory contains an invalid subdirectory: {subdirectory}"
                )

        # Check the contents of the subdirectories, if they exist.
        seq_file_path = ""
        seq_file_format = ""
        for subdirectory in valid_subdirectories:
            subdirectory_path = os.path.join(input, genome_dir, subdirectory)
            if os.path.isdir(subdirectory_path):
                if subdirectory == "genbank":
                    genbank_files = []
                    for file in os.listdir(subdirectory_path):
                        if (
                            file.endswith(".gbk")
                            or file.endswith(".gbff")
                            or file.endswith(".gb")
                        ):
                            genbank_files.append(file)
                            # Parse the GenBank file with BioPython to check if the file format is valid.
                            try:
                                SeqIO.parse(
                                    os.path.join(subdirectory_path, file), "genbank"
                                ).__next__()  # type: ignore
                            except Exception as e:
                                handle_error_exit(
                                    f"Error parsing the GenBank file {file}: {e}. Please check that the file is in GenBank format."
                                )
                    if len(genbank_files) == 0:
                        handle_error_exit(
                            f"The input directory {subdirectory_path} does not contain any GenBank files. Valid file extensions are .gbk, .gbff, and .gb."
                        )
                    elif len(genbank_files) > 1:
                        handle_error_exit(
                            f"The input directory {subdirectory_path} contains more than one GenBank file. Please ensure that there is exactly one GenBank file in the directory."
                        )
                    seq_file_path = os.path.join(subdirectory_path, genbank_files[0])
                    seq_file_format = "genbank"
                elif subdirectory == "fasta":
                    fasta_files = []
                    for file in os.listdir(subdirectory_path):
                        if (
                            file.endswith(".fna")
                            or file.endswith(".fa")
                            or file.endswith(".fasta")
                        ):
                            fasta_files.append(file)
                            # Parse the FASTA file with BioPython to check if the file format is valid.
                            try:
                                SeqIO.parse(
                                    os.path.join(subdirectory_path, file), "fasta"
                                ).__next__()  # type: ignore
                            except Exception as e:
                                handle_error_exit(
                                    f"Error parsing the FASTA file {file}: {e}. Please check that the file is in FASTA format."
                                )
                            # Check that the sequences are nucleotide sequences, and not protein sequences.
                            with open(os.path.join(subdirectory_path, file)) as f:
                                for line in f:
                                    if not line.startswith(">") and line.strip() != "":
                                        if not set(line.strip().upper()).issubset(
                                            set("ATGCNatgcn")
                                        ):
                                            handle_error_exit(
                                                f"The FASTA file {file} contains one-letter codes that do not represent DNA nucleotide bases. Please check that the file contains nucleotide sequences only."
                                            )
                    if len(fasta_files) == 0:
                        handle_error_exit(
                            f"The input directory {subdirectory_path} does not contain any FASTA files. Valid file extensions are .fna and .fa."
                        )
                    elif len(fasta_files) > 1:
                        handle_error_exit(
                            f"The input directory {subdirectory_path} contains more than one FASTA file. Please ensure that there is exactly one FASTA file in the directory."
                        )
                    if seq_file_path == "":
                        seq_file_path = os.path.join(subdirectory_path, fasta_files[0])
                        seq_file_format = "fasta"
                elif subdirectory == "blast":
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
                                            try:
                                                float(line_split[5])
                                            except:
                                                assert line_split[5] == "."
                                        except ValueError:
                                            handle_error_exit(
                                                f"The VCF file {file} does not have the correct format. Please check that the file is in VCF format."
                                            )
                            # Check that all sequence IDs in the VCF file are contigs in the GenBank file.
                            if not check_vcf_seq_ids(
                                os.path.join(subdirectory_path, file),
                                seq_file_path,
                                seq_file_format,
                            ):
                                handle_error_exit(
                                    f"The VCF file {file} contains sequence IDs that are not contigs in the GenBank file {genbank_files[0]}."
                                )
                            # Check that the genotypes in the genome in the GenBank file match the REF genotypes in the VCF file.
                            if not check_vcf_ref_vs_alt_genotypes(
                                os.path.join(subdirectory_path, file),
                                seq_file_path,
                                seq_file_format,
                            ):
                                handle_error_exit(
                                    f"The genotypes in the genome in the GenBank file {genbank_files[0]} do not match the REF genotypes in the VCF file {file}."
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
                                            try:
                                                float(line_split[5])
                                            except:
                                                assert line_split[5] == "."
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
        elif data_type == "fasta":
            fasta_paths = [
                os.path.join(data_dir, file)
                for file in os.listdir(data_dir)
                if file.endswith(".fna")
                or file.endswith(".fa")
                or file.endswith(".fasta")
            ]
            assert (
                len(fasta_paths) == 1
            ), f"The input directory {data_dir} does not contain exactly one FASTA file. Valid file extensions are .fna, .fa, and .fasta."
            return fasta_paths
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


def check_vcf_seq_ids(
    vcf_file_path: str, seq_file_path: str, seq_file_format: str
) -> bool:
    """
    Checks if all the sequence IDs in the first column of the VCF file are contigs in a GenBank or FASTA file.

    Args:
        vcf_file_path (str): The path to the VCF file.
        seq_file_path (str): The path to the GenBank or FASTA file.
        seq_file_format (str): The format of the GenBank or FASTA file. Valid values are "genbank" and "fasta".

    Returns:
        bool: True if all sequence IDs in the VCF file are contigs in the sequence file, False otherwise.
    """
    assert seq_file_format in ["genbank", "fasta"]
    vcf_seq_ids = []
    with open(vcf_file_path) as vcf_file:
        for line in vcf_file:
            if not line.startswith("#") and line.strip() != "":
                vcf_seq_ids.append(line.split("\t")[0])
    seq_ids = []
    for record in SeqIO.parse(seq_file_path, seq_file_format):  # type: ignore
        seq_ids.append(record.id)
    return all(seq_id in seq_ids for seq_id in vcf_seq_ids)


def check_vcf_ref_vs_alt_genotypes(
    vcf_file_path: str, genome_file_path: str, genome_file_type: str
) -> bool:
    """
    Checks if the genotypes in the genome in the GenBank file match the REF genotypes in the VCF file.
    Args:
        vcf_file_path (str): The path to the VCF file.
        genbank_file_path (str): The path to the GenBank file.
    Returns:
        bool: True if the genotypes in the genome in the GenBank file match the REF genotypes in the VCF file, False otherwise.
    """
    # Get the genome sequences from the GenBank file
    genome_sequences = {}
    if genome_file_type == "genbank":
        genome_sequences = {
            record.id: str(record.seq)
            for record in SeqIO.parse(genome_file_path, "genbank")  # type: ignore
        }
    elif genome_file_type == "fasta":
        genome_sequences = {
            record.id: str(record.seq)
            for record in SeqIO.parse(genome_file_path, "fasta")  # type: ignore
        }
    else:
        handle_error_exit(f"Invalid genome file type: {genome_file_type}")

    # Parse the VCF file content
    vcf_data = []
    with open(vcf_file_path) as vcf_file:
        for line in vcf_file:
            if not line.startswith("#") and line.strip():
                parts = line.split("\t")
                vcf_data.append((parts[0], int(parts[1]), parts[3]))

    # Check each genotype in the VCF file
    for chrom, pos, ref_genotype in vcf_data:
        actual_genotype = genome_sequences[chrom][pos - 1]
        if ref_genotype != actual_genotype:
            print(
                f"Genotype mismatch at {chrom}:{pos}: {actual_genotype} should be {ref_genotype}"
            )
            return False

    return True
