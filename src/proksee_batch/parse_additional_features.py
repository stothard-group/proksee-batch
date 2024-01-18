import json
import os
from typing import Any
from typing import Dict
from typing import List
from typing import Tuple


def parse_blast_files(
    blast_files: List[str],
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """
    Parses BLAST result files.

    Parameters:
    blast_files (list): A list of paths to BLAST result files.

    Returns:
    tuple: A tuple containing a list of parsed BLAST features and a list of parsed BLAST tracks.
    """
    blast_features = []
    blast_tracks = []
    for i, blast_file in enumerate(blast_files):
        num = i + 1

        # Define the BLAST track.
        blast_track = {
            "name": os.path.basename(blast_file),
            "separateFeaturesBy": "none",
            "position": "both",
            "thicknessRatio": 1,
            "dataType": "feature",
            "dataMethod": "source",
            "dataKeys": f"blast_{num}",
            "drawOrder": "score",
        }
        # Add the BLAST track to the list of BLAST tracks.
        blast_tracks.append(blast_track)

        # Parse the BLAST result file.
        with open(blast_file) as f:
            blast_results = f.readlines()
        # Parse the BLAST result file.
        for blast_result_line in blast_results:
            # Parse the BLAST result line.
            blast_result = blast_result_line.split("\t")

            # Skip the header line(s), if present.
            if blast_result[0] == "qseqid" or blast_result[0].startswith("#"):
                continue

            # Define the BLAST feature.
            blast_feature = {
                "name": blast_result[
                    0
                ],  # Here we assume that the query sequence is a sequence aligned to a contig in the genome/map.
                "type": "blast",
                "start": int(
                    blast_result[8]
                ),  # Here we assume that the subject sequence is a contig in the genome/map.
                "stop": int(
                    blast_result[9]
                ),  # Here we assume that the subject sequence is a contig in the genome/map.
                "strand": 1 if blast_result[8] < blast_result[9] else -1,
                "source": f"blast_{num}",
                "contig": blast_result[
                    1
                ],  # Here we assume that the subject sequence is a contig in the genome/map.
                "legend": os.path.basename(blast_file),
                "tags": [],
                "meta": {
                    "identity": float(blast_result[2]),
                    "mismatches": int(blast_result[4]),
                    "evalue": float(blast_result[10]),
                    "bit_score": int(blast_result[11]),
                },
            }
            # Add the BLAST feature to the list of BLAST features.
            blast_features.append(blast_feature)
    assert (
        len(blast_features) > 0 and len(blast_tracks) > 0
    ), "No BLAST features or tracks were obtained from input file {}.".format(
        blast_file
    )
    return (blast_features, blast_tracks)


def add_blast_features_and_tracks(
    blast_files: List[str], json_file: str, output_file: str
) -> None:
    """
    Parses BLAST result files, adds the parsed BLAST features and tracks to the cgview map JSON data structure, and writes the cgview map JSON data structure to a new file.

    Parameters:
    blast_files (list): A list of paths to BLAST result files.
    json_file (str): The path to a cgview map JSON file.
    output_file (str): The path to the output file.

    Returns:
    None
    """
    # Parse the BLAST result file.
    blast_features, blast_tracks = parse_blast_files(blast_files)
    # Read the cgview map JSON file.
    with open(json_file) as f:
        json_data = json.load(f)

    # Get list of contigs from JSON data.
    contigs = [contig["name"] for contig in json_data["cgview"]["sequence"]["contigs"]]

    # Check that all the blast features are assigned to contigs that are in the JSON data.
    for blast_feature in blast_features:
        assert (
            blast_feature["contig"] in contigs
        ), "The subject sequence {} listed in the BLAST result file {} is not among the contigs in the genome being mapped.".format(
            blast_feature["contig"], blast_feature["legend"]
        )

    # Add the parsed BLAST features and tracks to the cgview map JSON data structure.
    json_data["cgview"]["features"] += blast_features
    json_data["cgview"]["tracks"] += blast_tracks
    # Write the cgview map JSON data structure to a new file.
    with open(output_file, "w") as f:
        json.dump(json_data, f, indent=4)


def parse_bed_files(
    bed_files: List[str],
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """
    Parses BED files.

    Parameters:
    bed_files (list): A list of paths to BED files.

    Returns:
    tuple: A tuple containing a list of parsed BED features and a list of parsed BED tracks.
    """
    bed_features = []
    bed_tracks = []
    for i, bed_file in enumerate(bed_files):
        num = i + 1

        # Define the BED track.
        bed_track = {
            "name": os.path.basename(bed_file),
            "separateFeaturesBy": "none",
            "position": "both",
            "thicknessRatio": 1,
            "dataType": "feature",
            "dataMethod": "source",
            "dataKeys": f"bed_{num}",
            "drawOrder": "score",
        }
        # Add the BED track to the list of BED tracks.
        bed_tracks.append(bed_track)

        # Parse the BED file.
        with open(bed_file) as f:
            bed_results = f.readlines()
        # Parse the BED file.
        for bed_result_line in bed_results:
            # Parse the BED result line.
            bed_result = bed_result_line.strip().split()

            # Skip the header line(s), if present.
            if (
                bed_result[0].startswith("#")
                or bed_result[0].startswith("track ")
                or bed_result[0].startswith("browser ")
            ):
                continue

            assert (
                len(bed_result) >= 3
            ), f"BED file {bed_file} is not in the correct format."

            name = ""
            if len(bed_result) >= 4:
                name = bed_result[3]
            score = 0
            if len(bed_result) >= 5:
                score = int(bed_result[4])
            strand = 1
            if len(bed_result) >= 6:
                if bed_result[5] == "+":
                    strand = 1
                elif bed_result[5] == "-":
                    strand = -1
                else:
                    raise ValueError(
                        "The strand column of BED file {} is not in the correct format.".format(
                            bed_file
                        )
                    )

            # Define the BED feature.
            bed_feature = {
                "name": name,
                "type": "bed",
                "start": int(bed_result[1]) + 1,
                "stop": int(bed_result[2]),
                "strand": strand,
                "source": f"bed_{num}",
                "contig": bed_result[0],
                "legend": os.path.basename(bed_file),
                "tags": [],
                "meta": {
                    "score": score,
                },
            }
            # Add the BED feature to the list of BED features.
            bed_features.append(bed_feature)
    assert (
        len(bed_features) > 0 and len(bed_tracks) > 0
    ), f"No BED features or tracks were obtained from input file {bed_file}."
    return (bed_features, bed_tracks)


def add_bed_features_and_tracks(
    bed_files: List[str], json_file: str, output_file: str
) -> None:
    """
    Parses BED files, adds the parsed BED features and tracks to the cgview map JSON data structure, and writes the cgview map JSON data structure to a new file.

    Parameters:
    bed_files (list): A list of paths to BED files.
    json_file (str): The path to a cgview map JSON file.
    output_file (str): The path to the output file.

    Returns:
    None
    """
    # Parse the BED file.
    bed_features, bed_tracks = parse_bed_files(bed_files)
    # Read the cgview map JSON file.
    with open(json_file) as f:
        json_data = json.load(f)

    # Get list of contigs from JSON data.
    contigs = [contig["name"] for contig in json_data["cgview"]["sequence"]["contigs"]]

    # Check that all the bed features are assigned to contigs that are in the JSON data.
    for bed_feature in bed_features:
        assert (
            bed_feature["contig"] in contigs
        ), "The contig {} listed in the BED file {} is not among the contigs in the genome being mapped.".format(
            bed_feature["contig"], bed_feature["legend"]
        )

    # Add the parsed BED features and tracks to the cgview map JSON data structure.
    json_data["cgview"]["features"] += bed_features
    json_data["cgview"]["tracks"] += bed_tracks
    # Write the cgview map JSON data structure to a new file.
    with open(output_file, "w") as f:
        json.dump(json_data, f, indent=4)
