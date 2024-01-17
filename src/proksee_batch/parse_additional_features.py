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
    for blast_file in blast_files:
        # Parse the BLAST result file.
        with open(blast_file) as f:
            blast_results = f.readlines()
        # Parse the BLAST result file.
        for i, blast_result_line in enumerate(blast_results):
            num = i + 1
            # Parse the BLAST result line.
            blast_result = blast_result_line.split("\t")
            # Define the BLAST feature.
            blast_feature = {
                "name": "",
                "type": "blast",
                "start": int(blast_result[8]),
                "stop": int(blast_result[9]),
                "strand": 1 if blast_result[8] < blast_result[9] else -1,
                "source": f"blast_{num}",
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
    # Add the parsed BLAST features and tracks to the cgview map JSON data structure.
    json_data["cgview"]["features"] += blast_features
    json_data["cgview"]["tracks"] += blast_tracks
    # Write the cgview map JSON data structure to a new file.
    with open(output_file, "w") as f:
        json.dump(json_data, f, indent=4)
