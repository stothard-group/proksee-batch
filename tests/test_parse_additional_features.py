import json
import os
import tempfile
from importlib import resources

from proksee_batch.parse_additional_features import add_blast_features_and_tracks
from proksee_batch.parse_additional_features import parse_blast_files


def test_parse_blast_files() -> None:
    """Test the parse_blast_files function using a simple BLAST result file."""
    # Define path to an example BLAST result file (in oufmt 6 format).
    with resources.path("tests.data.blast", "U49845.1.txt") as blast_path:
        blast_file = str(blast_path)

        # Use the parse_blast_files function to parse the BLAST result file.
        blast_features, blast_tracks = parse_blast_files([blast_file])

        # Assert that the BLAST result file was parsed correctly.
        expected_blast_features = [
            {
                "name": "",
                "type": "blast",
                "start": 1,
                "stop": 5028,
                "strand": 1,
                "source": "blast_1",
                "legend": "U49845.1.txt",
                "tags": [],
                "meta": {
                    "identity": 100.0,
                    "mismatches": 0,
                    "evalue": 0.0,
                    "bit_score": 9286,
                },
            }
        ]
        expected_blast_tracks = [
            {
                "name": "U49845.1.txt",
                "separateFeaturesBy": "none",
                "position": "both",
                "thicknessRatio": 1,
                "dataType": "feature",
                "dataMethod": "source",
                "dataKeys": "blast_1",
                "drawOrder": "score",
            }
        ]

        assert blast_features == expected_blast_features
        assert blast_tracks == expected_blast_tracks


def test_add_blast_features_and_tracks() -> None:
    """Test the add_blast_features_and_tracks function using a simple BLAST result file, and a simple cgview map file in JSON format."""
    # Use the add_blast_features_and_tracks function to parse the BLAST result file, a cgview map JSON file, and write a new cgview map file in JSON format.
    with resources.path(
        "tests.data.blast", "U49845.1.txt"
    ) as blast_path, resources.path(
        "tests.data.json", "U49845.1.json"
    ) as json_path, tempfile.TemporaryDirectory() as temp_dir:
        blast_file = str(blast_path)
        json_file = str(json_path)
        output_file = os.path.join(temp_dir, "U49845.1.json")

        # Use the add_blast_features_and_tracks function to parse the BLAST result file, a cgview map JSON file, and write a new cgview map file in JSON format.
        add_blast_features_and_tracks([blast_file], json_file, output_file)

        # Define expected BLAST features and tracks.
        expected_blast_features = [
            {
                "name": "",
                "type": "blast",
                "start": 1,
                "stop": 5028,
                "strand": 1,
                "source": "blast_1",
                "legend": "U49845.1.txt",
                "tags": [],
                "meta": {
                    "identity": 100.0,
                    "mismatches": 0,
                    "evalue": 0.0,
                    "bit_score": 9286,
                },
            }
        ]
        expected_blast_tracks = [
            {
                "name": "U49845.1.txt",
                "separateFeaturesBy": "none",
                "position": "both",
                "thicknessRatio": 1,
                "dataType": "feature",
                "dataMethod": "source",
                "dataKeys": "blast_1",
                "drawOrder": "score",
            }
        ]

        # Check that the expected BLAST features and tracks were correctly added to the cgview map JSON data structure.
        with open(output_file) as f:
            json_data = json.load(f)
            for expected_feature in expected_blast_features:
                assert expected_feature in json_data["cgview"]["features"]
            for expected_track in expected_blast_tracks:
                assert expected_track in json_data["cgview"]["tracks"]
