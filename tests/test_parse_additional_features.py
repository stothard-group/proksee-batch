import json
import os
import tempfile
from importlib import resources

from proksee_batch.parse_additional_features import add_bed_features_and_tracks
from proksee_batch.parse_additional_features import add_blast_features_and_tracks
from proksee_batch.parse_additional_features import parse_bed_files
from proksee_batch.parse_additional_features import parse_blast_files


def test_parse_blast_files() -> None:
    """Test the parse_blast_files function using a simple BLAST result file."""
    # Define path to an example BLAST result file (in oufmt 6 format).
    with resources.path("tests.data.blast", "U49845.1_subject1.txt") as blast_path:
        blast_file = str(blast_path)

        # Use the parse_blast_files function to parse the BLAST result file.
        blast_features, blast_tracks = parse_blast_files([blast_file])

        # Assert that the BLAST result file was parsed correctly.
        expected_blast_features = [
            {
                "name": "subject_seq_1",
                "type": "blast",
                "start": 11,
                "stop": 70,
                # "strand": 1,
                "strand": ".",
                "source": "blast_1",
                "contig": "U49845.1",
                "legend": "U49845.1_subject1.txt",
                "tags": [],
                "meta": {
                    "query": "U49845.1",
                    "query_start": 11,
                    "query_stop": 70,
                    "subject": "subject_seq_1",
                    "subject_start": 1,
                    "subject_stop": 60,
                    "alignment_length": 60,
                    "identity": 100.0,
                    "mismatches": 0,
                    "evalue": 5.01e-29,
                    "bit_score": 111.0,
                },
            }
        ]
        expected_blast_tracks = [
            {
                "name": "U49845.1_subject1.txt",
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
    """Test the add_blast_features_and_tracks function using a simple BLAST
    result file, and a simple cgview map file in JSON format."""
    # Use the add_blast_features_and_tracks function to parse the BLAST result
    # file, a cgview map JSON file, and write a new cgview map file in JSON
    # format.
    with resources.path(
        "tests.data.blast", "U49845.1_subject1.txt"
    ) as blast_path, resources.path(
        "tests.data.json", "U49845.1.json"
    ) as json_path, tempfile.TemporaryDirectory() as temp_dir:
        blast_file = str(blast_path)
        json_file = str(json_path)
        output_file = os.path.join(temp_dir, "U49845.1.json")

        # Use the add_blast_features_and_tracks function to parse the BLAST
        # result file, a cgview map JSON file, and write a new cgview map file
        # in JSON format.
        add_blast_features_and_tracks([blast_file], json_file, output_file)

        ## Define expected BLAST features and tracks.
        expected_blast_features = [
            {
                "name": "subject_seq_1",
                "type": "blast",
                "start": 11,
                "stop": 70,
                # "strand": 1,
                "strand": ".",
                "source": "blast_1",
                "contig": "U49845.1",
                "legend": "U49845.1_subject1.txt",
                "tags": [],
                "meta": {
                    "query": "U49845.1",
                    "query_start": 11,
                    "query_stop": 70,
                    "subject": "subject_seq_1",
                    "subject_start": 1,
                    "subject_stop": 60,
                    "alignment_length": 60,
                    "identity": 100.0,
                    "mismatches": 0,
                    "evalue": 5.01e-29,
                    "bit_score": 111.0,
                },
            }
        ]
        expected_blast_tracks = [
            {
                "name": "U49845.1_subject1.txt",
                "separateFeaturesBy": "none",
                "position": "both",
                "thicknessRatio": 1,
                "dataType": "feature",
                "dataMethod": "source",
                "dataKeys": "blast_1",
                "drawOrder": "score",
            }
        ]

        # Check that the expected BLAST features and tracks were correctly added
        # to the cgview map JSON data structure.
        with open(output_file) as f:
            json_data = json.load(f)
            for expected_feature in expected_blast_features:
                assert expected_feature in json_data["cgview"]["features"]
            for expected_track in expected_blast_tracks:
                assert expected_track in json_data["cgview"]["tracks"]


def test_add_blast_features_and_tracks_multiple_features() -> None:
    """Test the add_blast_features_and_tracks function using a simple BLAST
    result file, and a simple cgview map file in JSON format. There are multiple
    features to map on a single track."""
    # Use the add_blast_features_and_tracks function to parse the BLAST result
    # file, a cgview map JSON file, and write a new cgview map file in JSON
    # format.
    with resources.path(
        "tests.data.blast", "U49845.1_subject2.txt"
    ) as blast_path, resources.path(
        "tests.data.json", "U49845.1.json"
    ) as json_path, tempfile.TemporaryDirectory() as temp_dir:
        blast_file = str(blast_path)
        json_file = str(json_path)
        output_file = os.path.join(temp_dir, "U49845.1.json")

        # Use the add_blast_features_and_tracks function to parse the BLAST
        # result file, a cgview map JSON file, and write a new cgview map file
        # in JSON format.
        add_blast_features_and_tracks([blast_file], json_file, output_file)

        ## Define expected BLAST features and tracks.
        # expected_blast_features = [
        #    {
        #        "name": "U49845.1",
        #        "type": "blast",
        #        "start": 1,
        #        "stop": 5028,
        #        # "strand": 1,
        #        "strand": ".",
        #        "source": "blast_1",
        #        "contig": "U49845.1",
        #        "legend": "U49845.1_v2.txt",
        #        "tags": [],
        #        "meta": {
        #            "query": "U49845.1",
        #            "query_start": 1,
        #            "query_stop": 5028,
        #            "alignment_length": 5028,
        #            "identity": 100.0,
        #            "mismatches": 0,
        #            "evalue": 0.0,
        #            "bit_score": 9286,
        #        },
        #    },
        #    {
        #        "name": "XXXXXX.1",
        #        "type": "blast",
        #        "start": 1,
        #        "stop": 100,
        #        # "strand": 1,
        #        "strand": ".",
        #        "source": "blast_1",
        #        "contig": "U49845.1",
        #        "legend": "U49845.1_v2.txt",
        #        "tags": [],
        #        "meta": {
        #            "query": "XXXXXX.1",
        #            "query_start": 1,
        #            "query_stop": 100,
        #            "alignment_length": 100,
        #            "identity": 100.0,
        #            "mismatches": 0,
        #            "evalue": 0.0,
        #            "bit_score": 9287,
        #        },
        #    },
        # ]
        # expected_blast_tracks = [
        #    {
        #        "name": "U49845.1_v2.txt",
        #        "separateFeaturesBy": "none",
        #        "position": "both",
        #        "thicknessRatio": 1,
        #        "dataType": "feature",
        #        "dataMethod": "source",
        #        "dataKeys": "blast_1",
        #        "drawOrder": "score",
        #    }
        # ]
        expected_blast_features = [
            {
                "name": "subject_seq_1",
                "type": "blast",
                "start": 11,
                "stop": 70,
                # "strand": 1,
                "strand": ".",
                "source": "blast_1",
                "contig": "U49845.1",
                "legend": "U49845.1_subject2.txt",
                "tags": [],
                "meta": {
                    "query": "U49845.1",
                    "query_start": 11,
                    "query_stop": 70,
                    "subject": "subject_seq_1",
                    "subject_start": 1,
                    "subject_stop": 60,
                    "alignment_length": 60,
                    "identity": 100.0,
                    "mismatches": 0,
                    "evalue": 9.81e-29,
                    "bit_score": 111.0,
                },
            },
            {
                "name": "subject_seq_2",
                "type": "blast",
                "start": 81,
                "stop": 140,
                # "strand": 1,
                "strand": ".",
                "source": "blast_1",
                "contig": "U49845.1",
                "legend": "U49845.1_subject2.txt",
                "tags": [],
                "meta": {
                    "query": "U49845.1",
                    "query_start": 81,
                    "query_stop": 140,
                    "subject": "subject_seq_2",
                    "subject_start": 1,
                    "subject_stop": 60,
                    "alignment_length": 60,
                    "identity": 100.0,
                    "mismatches": 0,
                    "evalue": 9.81e-29,
                    "bit_score": 111.0,
                },
            },
        ]
        expected_blast_tracks = [
            {
                "name": "U49845.1_subject2.txt",
                "separateFeaturesBy": "none",
                "position": "both",
                "thicknessRatio": 1,
                "dataType": "feature",
                "dataMethod": "source",
                "dataKeys": "blast_1",
                "drawOrder": "score",
            }
        ]

        # Check that the expected BLAST features and tracks were correctly added
        # to the cgview map JSON data structure.
        with open(output_file) as f:
            json_data = json.load(f)
            for expected_feature in expected_blast_features:
                print()
                print()
                print("Expected feature:")
                print(expected_feature)
                print()
                for feature in json_data["cgview"]["features"]:
                    print("Actual feature:")
                    print(feature)
                print()
                print()
                assert (
                    expected_feature in json_data["cgview"]["features"]
                ), f"Expected feature not found in actual features"
            for expected_track in expected_blast_tracks:
                assert (
                    expected_track in json_data["cgview"]["tracks"]
                ), f"Expected track not found in actual tracks"


def test_parse_bed_files() -> None:
    """Test the parse_bed_files function using a simple BED file."""
    # Define path to an example BED file.
    with resources.path("tests.data.bed", "minimal_example_4.bed") as bed_path:
        bed_file = str(bed_path)

        # Use the parse_bed_files function to parse the BED file.
        bed_features, bed_tracks = parse_bed_files([bed_file])

        # Assert that the BED file was parsed correctly.
        expected_bed_features = [
            {
                "name": "FirstFeature",
                "type": "bed",
                "start": 1,
                "stop": 100,
                "strand": 1,
                "source": "bed_1",
                "contig": "U49845.1",
                "legend": "minimal_example_4.bed",
                "tags": [],
                "meta": {
                    "score": 0,
                },
            },
            {
                "name": "SecondFeature",
                "type": "bed",
                "start": 101,
                "stop": 200,
                "strand": 1,
                "source": "bed_1",
                "contig": "U49845.1",
                "legend": "minimal_example_4.bed",
                "tags": [],
                "meta": {
                    "score": 0,
                },
            },
            {
                "name": "ThirdFeature",
                "type": "bed",
                "start": 201,
                "stop": 300,
                "strand": -1,
                "source": "bed_1",
                "contig": "U49845.1",
                "legend": "minimal_example_4.bed",
                "tags": [],
                "meta": {
                    "score": 0,
                },
            },
            {
                "name": "FourthFeature",
                "type": "bed",
                "start": 301,
                "stop": 400,
                "strand": -1,
                "source": "bed_1",
                "contig": "U49845.1",
                "legend": "minimal_example_4.bed",
                "tags": [],
                "meta": {
                    "score": 0,
                },
            },
        ]
        expected_bed_tracks = [
            {
                "name": "minimal_example_4.bed",
                "separateFeaturesBy": "none",
                "position": "both",
                "thicknessRatio": 1,
                "dataType": "feature",
                "dataMethod": "source",
                "dataKeys": "bed_1",
                "drawOrder": "score",
            }
        ]

        assert bed_features == expected_bed_features
        assert bed_tracks == expected_bed_tracks


def test_add_bed_features_and_tracks() -> None:
    """Test the add_bed_features_and_tracks function using a simple BED file,
    and a simple cgview map file in JSON format."""
    # Use the add_bed_features_and_tracks function to parse the BED file, a
    # cgview map JSON file, and write a new cgview map file in JSON format.
    with resources.path(
        "tests.data.bed", "minimal_example_4.bed"
    ) as bed_path, resources.path(
        "tests.data.json", "U49845.1.json"
    ) as json_path, tempfile.TemporaryDirectory() as temp_dir:
        bed_file = str(bed_path)
        json_file = str(json_path)
        output_file = os.path.join(temp_dir, "U49845.1.json")

        # Use the add_bed_features_and_tracks function to parse the BED file, a
        # cgview map JSON file, and write a new cgview map file in JSON format.
        add_bed_features_and_tracks([bed_file], json_file, output_file)

        # Define expected BED features and tracks.
        expected_bed_features = [
            {
                "name": "FirstFeature",
                "type": "bed",
                "start": 1,
                "stop": 100,
                "strand": 1,
                "source": "bed_1",
                "contig": "U49845.1",
                "legend": "minimal_example_4.bed",
                "tags": [],
                "meta": {
                    "score": 0,
                },
            },
            {
                "name": "SecondFeature",
                "type": "bed",
                "start": 101,
                "stop": 200,
                "strand": 1,
                "source": "bed_1",
                "contig": "U49845.1",
                "legend": "minimal_example_4.bed",
                "tags": [],
                "meta": {
                    "score": 0,
                },
            },
            {
                "name": "ThirdFeature",
                "type": "bed",
                "start": 201,
                "stop": 300,
                "strand": -1,
                "source": "bed_1",
                "contig": "U49845.1",
                "legend": "minimal_example_4.bed",
                "tags": [],
                "meta": {
                    "score": 0,
                },
            },
            {
                "name": "FourthFeature",
                "type": "bed",
                "start": 301,
                "stop": 400,
                "strand": -1,
                "source": "bed_1",
                "contig": "U49845.1",
                "legend": "minimal_example_4.bed",
                "tags": [],
                "meta": {
                    "score": 0,
                },
            },
        ]
        expected_bed_tracks = [
            {
                "name": "minimal_example_4.bed",
                "separateFeaturesBy": "none",
                "position": "both",
                "thicknessRatio": 1,
                "dataType": "feature",
                "dataMethod": "source",
                "dataKeys": "bed_1",
                "drawOrder": "score",
            }
        ]

        # Check that the expected BED features and tracks were correctly added
        # to the cgview map JSON data structure.
        with open(output_file) as f:
            json_data = json.load(f)
            for expected_feature in expected_bed_features:
                assert expected_feature in json_data["cgview"]["features"]
            for expected_track in expected_bed_tracks:
                assert expected_track in json_data["cgview"]["tracks"]
