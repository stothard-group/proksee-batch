import json
import os
import tempfile
from importlib import resources

from proksee_batch.genbank_to_cgview_json import genbank_to_cgview_json


def test_genbank_to_cgview_json() -> None:
    # Define path to an example genbank file (it is in the tests/data directory).
    with resources.path("tests.data.genbank", "minimal_example_1.gbk") as genbank_path:
        genbank_file = str(genbank_path)

        # Define path for a temporary output directory.
        with tempfile.TemporaryDirectory() as temp_output:
            # Define path to temporary json file.
            json_file = os.path.join(temp_output, "temp.json")

            # Use the genbank_to_cgview_json function to convert the genbank file to a CGView JSON file.
            genbank_to_cgview_json(genbank_file, json_file)

            # Check that the JSON file was created.
            assert os.path.exists(json_file)

            # Check that the JSON file is not empty.
            assert os.path.getsize(json_file) > 0

            # Check that there is one, and only one, CDS feature in the JSON file.
            with open(json_file) as json_fh:
                json_data = json.load(json_fh)
                num_cds_features = 0
                for feature in json_data["cgview"]["features"]:
                    if feature["type"] == "CDS":
                        num_cds_features += 1
                assert num_cds_features == 1


def test_genbank_to_cgview_json_on_compound_location() -> None:
    # Define path to an example genbank file (it is in the tests/data directory).
    with resources.path("tests.data.genbank", "minimal_example_2.gbk") as genbank_path:
        genbank_file = str(genbank_path)

        # Define path for a temporary output directory.
        with tempfile.TemporaryDirectory() as temp_output:
            # Define path to temporary json file.
            json_file = os.path.join(temp_output, "temp.json")

            # Use the genbank_to_cgview_json function to convert the genbank file to a CGView JSON file.
            genbank_to_cgview_json(genbank_file, json_file)

            # Check that the JSON file was created.
            assert os.path.exists(json_file)

            # Check that the JSON file is not empty.
            assert os.path.getsize(json_file) > 0

            # Check that there is one, and only one, CDS feature in the JSON file.
            with open(json_file) as json_fh:
                json_data = json.load(json_fh)
                num_cds_features = 0
                for feature in json_data["cgview"]["features"]:
                    if feature["type"] == "CDS":
                        num_cds_features += 1
                assert (
                    num_cds_features == 2
                ), "There should be two CDS features in the JSON file, because the CDS feature location is a CompoundLocation with two exons."
