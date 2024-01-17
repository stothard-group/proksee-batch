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


def test_genbank_to_cgview_on_origin_in_codon() -> None:
    # Define path to an example genbank file (it is in the tests/data directory).
    with resources.path("tests.data.genbank", "minimal_example_3.gbk") as genbank_path:
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

            # Check that there are the correct number of CDS features in the JSON file.
            with open(json_file) as json_fh:
                json_data = json.load(json_fh)
                num_cds_features = 0
                for feature in json_data["cgview"]["features"]:
                    if feature["type"] == "CDS":
                        num_cds_features += 1
                assert (
                    num_cds_features == 4
                ), """There should be four CDS features in the JSON file,
                because one CDS feature location is split up by the origin."""

            ## Check that the length of each coding sequence is a multiple of three.
            # for feature in json_data["cgview"]["features"]:
            #    if feature["type"] == "CDS":
            #        assert (
            #            abs(feature["stop"] - feature["start"]) + 1
            #        ) % 3 == 0, """The length of each coding sequence should be a multiple
            #        of three, otherwise the reading frame may be disrupted by
            #        splitting the original coding sequence at the origin."""

            # Check that specific expected coding sequences are present.
            expected_cds_features = [
                {
                    "type": "CDS",
                    "name": "VHWIEOXL_CDS249",
                    "start": 1,
                    "stop": 239,
                    "strand": 1,
                    "source": "genbank-features",
                    "contig": "AB123456",
                    "legend": "CDS",
                    "codonStart": 3,
                },
                {
                    "type": "CDS",
                    "name": "VHWIEOXL_CDS249",
                    "start": 140888,
                    "stop": 141098,
                    "strand": 1,
                    "source": "genbank-features",
                    "contig": "AB123456",
                    "legend": "CDS",
                    "codonStart": 1,
                },
            ]
            found_expected_cds_features = []
            for feature in json_data["cgview"]["features"]:
                if feature["type"] == "CDS":
                    print(feature)
                    if feature in expected_cds_features:
                        found_expected_cds_features.append(feature)
            assert len(found_expected_cds_features) == len(
                expected_cds_features
            ), """The expected coding sequences are not present in the JSON file."""


def test_genbank_to_cgview_on_reverse_strand() -> None:
    """This uses the sample GenBank file from NCBI
    (https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html) to confirm the
    manner in which features on the reverse strand are parsed. This is important
    because the CGView JSON format requires that the start and stop positions of
    features on the reverse strand be specified in the reverse order from the
    way they are specified in the GenBank file."""
    # Define path to an example genbank file (it is in the tests/data directory).
    with resources.path("tests.data.genbank", "minimal_example_4.gbk") as genbank_path:
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

            # Check that there are the correct number of CDS features in the JSON file.
            with open(json_file) as json_fh:
                json_data = json.load(json_fh)
                num_cds_features = 0
                for feature in json_data["cgview"]["features"]:
                    if feature["type"] == "CDS":
                        num_cds_features += 1
                assert num_cds_features == 3

            # Check that the start and stop positions of the CDS for the REV7 gene are correct.
            rev7_gene_found = False
            for feature in json_data["cgview"]["features"]:
                if feature["type"] == "CDS":
                    if feature["name"] == "REV7":
                        rev7_gene_found = True
                        assert feature["start"] == 3300
                        assert feature["stop"] == 4037
                        assert feature["strand"] == -1
            assert rev7_gene_found
