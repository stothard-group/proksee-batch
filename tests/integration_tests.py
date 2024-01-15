"""Test functionality involving access to external resources including the
Proksee server."""
import json
import os
import tempfile
from importlib import resources
from pathlib import Path

import pytest
from Bio import SeqIO
from click.testing import CliRunner

from proksee_batch import __main__
from proksee_batch.generate_proksee_link import generate_proksee_link
from proksee_batch.scrape_proksee_image import scrape_proksee_image


@pytest.fixture
def runner() -> CliRunner:
    """Fixture for invoking command-line interfaces."""
    return CliRunner()


def test_download_example(runner: CliRunner) -> None:
    """
    Downloads example GenBank files from NCBI FTP site for running tests to test
    the proksee-batch CLI.
    """
    # Open a temporary directory for the example GenBank files
    with tempfile.TemporaryDirectory() as temp_dir:
        genome_dir = os.path.join(temp_dir, "genbank_files")

        # Run the CLI command
        result = runner.invoke(
            __main__.main,
            [
                "--download-example",
                genome_dir,
            ],
        )

        # Check that the command ran successfully
        assert result.exit_code == 0
        assert os.path.exists(genome_dir)

        # List all the .gbff files in the directory
        gbff_files = [f for f in os.listdir(genome_dir) if f.endswith(".gbff")]

        # Check that there is at least one .gbff file in the directory
        assert len(gbff_files) > 0

        # Parse the first .gbff file with BioPython to check that it is a valid GenBank file.
        records = SeqIO.parse(os.path.join(genome_dir, gbff_files[0]), "genbank")
        assert len(list(records)) > 0


def test_default_proksee_json() -> None:
    """Test that the default Proksee JSON template yields a valid Proksee project,
    by checking the generation of a Proksee link and the subsequent scraping of
    an SVG image, ensuring that the SVG contains the expected data."""

    # Accessing the JSON template using importlib.resources
    with resources.files("proksee_batch.data").joinpath(
        "default_proksee_template.json"
    ).open() as file:
        template_json = json.load(file)
    cgview_name = template_json["cgview"]["name"]

    # Use temporary files for Proksee link and SVG image
    with tempfile.TemporaryDirectory() as temp_dir, resources.files(
        "proksee_batch.data"
    ).joinpath("default_proksee_template.json") as template_file_path:
        temp_proksee_link = Path(temp_dir, "proksee_link.txt")
        temp_svg_file = Path(temp_dir, "proksee_image.svg")

        # Proceed with the test
        try:
            generate_proksee_link(str(template_file_path), str(temp_proksee_link))
            assert temp_proksee_link.exists(), "Proksee link file was not created."

            scrape_proksee_image(str(temp_proksee_link), str(temp_svg_file))
            assert temp_svg_file.exists(), "SVG file was not created."

            with open(temp_svg_file) as file:
                svg_contents = file.read()
            assert cgview_name in svg_contents, (
                "SVG file does not contain the cgview name. This suggests that the "
                "Proksee image was not generated correctly."
            )
        except Exception as e:
            raise AssertionError(f"An error occurred during the test: {e}")
