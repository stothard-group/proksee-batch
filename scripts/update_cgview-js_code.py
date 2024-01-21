# Script for updating the javascript code from the CGView.js project
# https://github.com/stothard-group/cgview-js. This code is open source, and is
# used to generate interactive visualizations in the HTML report output by
# proksee-batch, similar to the maps on the proksee.ca server (which also uses
# CGView.js).
#
# Usage:
# poetry run python3 scripts/update_cgview-js_code.py

# Steps:
#
# Download the latest CGView.js release from GitHub.
#
# Extract the zip file.
#
# Check whether the expected files are present.
#
# Check whether any of them have been modified since the last update (by
# comparing the md5sums to those of the previous versions of the files that are
# already in the src/proksee-batch/data/cgview-js directory).
#
# If any of the files have been modified, then copy them to the
# src/proksee-batch/data/cgview-js directory.

# Note: the directory structure of the CGView.js project should be preserved
# within the src/proksee-batch/data/cgview-js directory.

# Note: These are the relevant files from the CGView.js project:
#
# ├── LICENSE
# ├── docs
#     ├── dist
#     │   ├── cgview.css
#     │   └── cgview.min.js
#     ├── scripts
#     │   ├── bootstrap.min.js
#     │   ├── controls.js
#     │   ├── d3.min.js
#     │   ├── general.js
#     │   ├── marked.min.js
#     │   ├── prism.js
#     │   └── svgcanvas.iife.js
#     ├── styles
#         ├── bootstrap.min.css
#         ├── controls.css
#         ├── general.css
#         ├── prism.css

# Note: The script should download the latest version of the CGView.js code to a temporary directory

import hashlib
import os
import shutil
import sys
import zipfile
from pathlib import Path
from typing import Any
from typing import Dict
from typing import List
from typing import Optional
from typing import Union

import requests


# Define a list of the expected files
expected_files = [
    "LICENSE",
    "docs",
    "docs/dist",
    "docs/dist/cgview.css",
    "docs/dist/cgview.min.js",
    "docs/scripts",
    "docs/scripts/bootstrap.min.js",
    "docs/scripts/controls.js",
    "docs/scripts/d3.min.js",
    "docs/scripts/general.js",
    "docs/scripts/marked.min.js",
    "docs/scripts/prism.js",
    "docs/scripts/svgcanvas.iife.js",
    "docs/styles",
    "docs/styles/bootstrap.min.css",
    "docs/styles/controls.css",
    "docs/styles/general.css",
    "docs/styles/prism.css",
]

license_ending = """
   Copyright [yyyy] [name of copyright owner]

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

license_to_add = """/*
   Copyright 2018 Jason R. Grant

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.


   See the source repository for more information:
   https://github.com/stothard-group/cgview-js
*/
"""


def download_cgview_js(url: str, output_dir: str) -> None:
    # Define the path to the zip file
    zip_file = output_dir + ".zip"

    # Download the zip file
    print("Downloading CGView.js from GitHub...")
    r = requests.get(url)
    with open(zip_file, "wb") as f:
        f.write(r.content)

    # Make a directory to store the CGView.js code
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.mkdir(output_dir)

    # Extract the zip file
    print("Extracting CGView.js...")
    with zipfile.ZipFile(zip_file, "r") as zip_ref:
        zip_ref.extractall(output_dir)

    # Move the contents of the output_dir/cgview-js-main directory to the
    # output_dir, and delete the output_dir/cgview-js-main directory.
    for item in os.listdir(os.path.join(output_dir, "cgview-js-main")):
        shutil.move(os.path.join(output_dir, "cgview-js-main", item), output_dir)
    os.rmdir(os.path.join(output_dir, "cgview-js-main"))

    # Delete the zip file
    os.remove(zip_file)

    # Check whether the expected files are present
    check_files(output_dir)

    # Check whether the LICENSE file has the expected ending
    check_license_ending(output_dir)

    # Prefix the cgview.min.js file with the license (license_to_add).
    with open(os.path.join(output_dir, "docs", "dist", "cgview.min.js")) as f:
        cgview_js = f.read()
    with open(os.path.join(output_dir, "docs", "dist", "cgview.min.js"), "w") as f:
        f.write(license_to_add + cgview_js)


def check_files(output_dir: str) -> None:
    for expected_file in expected_files:
        if not os.path.exists(os.path.join(output_dir, expected_file)):
            raise Exception(f"Missing expected file: {expected_file}")


def check_license_ending(output_dir: str) -> None:
    with open(os.path.join(output_dir, "LICENSE")) as f:
        license_file = f.read()
    if not license_file.endswith(license_ending):
        raise Exception(f"The LICENSE file has an unexpected ending:\n\n{license_file}")


# Define a function to update the CGView.js code if necessary.
def update_cgview_js_code(downloaded_code_dir: str, existing_code_dir: str) -> None:
    # Check whether any of the files have been modified since the last update
    # (by comparing the md5sums to those of the previous versions of the files
    # that are already in the src/proksee-batch/data/cgview-js_code directory).
    for expected_file in expected_files:
        downloaded_file = os.path.join(downloaded_code_dir, expected_file)
        existing_file = os.path.join(existing_code_dir, expected_file)
        # Skip directories.
        if not os.path.isdir(downloaded_file):
            if not os.path.exists(existing_file):
                raise Exception(f"Missing expected file: {existing_file}")
            if not are_files_same(downloaded_file, existing_file):
                print(f"Updating {existing_file}...")
                shutil.copyfile(downloaded_file, existing_file)
            else:
                print(f"{existing_file} is up to date.")


def file_md5(file_path: str) -> str:
    """
    Compute the MD5 checksum of a file specified by the file path.

    Args:
    file_path (str): The path to the file whose MD5 checksum is to be calculated.

    Returns:
    str: The MD5 checksum of the file.
    """
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def are_files_same(file_path1: str, file_path2: str) -> bool:
    """
    Check if two files are the same based on their MD5 checksums.

    Args:
    file_path1 (str): The path to the first file.
    file_path2 (str): The path to the second file.

    Returns:
    bool: True if the files are the same, False otherwise.
    """
    return file_md5(file_path1) == file_md5(file_path2)


if __name__ == "__main__":
    # Define the URL of the latest CGView.js release
    url = "https://github.com/stothard-group/cgview-js/archive/main.zip"

    # Define the path to the directory to contain the new code version.
    output_dir = os.path.join(
        os.path.dirname(os.path.dirname(os.path.realpath(__file__))),
        "src",
        "proksee_batch",
        "data",
        "cgview-js_code_temp",
    )

    # Download the latest CGView.js code from GitHub
    download_cgview_js(url, output_dir)

    # Define the path to the directory containing the existing code version.
    existing_code_dir = os.path.join(
        os.path.dirname(os.path.dirname(os.path.realpath(__file__))),
        "src",
        "proksee_batch",
        "data",
        "cgview-js_code",
    )

    # Update the CGView.js code if necessary.
    update_cgview_js_code(output_dir, existing_code_dir)

    # Delete the temporary directory
    shutil.rmtree(output_dir)
