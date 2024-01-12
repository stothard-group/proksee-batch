import json
from typing import Any
from typing import Dict


def merge_cgview_json_with_template(
    basic_json_file: str, template_file: str, output_file: str
) -> None:
    """
    Merge a basic cgview map in JSON format with a Proksee configuration file in JSON format.
    """
    # Load the basic cgview map
    with open(basic_json_file) as basic_json:
        basic_cgview_map = json.load(basic_json)

    # Load the Proksee configuration file
    with open(template_file) as template_json:
        proksee_configuration = json.load(template_json)

    # Define keys for formatting and data components
    formatting_keys = [
        "version",
        "created",
        "updated",
        "settings",
        "backbone",
        "ruler",
        "annotation",
        "dividers",
        "highlighter",
        "captions",
        "legend",
    ]
    data_keys = ["name", "sequence", "features", "plots", "bookmarks", "tracks"]

    # Create a new JSON structure
    merged_data: Dict[str, Any] = {"cgview": {}}

    # Add formatting components
    for key in formatting_keys:
        if key in proksee_configuration["cgview"]:
            merged_data["cgview"][key] = proksee_configuration["cgview"][key]

    # Add data components
    for key in data_keys:
        if key in basic_cgview_map["cgview"]:
            merged_data["cgview"][key] = basic_cgview_map["cgview"][key]

    # Assert that there is only one element in the "captions" list.
    assert (
        len(merged_data["cgview"]["captions"]) == 1
    ), "The 'captions' list must contain exactly one element."
    # Change the value of the "name" key in the first element of the "captions"
    # list to the value of the "name" key in the main "cgview" dictionary.
    merged_data["cgview"]["captions"][0]["name"] = merged_data["cgview"]["name"]

    # Write the merged data to a new file
    with open(output_file, "w") as file:
        json.dump(merged_data, file, indent=4)
