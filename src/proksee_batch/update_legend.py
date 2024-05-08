from typing import Any
from typing import Dict

import seaborn as sns  # type: ignore


def update_legend(cgview_map_json_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Update the legend of the CGView map so that the colors of tracks for
    additional features are consistent among maps of different assemblies (with
    may have different numbers of tracks and feature types).
    """
    # Define the prefixes of data keys that correspond to additional feature
    # tracks (from GFF files, etc.).
    relevant_datakey_prefixes = ["gff", "blast", "bed", "vcf"]

    # Define the names of legend items that should not be modified (if present).
    names_of_legend_items_to_leave_unmodified = [
        "CG Content",
        "CG Skew+",
        "CG Skew-",
    ]

    # Define the common GenBank feature types that are typically present in
    # the maps, and should be assigned consistent colors.
    common_genbank_features = [
        "CDS",
        "tRNA",
        "rRNA",
        "tmRNA",
        "ncRNA",
    ]

    # Ensure that the 'legend' key is initialized as a list in the JSON data.
    assert "legend" in cgview_map_json_data["cgview"]

    # Extract the items in the legend that should not be modified. These will be
    # appended to the new legend later (if necessary).
    legend_items_to_leave_unmodified = [
        item
        for item in cgview_map_json_data["cgview"]["legend"]["items"]
        if item["name"] in names_of_legend_items_to_leave_unmodified
    ]

    # Clear the legend items in the JSON data.
    cgview_map_json_data["cgview"]["legend"]["items"] = []

    # Iterate over all the features, and compile a comprehensive nonredundant
    # list of all legend names referenced by the features.
    assert "features" in cgview_map_json_data["cgview"], "No features found in map."
    all_legend_names = set()
    for feature in cgview_map_json_data["cgview"]["features"]:
        all_legend_names.add(feature["legend"])

    # Initialize the dictionary to hold additional feature track/legend names.
    additional_feature_tracks_by_type: Dict[str, Any] = {}

    # Iterate over each track/legend and group additional feature tracks by type.
    for track in cgview_map_json_data["cgview"]["tracks"]:
        for prefix in relevant_datakey_prefixes:
            if track["dataKeys"].startswith(prefix + "_"):
                # Use setdefault to initialize as a list if not already and append the track info.
                additional_feature_tracks_by_type.setdefault(prefix, []).append(
                    track["name"]
                )
                break

    # Get a list of all the additional feature track names (joined values from the dict).
    all_additional_feature_track_names = set().union(
        *additional_feature_tracks_by_type.values()
    )

    # The track names for these additional features should be among the legend
    # item names.
    for track_name in all_additional_feature_track_names:
        assert (
            track_name in all_legend_names
        ), f"Track name '{track_name}' is not among the legend item names."

    # The names_of_legend_items_to_leave_unmodified should not be among the
    # legend item names.
    for name in names_of_legend_items_to_leave_unmodified:
        assert (
            name not in all_legend_names
        ), f"Name '{name}' is among the legend item names."

    # Count the number of legend items that need to be assigned a color.
    num_legends_to_color = len(all_legend_names)

    # Set a default number of colors to use in a palette.
    default_num_colors = 12

    # Obtain a color palette of an appropriate length from Seaborn
    # (https://seaborn.pydata.org/tutorial/color_palettes.html).
    palette = sns.color_palette("husl", max(num_legends_to_color, default_num_colors))

    ## Check that none of the legend names contain a comma.
    # for legend_name in all_legend_names:
    #    assert (
    #        "," not in legend_name
    #    ), f"Legend name '{legend_name}' contains a comma, which is not allowed."

    # Sort the legend names for additional feature tracks alphabetically, assign
    # them colors from the palette (starting from the beginning), and add
    # corresponding legend items to the legend.
    # Legend items should be of the form: {"name": "legend name here", "swatchColor": "rgba(0,0,0,1)", "decoration": "arc"}
    for legend_name in sorted(all_additional_feature_track_names):
        color = palette.pop(0)
        cgview_map_json_data["cgview"]["legend"]["items"].append(
            {
                "name": legend_name,
                "swatchColor": f"rgba({int(color[0]*255)},{int(color[1]*255)},{int(color[2]*255)},1)",
                "decoration": "arc",
            }
        )

    # Reverse the remaining palette colors.
    palette.reverse()

    # Make a list of the common genbank features (legends) that are present in
    # the complete list.
    common_genbank_legend_names = [
        legend_name
        for legend_name in common_genbank_features
        if legend_name in all_legend_names
    ]

    # For each common_genbank_features, assign a color from the palette and add
    # corresponding legend items to the legend.
    for legend_name in common_genbank_legend_names:
        color = palette.pop(0)
        opacity = 1.0
        if legend_name == "CDS":
            opacity = 0.5
        cgview_map_json_data["cgview"]["legend"]["items"].append(
            {
                "name": legend_name,
                "swatchColor": f"rgba({int(color[0]*255)},{int(color[1]*255)},{int(color[2]*255)},{opacity})",
                "decoration": "arrow",
            }
        )

    # For the remaining legend items, assign a color from the palette and add
    # corresponding legend items to the legend.
    remaining_legends = list(
        set(all_legend_names)
        - set(common_genbank_legend_names)
        - set(all_additional_feature_track_names)
    )
    for legend_name in remaining_legends:
        color = palette.pop(0)
        cgview_map_json_data["cgview"]["legend"]["items"].append(
            {
                "name": legend_name,
                "swatchColor": f"rgba({int(color[0]*255)},{int(color[1]*255)},{int(color[2]*255)},1)",
                "decoration": "arc",
            }
        )

    # Append the legend items that should not be modified to the new legend.
    cgview_map_json_data["cgview"]["legend"]["items"].extend(
        legend_items_to_leave_unmodified
    )

    # Return the updated JSON data.
    return cgview_map_json_data
