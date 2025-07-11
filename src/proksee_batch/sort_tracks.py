from typing import Any
from typing import Dict
from typing import List


def sort_tracks(cgview_map_json_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Update the track list in the CGView map so that the tracks for additional
    features from GFF files, etc. are sorted alphabetically.
    """
    # Define the prefixes of data keys that correspond to additional feature
    # tracks (from GFF files, etc.).
    relevant_datakey_prefixes = ["gff", "blast", "bed", "vcf"]

    # Ensure that the 'tracks' key is initialized as a list in the JSON data.
    assert "tracks" in cgview_map_json_data["cgview"]

    # Initialize a list to hold the tracks for additional features.
    additional_feature_tracks: List[Dict[str, Any]] = []

    # Iterate over each track and add it to the list if it's name starts with
    # one of the relevant prefixes.
    for track in cgview_map_json_data["cgview"]["tracks"]:
        for prefix in relevant_datakey_prefixes:
            if track["dataKeys"].startswith(prefix + "_"):
                additional_feature_tracks.append(track)
                break

    # Sort the additional feature tracks alphabetically by name.
    additional_feature_tracks.sort(key=lambda x: str(x["name"]))

    # Make the complete list of tracks by concatenating the sorted additional
    # feature tracks with the tracks that are not additional features.
    cgview_map_json_data["cgview"]["tracks"] = [
        track
        for track in cgview_map_json_data["cgview"]["tracks"]
        if track not in additional_feature_tracks
    ] + additional_feature_tracks

    # Return the updated JSON data.
    return cgview_map_json_data
