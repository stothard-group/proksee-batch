def update_legend(cgview_map_json_data: dict) -> dict:
    """
    Update the legend of the CGView map so that the colors of tracks for
    additional features are consistent among maps of different assemblies (with
    may have different numbers of tracks and feature types).
    """
    relevant_datakey_prefixes = ["gff", "blast", "bed", "vcf"]

    # Initialize the dictionary to hold additional feature tracks.
    additional_feature_tracks_by_type = {}

    # Iterate over each track and group additional feature tracks by type.
    track_num = 0
    for track in cgview_map_json_data["cgview"]["tracks"]:
        if track["name"] in ["CG Content", "CG Skew", "Features"]:
            continue
        track_num += 1
        for prefix in relevant_datakey_prefixes:
            if track["dataKeys"].startswith(prefix + "_"):
                # Use setdefault to initialize as a list if not already and append the track info.
                additional_feature_tracks_by_type.setdefault(prefix, []).append(
                    {"name": track["name"], "track_num": track_num}
                )
                break

    # Collect the names of the tracks in order they appear.
    track_names = [
        track["name"]
        for track in cgview_map_json_data["cgview"]["tracks"]
        if track["name"] not in ["CG Content", "CG Skew", "Features"]
    ]

    # Ensure that the 'legend' key is initialized as a list in the JSON data.
    assert "legend" in cgview_map_json_data["cgview"]

    # Update the legend, alternating colors between black and blue.
    color = "rgba(0,0,0,1)"
    for track in track_names:
        cgview_map_json_data["cgview"]["legend"]["items"].append(
            {"name": track, "swatchColor": color, "decoration": "arc"}
        )
        # Switch color for the next track.
        color = "rgba(91,91,91,1)" if color == "rgba(0,0,0,1)" else "rgba(0,0,0,1)"

    return cgview_map_json_data
