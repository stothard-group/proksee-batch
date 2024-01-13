import json

import requests


def generate_proksee_link(json_data_file: str, output_file: str) -> None:
    """
    Generate a Proksee link from a JSON file.
    """
    url = "https://proksee.ca/api/v1/projects.json"

    json_data = None
    with open(json_data_file) as file:
        json_data = json.load(file)
    data = {"data": json.dumps(json_data), "origin": "proksee-batch"}

    try:
        response = requests.post(
            url, headers={"Content-Type": "application/json"}, json=data
        )
        response.raise_for_status()
        data = response.json()
        if data.get("status") == "success" and data.get("url"):
            with open(output_file, "w") as file:
                file.write(data["url"])
        else:
            print(f"Failed to create Proksee project: {data.get('error')}")
    except requests.RequestException as e:
        print(f"Error: {e}")
