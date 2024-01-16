import datetime
import json
import os
import random
import string
from typing import Any
from typing import Dict
from typing import List

from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord


# Global variables (which should later be defined in a config file).
genetic_code: str = "11"
orientation: str = "+"


def genbank_to_cgview_json(genbank_file: str, json_file: str) -> None:
    genbank_records: List[SeqRecord] = list(SeqIO.parse(genbank_file, "genbank"))  # type: ignore

    now: datetime.datetime = datetime.datetime.now()

    map_id: str = "".join(random.choices(string.ascii_lowercase + string.digits, k=40))

    json_data: Dict[str, Any] = {
        "cgview": {
            "version": "1.0",
            "created": now.strftime("%Y-%m-%d %H:%M:%S"),
            "id": map_id,
            "name": os.path.basename(genbank_file.rsplit(".", 1)[0]),
            "geneticCode": genetic_code,
            "settings": {},
            "backbone": {},
            "ruler": {},
            "dividers": {},
            "annotation": {},
            "sequence": {"contigs": []},
            "captions": [],
            "legend": {"items": []},
            "features": [],
            "tracks": [
                {
                    "name": "Features",
                    "separateFeaturesBy": "strand",
                    "position": "both",
                    "dataType": "feature",
                    "dataMethod": "source",
                    "dataKeys": "genbank-features",
                }
            ],
        }
    }

    all_contig_names: List[str] = []
    for record in genbank_records:
        # Handle potentially problematic characters in contig names.
        assert record.name is not None
        assert (
            record.name not in all_contig_names
        ), f"Duplicate contig name {record.name} in {genbank_file}."
        all_contig_names.append(record.name)
        contig_name = record.name.replace("|", "_").replace(";", "_")

        contig_data: Dict[str, Any] = {
            "name": contig_name,
            "orientation": orientation,
            "length": len(record.seq),
            "seq": str(record.seq),
        }
        json_data["cgview"]["sequence"]["contigs"].append(contig_data)

        # Define feature types to skip.
        features_to_skip = ["source", "gene", "exon"]

        for feature in record.features:
            if feature.type in features_to_skip:
                continue

            if feature.location is None:
                raise ValueError("Location is None")

            # Determine the location type and process accordingly
            locations = (
                feature.location.parts
                if isinstance(feature.location, CompoundLocation)
                else [feature.location]
                if isinstance(feature.location, FeatureLocation)
                else None
            )

            # Determine strand.
            strand = 0
            if {-1} == {loc.strand for loc in locations}:
                strand = -1
            elif {1} == {loc.strand for loc in locations}:
                strand = 1
            else:
                raise ValueError(
                    "Strand is not consistent among locations within feature {}".format(
                        feature
                    )
                )

            # Determine whether any features start or end at the origin.
            spans_origin = False
            for loc in locations:
                if (
                    loc.start == 0
                    or loc.start == len(record.seq)
                    or loc.end == 0
                    or loc.end == len(record.seq)
                ):
                    spans_origin = True
                    break

            for loc in locations:
                start = (
                    int(loc.start) + 1
                )  # Biopython converts to 0-based indexing, but CGView uses 1-based indexing (like GenBank format).
                stop = int(loc.end)

                # Re-order start and end positions if necessary to conform to
                # the CGView convention.
                if strand == -1:
                    start, stop = stop, start

                # If the feature is a CDS and spans the origin (has locations
                # that start or stop there), then process it accordingly.
                if feature.type == "CDS" and spans_origin:
                    if strand == 1:
                        # Sum the length of all preceding locations.
                        preceding_length = 0
                        for preceding_loc_start, preceding_loc_stop in [
                            (loc.start + 1, loc.end) for loc in locations
                        ]:
                            if preceding_loc_stop < start:
                                preceding_length += (
                                    abs(preceding_loc_stop - preceding_loc_start) + 1
                                )
                        # Determine if the location ends at the origin.
                        if stop == len(record.seq):
                            # Determine if the length of the location is a multiple of three.
                            if (
                                abs(
                                    (preceding_length + stop)
                                    - (preceding_length + start)
                                )
                                + 1
                            ) % 3 != 0:
                                # Reduce the value of the stop position until the length of the location is a multiple of three.
                                while (
                                    abs(
                                        (preceding_length + stop)
                                        - (preceding_length + start)
                                    )
                                    + 1
                                ) % 3 != 0:
                                    stop -= 1
                        # Determine if the location starts at the origin.
                        if start == 1:
                            # Determine if the length of the location is a multiple of three.
                            if (
                                abs(
                                    (preceding_length + stop)
                                    - (preceding_length + start)
                                )
                                + 1
                            ) % 3 != 0:
                                # Increase the value of the start position until the length of the location is a multiple of three.
                                while (
                                    abs(
                                        (preceding_length + stop)
                                        - (preceding_length + start)
                                    )
                                    + 1
                                ) % 3 != 0:
                                    start += 1
                    elif strand == -1:
                        # Sum the length of all preceding locations.
                        preceding_length = 0
                        for preceding_loc_start, preceding_loc_stop in [
                            (loc.end, loc.start + 1) for loc in locations
                        ]:
                            if preceding_loc_stop > start:
                                preceding_length += (
                                    abs(preceding_loc_stop - preceding_loc_start) + 1
                                )
                        # Determine if the location starts at the origin.
                        if start == len(record.seq):
                            # Determine if the length of the location is a multiple of three.
                            if (
                                abs(
                                    (preceding_length + stop)
                                    - (preceding_length + start)
                                )
                                + 1
                            ) % 3 != 0:
                                # Reduce the value of the start position until the length of the location is a multiple of three.
                                while (
                                    abs(
                                        (preceding_length + stop)
                                        - (preceding_length + start)
                                    )
                                    + 1
                                ) % 3 != 0:
                                    start -= 1
                        # Determine if the location ends at the origin.
                        if stop == 1:
                            # Determine if the length of the location is a multiple of three.
                            if (
                                abs(
                                    (preceding_length + stop)
                                    - (preceding_length + start)
                                )
                                + 1
                            ) % 3 != 0:
                                # Increase the value of the stop position until the length of the location is a multiple of three.
                                while (
                                    abs(
                                        (preceding_length + stop)
                                        - (preceding_length + start)
                                    )
                                    + 1
                                ) % 3 != 0:
                                    stop += 1

                # Assign name based on availability of attributes.
                name = ""
                locus_tag = feature.qualifiers.get("locus_tag")
                gene = feature.qualifiers.get("gene")
                product = feature.qualifiers.get("product")

                if locus_tag is not None and locus_tag:
                    name = locus_tag[0]
                elif gene is not None and gene:
                    name = gene[0]
                elif product is not None and product:
                    name = product[0]

                feature_data = {
                    "type": feature.type,
                    "name": name,
                    "start": start,
                    "stop": stop,
                    "strand": strand,
                    "source": "genbank-features",
                    "contig": record.name,
                    "legend": feature.type,
                }
                json_data["cgview"]["features"].append(feature_data)

    with open(json_file, "w") as outfile:
        json.dump(json_data, outfile, indent=4)
