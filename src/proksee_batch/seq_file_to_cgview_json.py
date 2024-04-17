import datetime
import json
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


def genbank_to_cgview_json(genome_name: str, genbank_file: str, json_file: str) -> None:
    """Convert a GenBank file to a CGView JSON file."""
    genbank_records: List[SeqRecord] = list(SeqIO.parse(genbank_file, "genbank"))  # type: ignore

    now: datetime.datetime = datetime.datetime.now()

    map_id: str = "".join(random.choices(string.ascii_lowercase + string.digits, k=40))

    json_data: Dict[str, Any] = {
        "cgview": {
            "version": "1.0",
            "created": now.strftime("%Y-%m-%d %H:%M:%S"),
            "id": map_id,
            "name": genome_name,
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
        # Determine contig name based on availability of attributes.
        genbank_contig_name = ""
        if record.id is not None and record.id:
            genbank_contig_name = record.id
        elif record.name is not None and record.name:
            genbank_contig_name = record.name

        assert (
            genbank_contig_name not in all_contig_names
        ), f"Duplicate contig name {genbank_contig_name} in {genbank_file}."
        all_contig_names.append(genbank_contig_name)

        contig_name = remove_problematic_characters_from_contig_name(
            genbank_contig_name
        )

        contig_data = seq_to_json_contig(contig_name, str(record.seq))

        json_data["cgview"]["sequence"]["contigs"].append(contig_data)

        features_to_skip = ["source", "gene", "exon"]

        for feature in record.features:
            if feature.type in features_to_skip:
                continue

            if feature.location is None:
                raise ValueError("Location is None")

            # If the feature is a CDS, determine the reading frame.
            feature_codon_start = 1
            if feature.type == "CDS":
                codon_start_qualifier = feature.qualifiers.get("codon_start")
                if codon_start_qualifier is not None and codon_start_qualifier:
                    feature_codon_start = int(codon_start_qualifier[0])
            assert feature_codon_start in [
                1,
                2,
                3,
            ], f"codon_start must be 1, 2, or 3, not {feature_codon_start}"

            # Determine the location type, and process accordingly
            locations = (
                feature.location.parts
                if isinstance(feature.location, CompoundLocation)
                else (
                    [feature.location]
                    if isinstance(feature.location, FeatureLocation)
                    else None
                )
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

            if feature.type == "CDS":
                if strand == 1:
                    # Assume that the locations are already sorted by start position from 5' to 3' on the forward strand.
                    pass
                elif strand == -1:
                    # Sort the locations by start position from 3' to 5' on the reverse strand by reversing the list.
                    locations = list(reversed(locations))

            # Initiate a counter for the number of base pairs from the start of
            # the first codon (for use in the case of CDS features).
            cumulative_bp_from_first_codon_start = 0

            # Iterate over the locations within this feature.
            loc_num = 0
            for loc in locations:
                loc_num += 1

                # Switch back to 1-based indexing (Biopython converts to 0-based indexing when parsing GenBank format).
                start = int(loc.start) + 1
                stop = int(loc.end)

                # Assign name based on availability of attributes.
                name = ""
                product = feature.qualifiers.get("product")
                locus_tag = feature.qualifiers.get("locus_tag")
                gene = feature.qualifiers.get("gene")

                if gene is not None and gene:
                    name = gene[0]
                elif locus_tag is not None and locus_tag:
                    name = locus_tag[0]
                elif product is not None and product:
                    name = product[0]

                codon_start = 1
                if feature.type == "CDS":
                    # Figure out what the reading frame is for this location
                    # within the CDS feature.
                    if loc_num == 1:
                        # The codon start for this location is the same as the
                        # codon start for the original CDS feature.
                        codon_start = feature_codon_start

                        # Start counting at the first bp of the first codon.
                        cumulative_bp_from_first_codon_start = (
                            abs(stop - start) + 1 - (feature_codon_start - 1)
                        )

                    else:
                        remainder_after_last_full_codon = (
                            cumulative_bp_from_first_codon_start % 3
                        )
                        if remainder_after_last_full_codon == 0:
                            codon_start = 1
                        elif remainder_after_last_full_codon == 1:
                            codon_start = 3
                        elif remainder_after_last_full_codon == 2:
                            codon_start = 2

                        # Add the length of the current location to the total.
                        cumulative_bp_from_first_codon_start += abs(stop - start) + 1

                feature_data = {
                    "type": feature.type,
                    "name": name,
                    "start": start,
                    "stop": stop,
                    "strand": strand,
                    "source": "genbank-features",
                    "contig": contig_name,
                    "legend": feature.type,
                }

                if feature.type == "CDS":
                    # Define the reading frame in the feature data using the JSON feature "codonStart".
                    feature_data["codonStart"] = codon_start

                # Add metadata.
                if (
                    (gene is not None and gene)
                    or (locus_tag is not None and locus_tag)
                    or (product is not None and product)
                ):
                    feature_data["meta"] = {}
                if locus_tag is not None and locus_tag:
                    feature_data["meta"]["locus_tag"] = locus_tag[0]
                if gene is not None and gene:
                    feature_data["meta"]["gene"] = gene[0]
                if product is not None and product:
                    feature_data["meta"]["product"] = product[0]

                json_data["cgview"]["features"].append(feature_data)

    with open(json_file, "w") as outfile:
        json.dump(json_data, outfile, indent=4)


def fasta_to_cgview_json(genome_name: str, fasta_file: str, json_file: str) -> None:
    """Convert a FASTA file to a CGView JSON file. The JSON file will be in the
    same format as generated by the `genbank_to_cgview_json` function. There
    will be no features in the JSON file, only sequences/contigs.
    """
    fasta_records: List[SeqRecord] = list(SeqIO.parse(fasta_file, "fasta"))  # type: ignore

    now: datetime.datetime = datetime.datetime.now()

    map_id: str = "".join(random.choices(string.ascii_lowercase + string.digits, k=40))

    json_data: Dict[str, Any] = {
        "cgview": {
            "version": "1.0",
            "created": now.strftime("%Y-%m-%d %H:%M:%S"),
            "id": map_id,
            "name": genome_name,
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
            "tracks": [],
        }
    }

    all_contig_names: List[str] = []

    contigs_json = []
    for record in fasta_records:
        # Determine contig name based on availability of attributes.
        fasta_contig_name = ""
        if record.id is not None and record.id:
            fasta_contig_name = record.id
        elif record.description is not None and record.description:
            fasta_contig_name = record.description

        assert (
            fasta_contig_name not in all_contig_names
        ), f"Duplicate contig name {fasta_contig_name} in {fasta_file}."
        all_contig_names.append(fasta_contig_name)

        contig_name = remove_problematic_characters_from_contig_name(fasta_contig_name)

        contigs_json.append(seq_to_json_contig(contig_name, str(record.seq)))

    json_data["cgview"]["sequence"]["contigs"] = contigs_json

    with open(json_file, "w") as json_fh:
        json.dump(json_data, json_fh, indent=4)


def seq_to_json_contig(seq_id: str, seq: str) -> Dict[str, Any]:
    """Convert a sequence to a dictionary with the sequence ID, sequence length,
    and sequence. The dictionary will be in the format expected for a sequence
    in a CGView JSON file.
    """
    # Define a dictionary to hold the sequence information.
    seq_dict = {}

    # Add the sequence ID to the dictionary.
    seq_dict["name"] = seq_id

    # Add the sequence length to the dictionary.
    seq_dict["length"] = str(len(seq))

    # Add the sequence to the dictionary.
    seq_dict["seq"] = seq

    # Add the sequence orientation to the dictionary.
    seq_dict["orientation"] = orientation

    # Return the dictionary.
    return seq_dict


def remove_problematic_characters_from_contig_name(contig_name: str) -> str:
    """Remove problematic characters from a contig name. This is necessary
    some downstream software may not be able to handle contig names with
    certain characters.
    """
    return contig_name.replace("|", "_").replace(";", "_")
