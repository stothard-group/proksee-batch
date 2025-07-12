import json
import os
from typing import Any
from typing import Dict
from typing import List
from typing import Tuple
from typing import TypedDict
from typing import Union

import gffutils  # type: ignore

from .remove_covered_features import remove_covered_features_optimized


class FeatureDecorationDict(TypedDict):
    color: str
    opacity: float


class BlastMetaDict(TypedDict):
    query: str
    query_start: int
    query_stop: int
    subject: str
    subject_start: int
    subject_stop: int
    alignment_length: int
    identity: float
    mismatches: int
    evalue: float
    bit_score: float


class BlastFeatureDict(TypedDict):
    contig: str
    type: str
    start: int
    stop: int
    strand: Union[int, str]  # Can be "." or int
    source: str
    name: str
    legend: str
    tags: List[str]
    meta: BlastMetaDict


class BedMetaDict(TypedDict):
    score: float


class BedFeatureDict(TypedDict):
    contig: str
    type: str
    start: int
    stop: int
    strand: int
    source: str
    name: str
    legend: str
    tags: List[str]
    meta: BedMetaDict


class VcfMetaDict(TypedDict):
    ref: str
    alt: str


class VcfFeatureDict(TypedDict):
    contig: str
    type: str
    start: int
    stop: int
    strand: int
    source: str
    name: str
    legend: str
    tags: List[str]
    meta: VcfMetaDict


class GffMetaDict(TypedDict, total=False):
    score: float
    gene: str
    locus_tag: str
    product: str
    # Additional dynamic fields from GFF attributes


class GffFeatureDict(TypedDict):
    contig: str
    type: str
    start: int
    stop: int
    strand: Union[int, str]
    source: str
    name: str
    legend: str
    tags: List[str]
    meta: Dict[str, Any]  # GFF has dynamic meta fields


class TrackDict(TypedDict):
    name: str
    separateFeaturesBy: str
    position: str
    thicknessRatio: Union[int, float]
    dataType: str
    dataMethod: str
    dataKeys: str
    drawOrder: str


def parse_blast_files(
    blast_files: List[str],
) -> Tuple[List[BlastFeatureDict], List[TrackDict]]:
    """
    Parses BLAST result files.

    Parameters:
    blast_files (list): A list of paths to BLAST result files.

    Returns:
    tuple: A tuple containing a list of parsed BLAST features and a list of parsed BLAST tracks.
    """
    blast_features: List[BlastFeatureDict] = []
    blast_tracks: List[TrackDict] = []
    blast_file: str = ""  # Initialize to avoid unbound variable
    for i, blast_file in enumerate(blast_files):
        num = i + 1

        file_specific_blast_features: List[BlastFeatureDict] = []

        # Define the BLAST track.
        blast_track: TrackDict = {
            "name": os.path.basename(blast_file),
            "separateFeaturesBy": "none",
            "position": "both",
            "thicknessRatio": 1.0,
            "dataType": "feature",
            "dataMethod": "source",
            "dataKeys": f"blast_{num}",
            "drawOrder": "score",
        }
        # Add the BLAST track to the list of BLAST tracks.
        blast_tracks.append(blast_track)

        # Parse the BLAST result file.
        with open(blast_file) as f:
            blast_results = f.readlines()
        # Parse the BLAST result file.
        for blast_result_line in blast_results:
            # Parse the BLAST result line.
            blast_result = blast_result_line.strip().split("\t")

            # Skip the header line(s), or blank lines, if present.
            if (
                blast_result[0] == "qseqid"
                or blast_result[0].startswith("#")
                or blast_result[0].startswith("\n")
            ):
                continue

            # Define the BLAST feature.
            blast_feature: BlastFeatureDict = {
                "name": [
                    blast_result[1] if len(blast_result[1]) <= 20 else "blast_hit"
                ][
                    0
                ],  # "blast_hit" if the name is too long
                "type": "blast",
                "start": int(
                    # blast_result[8]
                    min([int(blast_result[6]), int(blast_result[7])])
                ),  # Here we assume that the subject sequence is a contig in the genome/map.
                "stop": int(
                    # blast_result[9]
                    max([int(blast_result[6]), int(blast_result[7])])
                ),  # Here we assume that the subject sequence is a contig in the genome/map.
                # "strand": 1 if blast_result[8] < blast_result[9] else -1,
                "strand": ".",
                "source": f"blast_{num}",
                "contig": blast_result[
                    0
                ],  # Here we assume that the subject sequence is a contig in the genome/map.
                "legend": os.path.basename(blast_file),
                "tags": [],
                "meta": {
                    "query": blast_result[0],
                    "query_start": int(blast_result[6]),
                    "query_stop": int(blast_result[7]),
                    "subject": blast_result[1],
                    "subject_start": int(blast_result[8]),
                    "subject_stop": int(blast_result[9]),
                    "alignment_length": int(blast_result[3]),
                    "identity": float(blast_result[2]),
                    "mismatches": int(blast_result[4]),
                    "evalue": float(blast_result[10]),
                    "bit_score": float(blast_result[11]),
                },
            }
            # Add the BLAST feature to the list of BLAST features.
            file_specific_blast_features.append(blast_feature)

        # Remove BLAST features that are completely covered by another BLAST feature
        # with an equal or higher score. Apply filtering per contig to avoid
        # features from different contigs competing with each other.

        # Group features by contig
        features_by_contig: Dict[str, List[BlastFeatureDict]] = {}
        for feature in file_specific_blast_features:
            contig = str(feature["contig"])
            if contig not in features_by_contig:
                features_by_contig[contig] = []
            features_by_contig[contig].append(feature)

        # Apply filtering per contig and collect results
        filtered_features: List[BlastFeatureDict] = []
        for contig, contig_features in features_by_contig.items():
            # Apply remove_covered_features to features from this contig only
            contig_filtered = [
                feature
                for feature, is_not_covered in zip(
                    contig_features,
                    remove_covered_features_optimized(
                        get_feature_locations_and_scores_from_blast_features(
                            contig_features
                        )
                    ),
                )
                if is_not_covered
            ]
            filtered_features.extend(contig_filtered)

        file_specific_blast_features = filtered_features

        # Add the file-specific BLAST features to the list of BLAST features.
        blast_features += file_specific_blast_features

    assert (
        len(blast_features) > 0 and len(blast_tracks) > 0
    ), "No BLAST features or tracks were obtained from input files. Please check that the query sequences are sequences (contigs, chromosomes, etc.) of the genome being mapped."
    return (blast_features, blast_tracks)


def get_feature_locations_and_scores_from_blast_features(
    blast_features: List[BlastFeatureDict],
) -> List[Tuple[int, int, float]]:
    """
    Gets feature locations and scores from BLAST features.

    Parameters:
    blast_features (list): A list of parsed BLAST features.

    Returns:
    list: A list of tuples containing feature locations and scores.
    """
    feature_locations_and_scores: List[Tuple[int, int, float]] = []
    for blast_feature in blast_features:
        feature_locations_and_scores.append(
            (
                blast_feature["start"],
                blast_feature["stop"],
                blast_feature["meta"]["bit_score"],
            )
        )
    return feature_locations_and_scores


def add_blast_features_and_tracks(
    blast_files: List[str], json_file: str, output_file: str
) -> None:
    """
    Parses BLAST result files, adds the parsed BLAST features and tracks to the cgview map JSON data structure, and writes the cgview map JSON data structure to a new file.

    Parameters:
    blast_files (list): A list of paths to BLAST result files.
    json_file (str): The path to a cgview map JSON file.
    output_file (str): The path to the output file.

    Returns:
    None
    """
    # Parse the BLAST result file.
    blast_features, blast_tracks = parse_blast_files(blast_files)
    # Read the cgview map JSON file.
    with open(json_file) as f:
        json_data = json.load(f)

    # Get list of contigs from JSON data.
    contigs = [contig["name"] for contig in json_data["cgview"]["sequence"]["contigs"]]

    # Check that all the blast features are assigned to contigs that are in the JSON data.
    for blast_feature in blast_features:
        assert (
            blast_feature["contig"] in contigs
        ), "The subject sequence {} listed in the BLAST result file {} is not among the contigs in the genome being mapped.".format(
            blast_feature["contig"], blast_feature["legend"]
        )

    # Add the parsed BLAST features and tracks to the cgview map JSON data structure.
    json_data["cgview"]["features"] += blast_features
    json_data["cgview"]["tracks"] += blast_tracks
    # Write the cgview map JSON data structure to a new file.
    with open(output_file, "w") as f:
        json.dump(json_data, f, indent=4)


def parse_bed_files(
    bed_files: List[str],
) -> Tuple[List[BedFeatureDict], List[TrackDict]]:
    """
    Parses BED files.

    Parameters:
    bed_files (list): A list of paths to BED files.

    Returns:
    tuple: A tuple containing a list of parsed BED features and a list of parsed BED tracks.
    """
    bed_features: List[BedFeatureDict] = []
    bed_tracks: List[TrackDict] = []
    bed_file: str = ""  # Initialize to avoid unbound variable
    for i, bed_file in enumerate(bed_files):
        num = i + 1

        file_specific_bed_features: List[BedFeatureDict] = []

        # Define the BED track.
        bed_track: TrackDict = {
            "name": os.path.basename(bed_file),
            "separateFeaturesBy": "none",
            "position": "both",
            "thicknessRatio": 1,
            "dataType": "feature",
            "dataMethod": "source",
            "dataKeys": f"bed_{num}",
            "drawOrder": "score",
        }
        # Add the BED track to the list of BED tracks.
        bed_tracks.append(bed_track)

        # Parse the BED file.
        with open(bed_file) as f:
            bed_results = f.readlines()
        # Parse the BED file.
        for bed_result_line in bed_results:
            # Parse the BED result line.
            bed_result = bed_result_line.strip().split()

            # Skip the header line(s), if present.
            if (
                bed_result[0].startswith("#")
                or bed_result[0].startswith("track ")
                or bed_result[0].startswith("browser ")
            ):
                continue

            assert (
                len(bed_result) >= 3
            ), f"BED file {bed_file} is not in the correct format."

            name = ""
            if len(bed_result) >= 4:
                name = bed_result[3]
            score = 0.0
            if len(bed_result) >= 5:
                score = float(bed_result[4])
            strand = 1
            if len(bed_result) >= 6:
                if bed_result[5] == "+":
                    strand = 1
                elif bed_result[5] == "-":
                    strand = -1
                else:
                    raise ValueError(
                        "The strand column of BED file {} is not in the correct format.".format(
                            bed_file
                        )
                    )

            # Define the BED feature.
            bed_feature: BedFeatureDict = {
                "name": name,
                "type": "bed",
                "start": int(bed_result[1]) + 1,
                "stop": int(bed_result[2]),
                "strand": strand,
                "source": f"bed_{num}",
                "contig": bed_result[0],
                "legend": os.path.basename(bed_file),
                "tags": [],
                "meta": {
                    "score": score,
                },
            }
            # Add the BED feature to the list of BED features.
            file_specific_bed_features.append(bed_feature)

        # Remove BED features that are completely covered by another BED feature
        # with an equal or higher score. Apply filtering per contig to avoid
        # features from different contigs competing with each other.

        # Group features by contig
        features_by_contig: Dict[str, List[BedFeatureDict]] = {}
        for feature in file_specific_bed_features:
            contig = str(feature["contig"])
            if contig not in features_by_contig:
                features_by_contig[contig] = []
            features_by_contig[contig].append(feature)

        # Apply filtering per contig and collect results
        filtered_features: List[BedFeatureDict] = []
        for contig, contig_features in features_by_contig.items():
            # Apply remove_covered_features to features from this contig only
            contig_filtered = [
                feature
                for feature, is_not_covered in zip(
                    contig_features,
                    remove_covered_features_optimized(
                        get_feature_locations_and_scores_from_bed_features(
                            contig_features
                        )
                    ),
                )
                if is_not_covered
            ]
            filtered_features.extend(contig_filtered)

        file_specific_bed_features = filtered_features

        # Add the file-specific BED features to the list of BED features.
        bed_features += file_specific_bed_features

    assert (
        len(bed_features) > 0 and len(bed_tracks) > 0
    ), "No BED features or tracks were obtained from input files."
    return (bed_features, bed_tracks)


def get_feature_locations_and_scores_from_bed_features(
    bed_features: List[BedFeatureDict],
) -> List[Tuple[int, int, float]]:
    """
    Gets feature locations and scores from BED features.

    Parameters:
    bed_features (list): A list of parsed BED features.

    Returns:
    list: A list of tuples containing feature locations and scores.
    """
    feature_locations_and_scores: List[Tuple[int, int, float]] = []
    for bed_feature in bed_features:
        feature_locations_and_scores.append(
            (
                bed_feature["start"],
                bed_feature["stop"],
                bed_feature["meta"]["score"],
            )
        )
    return feature_locations_and_scores


def add_bed_features_and_tracks(
    bed_files: List[str], json_file: str, output_file: str
) -> None:
    """
    Parses BED files, adds the parsed BED features and tracks to the cgview map JSON data structure, and writes the cgview map JSON data structure to a new file.

    Parameters:
    bed_files (list): A list of paths to BED files.
    json_file (str): The path to a cgview map JSON file.
    output_file (str): The path to the output file.

    Returns:
    None
    """
    # Parse the BED file.
    bed_features, bed_tracks = parse_bed_files(bed_files)
    # Read the cgview map JSON file.
    with open(json_file) as f:
        json_data = json.load(f)

    # Get list of contigs from JSON data.
    contigs = [contig["name"] for contig in json_data["cgview"]["sequence"]["contigs"]]

    # Check that all the bed features are assigned to contigs that are in the JSON data.
    for bed_feature in bed_features:
        assert (
            bed_feature["contig"] in contigs
        ), "The contig {} listed in the BED file {} is not among the contigs in the genome being mapped.".format(
            bed_feature["contig"], bed_feature["legend"]
        )

    # Add the parsed BED features and tracks to the cgview map JSON data structure.
    json_data["cgview"]["features"] += bed_features
    json_data["cgview"]["tracks"] += bed_tracks
    # Write the cgview map JSON data structure to a new file.
    with open(output_file, "w") as f:
        json.dump(json_data, f, indent=4)


def parse_vcf_files(
    vcf_files: List[str],
) -> Tuple[List[VcfFeatureDict], List[TrackDict]]:
    """
    Parses VCF files.

    Parameters:
    vcf_files (list): A list of paths to VCF files.

    Returns:
    tuple: A tuple containing a list of parsed VCF features and a list of parsed VCF tracks.
    """
    vcf_features: List[VcfFeatureDict] = []
    vcf_tracks: List[TrackDict] = []
    vcf_file: str = ""  # Initialize to avoid unbound variable
    for i, vcf_file in enumerate(vcf_files):
        num = i + 1

        file_specific_vcf_features: List[VcfFeatureDict] = []

        # Define the VCF track.
        vcf_track: TrackDict = {
            "name": os.path.basename(vcf_file),
            "separateFeaturesBy": "none",
            "position": "both",
            "thicknessRatio": 1,
            "dataType": "feature",
            "dataMethod": "source",
            "dataKeys": f"vcf_{num}",
            "drawOrder": "score",
        }
        # Add the VCF track to the list of VCF tracks.
        vcf_tracks.append(vcf_track)

        # Parse the VCF file.
        with open(vcf_file) as f:
            vcf_results = f.readlines()
        # Parse the VCF file.
        for vcf_result_line in vcf_results:
            # Parse the VCF result line.
            vcf_result = vcf_result_line.strip().split()

            # Skip the header line(s), if present.
            if vcf_result[0].startswith("#") or vcf_result[0].startswith("\n"):
                continue

            assert (
                len(vcf_result) >= 5
            ), f"VCF file {vcf_file} is not in the correct format."

            name = ""
            if len(vcf_result) >= 8:
                name = vcf_result[7][20:]

            # Define the VCF feature.
            vcf_feature: VcfFeatureDict = {
                "name": name,
                "type": "vcf",
                "start": int(vcf_result[1]),
                "stop": int(vcf_result[1])
                + len(vcf_result[3])
                - 1,  # Assumes that the genome being mapped is the reference genome.
                "strand": 1,
                "source": f"vcf_{num}",
                "contig": vcf_result[0],
                "legend": os.path.basename(vcf_file),
                "tags": [],
                "meta": {
                    "ref": vcf_result[3],
                    "alt": vcf_result[4],
                },
            }

            # Add the VCF feature to the list of VCF features.
            file_specific_vcf_features.append(vcf_feature)

        # Remove VCF features that are completely covered by another VCF feature
        # with an equal or higher score. Apply filtering per contig to avoid
        # features from different contigs competing with each other.

        # Group features by contig
        features_by_contig: Dict[str, List[VcfFeatureDict]] = {}
        for feature in file_specific_vcf_features:
            contig = str(feature["contig"])
            if contig not in features_by_contig:
                features_by_contig[contig] = []
            features_by_contig[contig].append(feature)

        # Apply filtering per contig and collect results
        filtered_features: List[VcfFeatureDict] = []
        for contig, contig_features in features_by_contig.items():
            # Apply remove_covered_features to features from this contig only
            contig_filtered = [
                feature
                for feature, is_not_covered in zip(
                    contig_features,
                    remove_covered_features_optimized(
                        get_feature_locations_and_scores_from_vcf_features(
                            contig_features
                        )
                    ),
                )
                if is_not_covered
            ]
            filtered_features.extend(contig_filtered)

        file_specific_vcf_features = filtered_features

        # Add the file-specific VCF features to the list of VCF features.
        vcf_features += file_specific_vcf_features

    assert (
        len(vcf_features) > 0 and len(vcf_tracks) > 0
    ), "No VCF features or tracks were obtained from input files."
    return (vcf_features, vcf_tracks)


def get_feature_locations_and_scores_from_vcf_features(
    vcf_features: List[VcfFeatureDict],
) -> List[Tuple[int, int, float]]:
    """
    Gets feature locations and scores from VCF features.

    Parameters:
    vcf_features (list): A list of parsed VCF features.

    Returns:
    list: A list of tuples containing feature locations and scores.
    """
    feature_locations_and_scores: List[Tuple[int, int, float]] = []
    for vcf_feature in vcf_features:
        feature_locations_and_scores.append(
            (
                vcf_feature["start"],
                vcf_feature["stop"],
                1.0,
            )
        )
    return feature_locations_and_scores


def add_vcf_features_and_tracks(
    vcf_files: List[str], json_file: str, output_file: str
) -> None:
    """
    Parses VCF files, adds the parsed VCF features and tracks to the cgview map JSON data structure, and writes the cgview map JSON data structure to a new file.

    Parameters:
    vcf_files (list): A list of paths to VCF files.
    json_file (str): The path to a cgview map JSON file.
    output_file (str): The path to the output file.

    Returns:
    None
    """
    # Parse the VCF file.
    vcf_features, vcf_tracks = parse_vcf_files(vcf_files)
    # Read the cgview map JSON file.
    with open(json_file) as f:
        json_data = json.load(f)

    # Get list of contigs from JSON data.
    contigs = [contig["name"] for contig in json_data["cgview"]["sequence"]["contigs"]]

    # Check that all the vcf features are assigned to contigs that are in the JSON data.
    for vcf_feature in vcf_features:
        assert (
            vcf_feature["contig"] in contigs
        ), "The contig {} listed in the VCF file {} is not among the contigs in the genome being mapped.".format(
            vcf_feature["contig"], vcf_feature["legend"]
        )

    # Add the parsed VCF features and tracks to the cgview map JSON data structure.
    json_data["cgview"]["features"] += vcf_features
    json_data["cgview"]["tracks"] += vcf_tracks
    # Write the cgview map JSON data structure to a new file.
    with open(output_file, "w") as f:
        json.dump(json_data, f, indent=4)


# Use gffutils to parse GFF files.
def parse_gff_files(
    gff_files: List[str],
) -> Tuple[List[GffFeatureDict], List[TrackDict]]:
    """
    Parses GFF files.

    Parameters:
    gff_files (list): A list of paths to GFF files.

    Returns:
    tuple: A tuple containing a list of parsed GFF features and a list of parsed GFF tracks.
    """
    gff_features: List[GffFeatureDict] = []
    gff_tracks: List[TrackDict] = []
    for i, gff_file in enumerate(gff_files):
        num = i + 1

        file_specific_gff_features: List[GffFeatureDict] = []

        # Define the GFF track.
        gff_track: TrackDict = {
            "name": os.path.basename(gff_file),
            "separateFeaturesBy": "none",
            "position": "both",
            "thicknessRatio": 1,
            "dataType": "feature",
            "dataMethod": "source",
            "dataKeys": f"gff_{num}",
            "drawOrder": "score",
        }
        # Add the GFF track to the list of GFF tracks.
        gff_tracks.append(gff_track)

        # Parse the GFF file.
        db = gffutils.create_db(  # type: ignore[attr-defined]
            gff_file,
            ":memory:",
            force=True,
            keep_order=True,
            merge_strategy="merge",
            sort_attribute_values=True,
        )

        # Define feature types to exclude.
        feature_types_to_exclude = ["gene", "exon", "region"]

        # Parse the GFF file.
        for gff_result in db.all_features():  # type: ignore[attr-defined]
            # Skip the header line(s), if present, and feature types to exclude.
            if gff_result.featuretype in feature_types_to_exclude:  # type: ignore[attr-defined]
                continue

            # Define the name based on availability.
            name = gff_result.id  # type: ignore[attr-defined]
            if "gene" in gff_result.attributes:  # type: ignore[attr-defined]
                name = gff_result.attributes["gene"][0]  # type: ignore[attr-defined]

            # Define the GFF feature.
            gff_feature: GffFeatureDict = {
                "name": name,
                "type": "gff",
                "start": gff_result.start,  # type: ignore[attr-defined]
                "stop": gff_result.stop,  # type: ignore[attr-defined]
                "strand": gff_result.strand,  # type: ignore[attr-defined]
                "source": f"gff_{num}",
                "contig": gff_result.seqid,  # type: ignore[attr-defined]
                "legend": os.path.basename(gff_file),
                "tags": [],
                "meta": {
                    "score": 0.0,
                },
            }

            # Add any additional metadata, if present.
            for attribute, values in gff_result.attributes.items():  # type: ignore[attr-defined]
                if values:  # Checks if the list is non-empty
                    gff_feature["meta"][attribute] = (
                        float(values[0]) if attribute == "score" else values[0]
                    )

            # Add the GFF feature to the list of GFF features.
            file_specific_gff_features.append(gff_feature)

        # Close the GFF database.
        db.conn.close()  # type: ignore[attr-defined]

        assert (
            len(file_specific_gff_features) > 0
        ), f"No GFF features or tracks were obtained from input file {gff_file}."

        # Remove GFF features that are completely covered by another GFF feature
        # with an equal or higher score. Apply filtering per contig to avoid
        # features from different contigs competing with each other.

        # Group features by contig
        features_by_contig: Dict[str, List[GffFeatureDict]] = {}
        for feature in file_specific_gff_features:
            contig = str(feature["contig"])
            if contig not in features_by_contig:
                features_by_contig[contig] = []
            features_by_contig[contig].append(feature)

        # Apply filtering per contig and collect results
        filtered_features: List[GffFeatureDict] = []
        for contig, contig_features in features_by_contig.items():
            # Apply remove_covered_features to features from this contig only
            contig_filtered = [
                feature
                for feature, is_not_covered in zip(
                    contig_features,
                    remove_covered_features_optimized(
                        get_feature_locations_and_scores_from_gff_features(
                            contig_features
                        )
                    ),
                )
                if is_not_covered
            ]
            filtered_features.extend(contig_filtered)

        file_specific_gff_features = filtered_features

        assert (
            len(file_specific_gff_features) > 0
        ), f"No GFF features or tracks were obtained from input file {gff_file}."

        # Add the file-specific GFF features to the list of GFF features.
        gff_features += file_specific_gff_features

    return (gff_features, gff_tracks)


def get_feature_locations_and_scores_from_gff_features(
    gff_features: List[GffFeatureDict],
) -> List[Tuple[int, int, float]]:
    """
    Gets feature locations and scores from GFF features.

    Parameters:
    gff_features (list): A list of parsed GFF features.

    Returns:
    list: A list of tuples containing feature locations and scores.
    """
    feature_locations_and_scores: List[Tuple[int, int, float]] = []
    for gff_feature in gff_features:
        feature_locations_and_scores.append(
            (
                gff_feature["start"],
                gff_feature["stop"],
                gff_feature["meta"]["score"],
            )
        )
    return feature_locations_and_scores


def add_gff_features_and_tracks(
    gff_files: List[str], json_file: str, output_file: str
) -> None:
    """
    Parses GFF files, adds the parsed GFF features and tracks to the cgview map JSON data structure, and writes the cgview map JSON data structure to a new file.

    Parameters:
    gff_files (list): A list of paths to GFF files.
    json_file (str): The path to a cgview map JSON file.
    output_file (str): The path to the output file.

    Returns:
    None
    """
    # Parse the GFF file.
    gff_features, gff_tracks = parse_gff_files(gff_files)
    # Read the cgview map JSON file.
    with open(json_file) as f:
        json_data = json.load(f)

    # Get list of contigs from JSON data.
    contigs = [contig["name"] for contig in json_data["cgview"]["sequence"]["contigs"]]

    # Check that all the gff features are assigned to contigs that are in the JSON data.
    for gff_feature in gff_features:
        assert (
            gff_feature["contig"] in contigs
        ), "The contig {} listed in the GFF file {} is not among the contigs in the genome being mapped.".format(
            gff_feature["contig"], gff_feature["legend"]
        )

    # Add the parsed GFF features and tracks to the cgview map JSON data structure.
    for gff_feature in gff_features:
        assert gff_feature not in json_data["cgview"]["features"]
    json_data["cgview"]["features"] += gff_features
    json_data["cgview"]["tracks"] += gff_tracks
    # Write the cgview map JSON data structure to a new file.
    with open(output_file, "w") as f:
        json.dump(json_data, f, indent=4)
