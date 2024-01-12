from typing import Tuple

from Bio import SeqIO


def get_stats_from_genbank(genbank_file: str) -> Tuple[str, int, int, float]:
    """
    Get basic stats from a GenBank file.
    """
    description = ""
    total_size = 0
    number_of_contigs = 0
    gc_content = 0.0

    for record in SeqIO.parse(genbank_file, "genbank"):  # type: ignore
        # Use the description from the first record (contig) in the file.
        description = str(record.description)
        break

    for record in SeqIO.parse(genbank_file, "genbank"):  # type: ignore
        # Process each record (contig) in the file
        total_size += len(record.seq)
        number_of_contigs += 1
        gc_content += record.seq.count("G")
        gc_content += record.seq.count("C")

    gc_content = round(gc_content / total_size, 2)

    return (description, total_size, number_of_contigs, gc_content)
