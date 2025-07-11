from typing import Tuple

from Bio import SeqIO  # type: ignore
from Bio.SeqRecord import SeqRecord  # type: ignore


def get_stats_from_seq_file(
    seq_file: str, format: str
) -> Tuple[str, str, int, int, float]:
    """
    Get basic stats from a GenBank or FASTA file.
    """
    assert format in ["genbank", "fasta"]

    accession: str = ""
    description: str = ""
    total_size: int = 0
    number_of_contigs: int = 0
    gc_count: int = 0

    # Get description from first record
    records = list(SeqIO.parse(seq_file, format))  # type: ignore[arg-type, no-untyped-call]
    if records:
        first_record = records[0]
        description = str(first_record.description) if first_record.description else ""

    # Process all records
    for i, record in enumerate(records):
        # record is already a SeqRecord from the list
        if i == 0:
            accession = str(record.id) if record.id else ""
        # Process each record (contig) in the file
        if record.seq is not None:
            seq_str = str(record.seq)
            total_size += len(seq_str)
            gc_count += seq_str.count("G")
            gc_count += seq_str.count("C")
        number_of_contigs += 1

    gc_content: float = round(gc_count / total_size, 4) if total_size > 0 else 0.0

    return (accession, description, total_size, number_of_contigs, gc_content)
