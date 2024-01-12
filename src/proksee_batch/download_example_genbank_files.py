"""Downloads GenBank files from NCBI FTP site."""
import gzip
import os
import shutil

import requests


def download_file(url: str, local_filename: str) -> None:
    """
    Downloads a file from a given URL and saves it to the local file system.

    Parameters:
    url (str): The URL of the file to be downloaded.
    local_filename (str): The local path, including filename, where the file should be saved.

    Returns:
    None
    """
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_filename, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)


def download_example_genbank_files(genome_dir: str) -> None:
    """
    Downloads GenBank files from NCBI FTP site.

    Parameters:
    genome_dir (str): The path to the directory where the GenBank files should be saved.

    Returns:
    None
    """
    urls = [
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gbff.gz",
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/007/632/255/GCF_007632255.1_ASM763225v1/GCF_007632255.1_ASM763225v1_genomic.gbff.gz",
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.gbff.gz",
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/925/GCF_000006925.2_ASM692v2/GCF_000006925.2_ASM692v2_genomic.gbff.gz",
    ]

    # If the directory exists, delete it and create a new one.
    if os.path.exists(genome_dir):
        shutil.rmtree(genome_dir)
    os.mkdir(genome_dir)

    for url in urls:
        local_filename = os.path.join(genome_dir, url.split("/")[-1])
        download_file(url, local_filename)
        # Decompress the file.
        with gzip.open(local_filename, "rb") as f_in:
            with open(local_filename[:-3], "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        # Delete the compressed file.
        os.remove(local_filename)


if __name__ == "__main__":
    import sys

    # Download GenBank files.
    download_example_genbank_files(sys.argv[1])
