"""Downloads GenBank files from NCBI FTP site."""
import gzip
import os
import shutil
from importlib import resources

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


def download_genbank_file(output_dir: str, url: str) -> None:
    """
    Downloads GenBank files from NCBI FTP site.

    Parameters:
    output_dir (str): The path to the directory where the GenBank files should be saved.
    url (str): The URL of the file to be downloaded.

    Returns:
    None
    """
    assert os.path.exists(output_dir)
    assert os.path.isdir(output_dir)

    local_filename = os.path.join(output_dir, url.split("/")[-1])
    download_file(url, local_filename)
    # Decompress the file.
    with gzip.open(local_filename, "rb") as f_in:
        with open(local_filename[:-3], "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    # Delete the compressed file.
    os.remove(local_filename)


example_data_dict = {
    "Escherichia_coli_str._K-12_substr._MG1655": {
        "genbank": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gbff.gz",
        "gff": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz",
        "fasta": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz",
    },
    "Klebsiella_aerogenes": {
        "genbank": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/007/632/255/GCF_007632255.1_ASM763225v1/GCF_007632255.1_ASM763225v1_genomic.gbff.gz",
        "gff": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/007/632/255/GCF_007632255.1_ASM763225v1/GCF_007632255.1_ASM763225v1_genomic.gff.gz",
        "fasta": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/007/632/255/GCF_007632255.1_ASM763225v1/GCF_007632255.1_ASM763225v1_genomic.fna.gz",
    },
    "Salmonella_enterica_subsp._enterica_serovar_Typhimurium_str._LT2": {
        "genbank": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.gbff.gz",
        "gff": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.gff.gz",
        "fasta": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna.gz",
    },
    "Shigella_flexneri_2a_str._301": {
        "genbank": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/925/GCF_000006925.2_ASM692v2/GCF_000006925.2_ASM692v2_genomic.gbff.gz",
        "gff": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/925/GCF_000006925.2_ASM692v2/GCF_000006925.2_ASM692v2_genomic.gff.gz",
        "fasta": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/925/GCF_000006925.2_ASM692v2/GCF_000006925.2_ASM692v2_genomic.fna.gz",
    },
}


def download_example_input(output_dir: str) -> None:
    """
    Downloads example GenBank, BLAST, and BED files for several example bacterial genomes.
    """
    # If the provided directory exists, delete it.
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    # Retrieve the path to the local example input data directory using
    # importlib.
    with resources.path("proksee_batch.data", "example_input_dir") as template_path:
        template_path_str = str(template_path)
        assert os.path.exists(template_path_str)
        # Copy the example input data directory from the tests/data directory to
        # the output directory path.
        shutil.copytree(template_path_str, output_dir)

    # Check that the subdirectories in the output directory match the primary
    # keys in the example_data_dict.
    assert set(example_data_dict.keys()) == set(os.listdir(output_dir))

    # Download the example GenBank, GFF, and FASTA files using the URLs in the example_data_dict.
    for genome_name, genome_data in example_data_dict.items():
        output_subdir = os.path.join(output_dir, genome_name, "genbank")
        os.mkdir(output_subdir)
        download_genbank_file(output_subdir, genome_data["genbank"])
        output_subdir = os.path.join(output_dir, genome_name, "gff")
        os.mkdir(output_subdir)
        download_genbank_file(output_subdir, genome_data["gff"])
        output_subdir = os.path.join(output_dir, genome_name, "fasta")
        os.mkdir(output_subdir)
        download_genbank_file(output_subdir, genome_data["fasta"])


if __name__ == "__main__":
    import sys

    # Download GenBank files.
    download_example_genbank_files(sys.argv[1])
