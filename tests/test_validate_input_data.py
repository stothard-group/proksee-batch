import os
import tempfile
from importlib import resources

from proksee_batch.validate_input_data import check_vcf_ref_vs_alt_genotypes
from proksee_batch.validate_input_data import check_vcf_seq_ids


def test_check_vcf_seq_ids() -> None:
    """Test check_vcf_seq_ids function. It should take a path to a VCF file and
    a path to a GenBank file, and check that all the sequence IDs in the first
    column of the VCF file are contigs in the GenBank file. This test function
    should check that the function returns the expected result given the GenBank
    file `data/minimal_example_5.gbk` and some simple temporary VCF files.
    """
    # Define path to GenBank file
    with resources.path(
        "tests.data.genbank", "minimal_example_5.gbk"
    ) as path, tempfile.TemporaryDirectory() as temp_dir:
        genbank_file_path = str(path)
        assert os.path.exists(genbank_file_path)

        # First test: VCF file with correct sequence IDs
        vcf_file_path = os.path.join(temp_dir, "test.vcf")
        with open(vcf_file_path, "w") as vcf_file:
            vcf_file.write(
                """##fileformat=VCFv4.2
##fileDate=20190815
##source=proksee_batch
#CHROM	POS	ID	REF	ALT
XXXX0000001.1	1	.	A	C
XXXX0000002.1	1	.	A	C
XXXX0000003.1	1	.	A	C
"""
            )
        assert check_vcf_seq_ids(vcf_file_path, genbank_file_path) == True

        # Second test: VCF file with incorrect sequence IDs
        vcf_file_path = os.path.join(temp_dir, "test.vcf")
        with open(vcf_file_path, "w") as vcf_file:
            vcf_file.write(
                """##fileformat=VCFv4.2
##fileDate=20190815
##source=proksee_batch
#CHROM	POS	ID	REF	ALT
XXXX0000004.1	1	.	A	C
XXXX0000005.1	1	.	A	C
XXXX0000006.1	1	.	A	C
"""
            )
        assert check_vcf_seq_ids(vcf_file_path, genbank_file_path) == False


def test_check_vcf_ref_vs_alt_genotypes() -> None:
    """Test check_vcf_ref_vs_alt_genotypes function. It should take a path to a
    VCF file and a path to a GenBank file. It should check that the genotypes in
    the genome in the GenBank file match the REF genotypes in the VCF file.
    """
    # Define path to GenBank file
    with resources.path(
        "tests.data.genbank", "minimal_example_5.gbk"
    ) as path, tempfile.TemporaryDirectory() as temp_dir:
        genbank_file_path = str(path)
        assert os.path.exists(genbank_file_path)

        # First test: VCF file with correct genotypes
        vcf_file_path = os.path.join(temp_dir, "test.vcf")
        with open(vcf_file_path, "w") as vcf_file:
            vcf_file.write(
                """##fileformat=VCFv4.2
##fileDate=20190815
##source=proksee_batch
#CHROM	POS	ID	REF	ALT
XXXX0000001.1	1	.	G	C
XXXX0000002.1	1	.	A	T
XXXX0000003.1	1	.	C	G
"""
            )
        assert (
            check_vcf_ref_vs_alt_genotypes(vcf_file_path, genbank_file_path, "genbank")
            == True
        )

        # Second test: VCF file with incorrect genotypes
        vcf_file_path = os.path.join(temp_dir, "test.vcf")
        with open(vcf_file_path, "w") as vcf_file:
            vcf_file.write(
                """##fileformat=VCFv4.2
##fileDate=20190815
##source=proksee_batch
#CHROM	POS	ID	REF	ALT
XXXX0000001.1	1	.	C	G
XXXX0000002.1	1	.	T	A
XXXX0000003.1	1	.	G	C
"""
            )
        assert (
            check_vcf_ref_vs_alt_genotypes(vcf_file_path, genbank_file_path, "genbank")
            == False
        )
