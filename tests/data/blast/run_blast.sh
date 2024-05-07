#blastn -query tests/data/fasta/U49845.1.fna -subject tests/data/fasta/U49845.1.fna -outfmt 6 > tests/data/blast/U49845.1.txt

# Run this to generate the blast output files for tests.
blastn -query tests/data/fasta/U49845.1.fna -subject tests/data/blast/subject1.fna -outfmt 6 > tests/data/blast/U49845.1_subject1.txt
blastn -query tests/data/fasta/U49845.1.fna -subject tests/data/blast/subject2.fna -outfmt 6 > tests/data/blast/U49845.1_subject2.txt
