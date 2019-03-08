# Grading

+ Python script found in script/draw_motifs.py
+ SVG output file found in output/Fig_1

# motif-mark
Python script to plot protein binding motifs on an image of an exon and flanking introns

+ Handles cassette splicing
+ Input motif file must have individual motifs on new lines

# Minimum requirements:

+ Python3 compatible code
+ Use argparse
+ Useful help message
+ Input FASTA file and motifs filePreview the document
+ Multiple sequences
+ Multiple motifs
+ Ambiguous motif handling (only Y and N)
	+ Y --> Pyrimidines (C or T)
	+ N --> Unknown nucleic acid residue (A, T, C, or G)
+ svg output
+ Key/labeling

# Test files:

INSR.fasta

+ Exons capitalized, introns lowercase

motifs.txt

+ Contains one motif sequence of interest

# General organization

+ data contains all input data and temp fasta created from script
+ output contains all output svg
+ resources contains lecture materials and demos
+ script contains final script ran (in jupyter notebook form or .py)
