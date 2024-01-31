Bioinformatics Python Functions
This project contains Python functions for common bioinformatics tasks.

Functions
is_valid_dna(seq) - Checks if a DNA sequence contains only A, T, C, G nucleotides
get_gc_content(seq) - Calculates GC content percentage of a DNA sequence
reverse_compliment(seq) - Returns reverse compliment of a DNA sequence
change_t_to_u(seq) - Replaces T with U in an RNA sequence
get_transcription(seq, strand) - Returns RNA transcription of a DNA sequence on given strand
stats_amino_acids() - Returns count of each amino acid from genetic code dict
translate_seq(seq, reading_frame) - Translates DNA sequence into amino acids
count_codons(seq) - Counts occurrences of each codon in DNA sequence
Usage
Import the functions and call them by passing DNA/RNA sequences as strings:

from class_project_1 import get_gc_content, translate_seq

my_seq = "ACGTTGCA"
gc_content = get_gc_content(my_seq)

protein = translate_seq(my_seq)



Genetic Code
A dictionary gencode is defined with the standard genetic code for mapping codons to amino acids. This is used by translate_seq() and stats_amino_acids().

Requirements
Python 3.x
No external packages required