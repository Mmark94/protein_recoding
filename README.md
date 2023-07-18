# protein recoding
This is a set of simple python scripts to recode proteins. This script can be used to change codons in the DNA of protein sequences.

The "function change_codon_start_Dict_score" takes as input a DNA sequence of a protein and a dictionary of codons to change in the form {"codon_to_change":"codon_to_change_with"}.
The script will try to find the first codon, the ATG and it will start the find and replace.
In the script, I use as an  example the GFP protein.
For example, if the input codons are {"TCT":"TCG", "GCC":"GCT"}, the script will change all the TCT codons with TCG ad all the GCC codons with GCT.

The script is designed mostly to make synonymous changes, so if after the change the aa sequence changes, it will give a warning. However, this script can be used to perform all kinds of codon changes, also not synonymous.

The script will print the number of codons changed.

This script also calculates the change in the codon usage in the protein. I used the codon usage for the yeast genome (S. cerevisiae).

This is the codon usage for yeast genome (how much a codon is used to code an amino acid in proteins, the score is in percentage (from 0 to 1)):
Dic_codon_usage = {'TTT': 0.59, 'TTC': 0.41, 'TTA': 0.28, 'TTG': 0.29, 'CTT': 0.13, 'CTC': 0.06, 'CTA': 0.14,
                        'CTG': 0.11, 'ATT': 0.46, 'ATC': 0.26, 'ATA': 0.27, 'GTT': 0.39, 'GTC': 0.21, 'GTA': 0.21,
                        'GTG': 0.19, 'TCT': 0.26, 'TCC': 0.16, 'TCA': 0.21, 'TCG': 0.1, 'AGT': 0.16, 'AGC': 0.11,
                        'CCT': 0.31, 'CCC': 0.15, 'CCA': 0.41, 'CCG': 0.12, 'ACT': 0.35, 'ACC': 0.22, 'ACA': 0.3,
                        'ACG': 0.13, 'GCT': 0.38, 'GCC': 0.22, 'GCA': 0.29, 'GCG': 0.11, 'TAT': 0.56, 'TAC': 0.44,
                        'CAT': 0.64, 'CAC': 0.36, 'CAA': 0.69, 'CAG': 0.31, 'AAT': 0.59, 'AAC': 0.41, 'AAA': 0.58,
                        'AAG': 0.42, 'GAT': 0.65, 'GAC': 0.35, 'GAA': 0.71, 'GAG': 0.29, 'TGT': 0.63, 'TGC': 0.37,
                        'CGT': 0.15, 'CGC': 0.06, 'CGA': 0.07, 'CGG': 0.04, 'AGA': 0.48, 'AGG': 0.21, 'GGT': 0.47,
                        'GGC': 0.19, 'GGA': 0.22, 'GGG': 0.12, 'ATG': 1.0, 'TGG': 1.0, 'TAA': 0.61, 'TAG': 0.09, 'TGA': 0.3}

The script will calculate a score that reflects the changes in the codon usage. The score is: (score_previous_aa - score_new_aa) * number_of_aa_changed
The script will also calculate and print the previous usage codon mean and the new one.


I also add a new feature, if you want to not modify the first x aa in the sequence, you just put this number as the third argument in this way:
DIC= {"ATG":"CCC", "ATT":"CTC", "AAA":"GGG"}
change_codon_start_Dict_score(GFP,DIC, 5)
This will not change the first 5 aa. By default this is 0, so you can omit this third argument.
                     
