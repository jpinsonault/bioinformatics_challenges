from collections import defaultdict
from pprint import pprint
import re
import operator
from time import time
import sys

from utils import *

rna_string = "ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA"
tyrocidine = "VKLFPWFNQY"

peptide_string = "MA"


def main():
    # genome = read_fasta("Salmonella_enterica.fasta")

    # with progress_bar(len(rna_string)) as progress:

    # print(translate_rna(rna_string))

    rna_encoders = []

    len_peptide = len(peptide_string)

    word_size = len_peptide * 3

    for i in xrange(0, len(rna_string) - word_size):
        rna_substring = rna_string[i:i + word_size]

        if rna_encodes_peptide(rna_substring, peptide_string):
            rna_encoders.append(rna_substring)

        if rna_encodes_peptide(reverse_compliment(rna_substring), peptide_string):
            rna_encoders.append(reverse_compliment(rna_substring))

    pprint(rna_encoders)


if __name__ == '__main__':
    main()