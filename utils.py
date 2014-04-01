from collections import defaultdict
from pprint import pprint
import re
import itertools
from itertools import izip
from contextlib import contextmanager

from progressbar import Percentage
from progressbar import Bar
from progressbar import ProgressBar

codon_table_rna = {
    "AAA": "K", "AAC": "N", "AAG": "K", "AAU": "N",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACU": "T",
    "AGA": "R", "AGC": "S", "AGG": "R", "AGU": "S",
    "AUA": "I", "AUC": "I", "AUG": "M", "AUU": "I",
    "CAA": "Q", "CAC": "H", "CAG": "Q", "CAU": "H",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCU": "P",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGU": "R",
    "CUA": "L", "CUC": "L", "CUG": "L", "CUU": "L",
    "GAA": "E", "GAC": "D", "GAG": "E", "GAU": "D",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCU": "A",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGU": "G",
    "GUA": "V", "GUC": "V", "GUG": "V", "GUU": "V",
    "UAA": "",  "UAC": "Y", "UAG": "",  "UAU": "Y",
    "UCA": "S", "UCC": "S", "UCG": "S", "UCU": "S",
    "UGA": "",  "UGC": "C", "UGG": "W", "UGU": "C",
    "UUA": "L", "UUC": "F", "UUG": "L", "UUU": "F",
}

codon_table = {
    "AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AGA": "R", "AGC": "S", "AGG": "R", "AGT": "S",
    "ATA": "I", "ATC": "I", "ATG": "M", "ATT": "I",
    "CAA": "Q", "CAC": "H", "CAG": "Q", "CAT": "H",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
    "GAA": "E", "GAC": "D", "GAG": "E", "GAT": "D",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "TAA": "",  "TAC": "Y", "TAG": "",  "TAT": "Y",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TGA": "",  "TGC": "C", "TGG": "W", "TGT": "C",
    "TTA": "L", "TTC": "F", "TTG": "L", "TTT": "F",
}

def find_kmers(k, genome):
    positions = defaultdict(list)

    for index in range(len(genome) - k):
        word = genome[index:index + k]
        positions[word].append(index)

    return positions


def find_min_skew(genome):
    skew_map = {
        "C": -1,
        "G": 1,
        "A": 0,
        "T": 0
    }

    skew_array = [0]

    for x, letter in enumerate(genome):
        skew_array.append(skew_map[letter] + skew_array[x])

    minimum = 0
    min_array = []
    for x, skew in enumerate(skew_array):
        if skew == minimum:
            min_array.append(x)
        elif skew < minimum:
            minimum = skew
            min_array = [x]

    return min_array


def read_file(filename):
    with open(filename) as in_file:
        return in_file.readlines()[0].strip()


def read_fasta(filename):
    genome = ""

    print("Reading {}".format(filename))

    with open(filename) as in_file:
        for line in in_file:
            if line.startswith(">"):
                continue

            genome += line.strip()

    return genome


def hamming_distance(s1, s2):
    #Return the Hamming distance between equal-length sequences
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def mutations(word, hamming_distance, charset='ATCG'):
    for indices in itertools.combinations(range(len(word)), hamming_distance):
        for replacements in itertools.product(charset, repeat=hamming_distance):
            mutation = list(word)
            for index, replacement in zip(indices, replacements):
                mutation[index] = replacement
            yield "".join(mutation)


def reverse_compliment(string):
    letter_map = {
        "A": "T",
        "G": "C",
        "T": "A",
        "C": "G",
    }

    return "".join(reversed([letter_map[letter] for letter in string]))


def reverse_compliment_rna(string):
    letter_map = {
        "A": "U",
        "G": "C",
        "U": "A",
        "C": "G",
    }

    return "".join(reversed([letter_map[letter] for letter in string]))


@contextmanager
def progress_bar(max_number, title=""):
    progress_bar = ProgressBar(widgets=[title + ": " , Percentage(), Bar()], maxval=max_number).start()
    yield progress_bar

    progress_bar.finish()


def translate_rna(string):
    if len(string) % 3 != 0:
        raise TypeError("input string should be a multiple of 3")

    translated = []

    for i in xrange(0, len(string), 3):
        codon = string[i:i + 3]
        translated.append(codon_table[codon])

    return "".join(translated)

def count_rna_in_codons(codon_string):
    counts = defaultdict(int)

    for rna, codon in codon_table.iteritems():
        counts[codon] += 1

    count = 1
    for letter in codon_string:
        count *= counts[letter]
    print(count)


def rna_encodes_peptide(dna_string, peptide_string):
    if len(dna_string) % 3 != 0:
        return False

    encodings = [codon_table[dna_string[i:i + 3]] for i in xrange(0, len(dna_string), 3)]

    print(dna_string, "".join(encodings))
    for encoding, peptide in izip(encodings, peptide_string):

        if encoding != peptide:
            return False

    return True
