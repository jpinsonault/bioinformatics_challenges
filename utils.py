from collections import defaultdict
from pprint import pprint
import re
import itertools


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
