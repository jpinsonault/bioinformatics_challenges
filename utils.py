from collections import defaultdict
from pprint import pprint
import re


def find_kmers(k, genome):
    positions = defaultdict(list)

    for index in range(len(genome) - k):
        word = genome[index:index + k]
        positions[word].append(index)

    return positions


def read_file(filename):
    with open(filename) as o_file:
        return o_file.readlines()[0].strip()