from collections import defaultdict
from pprint import pprint
import re
import pickle

from utils import *

genome = "TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"

skew_map = {
    "C": -1,
    "G": 1,
    "A": 0,
    "T": 0
}


def main():
    skew_array = [0]

    genome = read_file("dataset_7_6.txt")

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

    print(min_array)


if __name__ == '__main__':
    main()