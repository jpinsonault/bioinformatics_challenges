from collections import defaultdict
from pprint import pprint
import re
import pickle

from utils import *

genome = "AACAAGCTGATAAACATTTAAAGAG"
search_string = "AAAAA"
d = 2

skew_map = {
    "C": -1,
    "G": 1,
    "A": 0,
    "T": 0
}


def main():
    positions = []

    search_len = len(search_string)

    for index in range(len(genome) - search_len + 1):
        word = genome[index:index + search_len]
        miss_count = 0
        for x, letter in enumerate(word):
            if letter != search_string[x]:
                miss_count += 1

        if miss_count <= d:
            positions.append(index)

    print(" ".join([str(position) for position in positions]))
    print(len(positions))


if __name__ == '__main__':
    main()