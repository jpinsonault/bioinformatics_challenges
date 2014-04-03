from collections import defaultdict
from pprint import pprint
import re
import pickle

from utils import *

genome = "AACAAGCTGATAAACATTTAAAGAG"
search_string = "AAAAA"
d = 1


def main():
    matches = set()
    positions = []

    search_len = len(search_string)

    kmers = find_kmers(search_len, genome)

    permutations = list(mutations(search_string, d))

    # with progress_bar(len(permutations)) as progress:
    #     for count, permutation in enumerate(permutations):
    #         progress.update(count)
    #         for position in kmers[permutation]:
    #             matches.add(position)

    # print(matches)


    for index in range(len(genome) - search_len + 1):
        word = genome[index:index + search_len]
        miss_count = 0
        for x, letter in enumerate(word):
            if letter != search_string[x]:
                miss_count += 1

        if miss_count <= d:
            positions.append(index)

    # print(" ".join(sorted([str(position) for position in positions])))
    print(len(positions))


if __name__ == '__main__':
    main()