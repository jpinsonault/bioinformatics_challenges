from collections import defaultdict
from pprint import pprint
import re
import pickle
import operator
from itertools import permutations
import sys
from suffix_tree import SuffixTree

from utils import *

genome = "GCGCTCTCGTGTCGGCTCCGCGGCTCTCGTTCGCGCGCGCGCGTTCTCGTGCGTCGTCCGTCCGTCGTCGCGTCCGTCGTGTTCTCTCGCGCTCGTGTTCTCTCGTTCCGGTTCGTTCCGTCGTGTTCGCGCTCGCGCGTCGGCTCGTTCTCCGTCTCTCTCTCGTGTGCGCTCTCGTCGTCTCTCTCGCCGTCCGGTGCTCGCTCGCGTCGTCTCTCTCGT"
k = 8
d = 2

def main():
    kmers = find_kmers(k, genome)
    fuzzy_kmers = defaultdict(int)
    num_kmers = len(kmers)

    # suffix_tree = SuffixTree(genome)

    print(genome[:10])
    with progress_bar(num_kmers, "Searching") as progress:

        for count, search_string in enumerate(kmers.keys()):
            progress.update(count)

            permuted_strings = mutations(search_string, d)
        
            for string in permuted_strings:
                fuzzy_kmers[string] = len(kmers[search_string])
            


    sorted_kmers = sorted(fuzzy_kmers.iteritems(), key=operator.itemgetter(1))
    pprint(sorted_kmers[-10:])


if __name__ == '__main__':
    main()