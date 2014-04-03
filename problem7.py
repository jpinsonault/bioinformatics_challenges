from collections import defaultdict
from pprint import pprint
import re
import pickle
import operator
from itertools import permutations
from timeit import timeit


from utils import *

genome = "CACAGTAGGCGCCGGCACACACAGCCCCGGGCCCCGGGCCGCCCCGGGCCGGCGGCCGCCGGCGCCGGCACACCGGCACAGCCGTACCGGCACAGTAGTACCGGCCGGCCGGCACACCGGCACACCGGGTACACACCGGGGCGCACACACAGGCGGGCGCCGGGCCCCGGGCCGTACCGGGCCGCCGGCGGCCCACAGGCGCCGGCACAGTACCGGCACACACAGTAGCCCACACACAGGCGGGCGGTAGCCGGCGCACACACACACAGTAGGCGCACAGCCGCCCACACACACCGGCCGGCCGGCACAGGCGGGCGGGCGCACACACACCGGCACAGTAGTAGGCGGCCGGCGCACAGCC"
k = 10
d = 2

def main():
    kmers = find_kmers(k, genome)
    counts = defaultdict(int)
    num_kmers = len(genome) - k + 1

    with progress_bar(num_kmers, "Searching") as progress:

        for index in range(num_kmers):
            word = genome[index:index + k]
            progress.update(index)

            permuted_strings = mutations(word, d)

            for permutation in permuted_strings:
                counts[permutation] += 1

    sorted_kmers = sorted(counts.iteritems(), key=lambda x: x[1])
    pprint(sorted_kmers[-10:])

def old_algorithm():
    kmers = find_kmers(k, genome)
    counts = defaultdict(int)
    num_kmers = len(kmers)


    with progress_bar(num_kmers, "Searching") as progress:
        for count, search_string in enumerate(kmers.keys()):
            progress.update(count)

            permuted_strings = mutations(search_string, d)
            # permuted_strings = permutations(search_string)
            # num_permutations = len(permuted_strings)
            # print(num_permutations)

            for pcount, permutation in enumerate(permuted_strings):
                # if pcount % 1000 == 0:
                    # print("{:.2f}".format(pcount/float(num_permutations)*100))
                if permutation not in counts:
                    # print(permutation)

                    for index in range(len(genome) - k + 1):
                        word = genome[index:index + k]
                        miss_count = 0
                        for x, letter in enumerate(word):
                            if letter != permutation[x]:
                                miss_count += 1
                                if miss_count > d:
                                    break

                        if miss_count <= d:
                            counts[permutation] += 1


    sorted_kmers = sorted(counts.iteritems(), key=operator.itemgetter(1))
    pprint(sorted_kmers[-10:])


if __name__ == '__main__':
    print("First: {:.2f} seconds".format(timeit(main, number=1)))
    print("Second: {:.2f} seconds".format(timeit(old_algorithm, number=1)))
