from collections import defaultdict
from pprint import pprint
import re
import pickle
import operator
from itertools import permutations


from utils import *

genome = "CACGTACAACACGTCGTCACGTCGTCGTCGTCACACACAACACGTACATTACATTCACGTCACGTTTTTCATTCACACGTTTCACATTCATTCACGTTTACATTTTTTTTACATTTTCACGTCAACATTCACGTACATTTTCACACATTCAACACGTTTCACGTCACACACGTCATTCACAACATTCATTCGTCAACACATTCGTACACGTCACACACACGTACACACAACACACGTCACACA"
k = 9
d = 2

def main():
    kmers = find_kmers(k, genome)
    fuzzy_kmers = defaultdict(int)
    num_kmers = len(kmers)

    for count, search_string in enumerate(kmers.keys()):
        if count % (num_kmers/20) == 0:
            print("{:.2f}%".format(count/float(num_kmers)*100))

        permuted_strings = mutations(search_string, d)
        # permuted_strings += mutations(reverse_compliment(search_string), d)
        # permuted_strings = permutations(search_string)
        # num_permutations = len(permuted_strings)
        # print(num_permutations)

        for pcount, permutation in enumerate(permuted_strings):
            # if pcount % 1000 == 0:
                # print("{:.2f}".format(pcount/float(num_permutations)*100))
            if permutation not in fuzzy_kmers:
                # print(permutation)
                reversed_permuation = reverse_compliment(permutation)
                reverse_not_equal = reversed_permuation != permutation

                for index in range(len(genome) - k + 1):
                    word = genome[index:index + k]
                    miss_count = 0
                    reversed_miss_count = 0
                    for x, letter in enumerate(word):
                        if letter != permutation[x]:
                            miss_count += 1

                        if letter != reversed_permuation[x]:
                            reversed_miss_count += 1

                    if miss_count <= d:
                        fuzzy_kmers[permutation] += 1

                    if reverse_not_equal:
                        if reversed_miss_count <= d:
                            fuzzy_kmers[permutation] += 1


    sorted_kmers = sorted(fuzzy_kmers.iteritems(), key=operator.itemgetter(1))
    pprint(sorted_kmers[-10:])


if __name__ == '__main__':
    main()