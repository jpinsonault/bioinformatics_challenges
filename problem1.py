from collections import defaultdict
from pprint import pprint
genome = "GACATATTTACGCTGAGGAATACCATGCCGAGGAATTTACGCTACCATGCCCCGGATCTTTACGCTGACATATTTACGCTGACATATGAGGAATTTACGCTCCGGATCTTTACGCTACCATGCCCCGGATCTTTACGCTGACATATGACATATACCATGCCACCATGCCTTACGCTCCGGATCTTTACGCTGAGGAATACCATGCCACCATGCCGAGGAATTTACGCTACCATGCCGACATATCCGGATCTTTACGCTGACATATCCGGATCTACCATGCCACCATGCCCCGGATCTGACATATTTACGCTGACATATTTACGCTTTACGCTGACATATACCATGCCGAGGAATTTACGCTGACATATGAGGAATTTACGCTACCATGCCTTACGCTGACATATGACATATTTACGCTCCGGATCTTTACGCTACCATGCCACCATGCCGACATATCCGGATCTGACATATTTACGCTGAGGAATGACATATACCATGCCCCGGATCTTTACGCTCCGGATCTGACATATCCGGATCTGACATATGAGGAATTTACGCTGACATATGACATATGAGGAATTTACGCTCCGGATCTCCGGATCTGAGGAATGACATATACCATGCCACCATGCCTTACGCTTTACGCTCCGGATCTTTACGCTACCATGCCGACATATACCATGCCTTACGCTGACATATGAGGAATGACATATACCATGCCGAGGAATTTACGCTCCGGATCTGAGGAATTTACGCTCCGGATCTGACATATTTACGCTGACATATCCGGATCTGACATATTTACGCTGAGGAATACCATGCCGAGGAAT"
k = 14


def main():
    counts = defaultdict(int)

    for index in range(len(genome) - k):
        counts[genome[index:index + k]] += 1

    pprint(counts)

    sorted_counts = sorted(list(counts.items()), key=lambda item: item[1])

    pprint(sorted_counts)

if __name__ == '__main__':
    main()