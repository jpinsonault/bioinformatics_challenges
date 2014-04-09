from collections import defaultdict
from pprint import pprint
import re
import operator
from time import time
import sys

from utils import *

rna_string = "TCCTTCCCTCGATTATACAAGTCTGGCTAGCCGACTTGCTCATGATTTCTCCATTCTGCAGGGCGCAAGAACCTTCTAGACGGTGCTGTACGTGCACAACGCCGCAGAGGGGGGCGGGGAGAAGATCCGTTTGGCAGCAGCACCATTACGAGTGTCATTATCCATGTTAAGCTGTTCAACGGTTCATTGAGCGATCGGGTCTGTCTGATCGACCAATAGGTCCCGGCAGTATCCATGTGCTAGGAAGACTGCACGTCAGATGAGGAAAAAATTCAATACGGTTGGTAATCGTTGGAGCCTTACGTTGTGGAACGCGTCATTCGAAGCAGGCCAATGGTGAAATCATGTCAAAGTCTCCGTTGTGCGAGAGAATTCGGCAAGATTAGCCCGTGTGTATGAGACCTGTACATAGCAGCCTTATTAACCGCTAAGTCAGTCAACTGGCACGTGCGCGTAGATGAGAACACCTGGTAAATACCGCGGTTACGCACGTTCGTACTCGGAGTGGTAGCGCCGGTATTTGAACCCGAGCTCAATCATGCTGTTGGCGGGACAAAAAGTACATGCCTTTATGTGGGGCAGGGCTTTGAATCTTATCCATGACATATTCGCCGAATCGCCCCGCTTAACTTTCCGATCGAGGTGTTGAGCTCAGGGCATAGATGAAAGTTACTATCGTCGACTTTGAACTCCCCCGCTAACCGCAACTTTTCTACCGCTTGAGCTCATAGGAATCGAGAGGAAACGGCGAAATTATGTCCAAGTCGTCGCCGTTACCCCGCAGATAGTGACAATATAATGATAAAGCCTTGTTTGACTAAGGCTGGGATACATTGTCTATTCGCATCAGAATGATCGAATTATCGCAAGTTGACATCTCGAACCATGTAACCAGTAGATGGCGGAGATATTTGTCCGCAATCTTCTACGTTAAACTCAAAAATCTTCTGTACAAAGATAGTTCTTCAAGGCCGGGGTGAAGGTACGATTGAAGTGGAACGGTGCGTATCATGGGTCTATTAACCCTACTGCCTAAAGCAGAATCCCGAAACATCCGTTAATAAGCCCCTTTAGAGGAACGGGTTCTATACAAGGCCCAAATCTACGGCATTACCACGCAACTCGACGCCAGGGGGAACCCGGTTTAGAGGACGTAGCGACACGATGTAGACTTCATGACGCCCAGTTATTTGAACCATTTCAATACGGAGTGTAGAATGATCAGGTGCACAAGAATTAGCGTAGGAAAATGTTGTTACACCCTACCCGGGGAGCGATCTGGCGTCCCAGGTCACATAAATCGGTCGGTCACTACACAGCCCGAGATAACTGGCGGGCTGCACATCGCGCAGCTGCCAGCCTAGTACTCAGTCGGTGGGATGCCCCCCTCTGATCAATTCTAAGTCCAGATAATATGGAGACGCTATATTACCCCTGATCAGCGATGAATTCTCAATTTTAAAGATATCCAAGCAGTAGCTTAAAACATAAATTCCTTCCGTTGATCCCATTTTCTCAAGTTTGGTCGCCGACGAGGTGGCGGATAAAACCCCGCGATTAAGGCTGCCTTGAACCCAGTGACTGGGGGGTATCCCGCGGCCGTCCATTTGTGCTGTCATGGCGGTACACTAGTGGTGCCACCTCAAACTCACACAATACCTTCTGGCCCCGATTCCTTGTACTCATTAGACTCACTGCGACGGTGATAACGGACATTGTGGATACCCGCCCAGAAGTCGCTAGCCGCTGCAGAGTCTTTGCACAAATGCGGGAACCAGCGAAATTGACTTCTGACGATGGGTGGTCTCTCCTTGCAGCTAAATGATGTGGAAATATTTCACTTCCGGTCCTTACTCATAACCCCCTGAGAGGGTCCGAAGAGTTTCGACGCGGGCAGGCCGAACTGTCTAGTGGAAACAACCAGCAGTACCAAGTCTGATTTCGACATTATCTCGCCGTTATGATTATGAATTGGGTAGCGACGGCCGATACCAGAAGGAGGTAATCATATTTTAGGGGGTGCCGGGACACAGATCGATACCTCAAAGGCCTAACGCAACGCCTCTTTCATCGGGCAGTTTCATTCTCGTAACTGCAGTGAGCAAAGTCACCCCTAAGCAGAAAATTGAGAACTCTAACTGTACCTGTTAACCACACGCTGGGAGGCACCGCTGTGTCCTCCTAACTGTATCTAGTACTCAGCTGGGCCTTCAAATGGTGAGATAATGTCGAAATCCGGGGACTCAACTGGTTCGAGATGTTGGTTTGGACTCACCCACGGCAAACGAGAAGTAACTTGGACTGGGATGACTAAAAGGACTAACTCTTTCATAAGGGCTACGTTCGGCATTTTATGCGGATTGGCATGTACAACTATATGGATCATCAGTATGTGAGTGGGCATCAGGTCGTGTCGCGGGACTTGCCGCGTTGAGGGGGTCCAGTTGCTCACAGGTTGACTGAATATTACGGTTAACTAGGCTCCAAAAGTGAGGTGGTCAGGCCCCGGATTTAGACATAATTTCACCGTTTGACATTATACTTCCCTGACTAAGTCTCTCTATTATTGATAACTGTGTGCGTAATAATTCACAGGCATGCCCCGTTACCTCCGCTGAGTGTTGAGGGGTTGTCCGTCCTACACGATGAGTACTTGCCGTTTAGATGGCTCGTAGTGTAATATTCCTTTGCGTAGTTGACCAGCAGCGGGCGGCGTAATGGTAGTTTGAGAACATAAACGTATTGTAACGAGACAACCAGGCACCTCTTCCATAGATGGTTTAGCTATCAACAGTAGGTGTGGCACGTGGGTCGTGCCATGCCTAGTAGAGGGTTATCACGCGGGAGGCCTCTATCCGACTCAAAAACCGCCTCAGTATCGAAATCCCGTTTGGTGAGGTCCGGTTTAAGCATGTTAGCCATGTAGTGAGTTGACCTCGACTCATCAGGTGCGCAAATAAAGAATGCCCGTTAAACGGCCCTAATACTTCAGATAGGACTGCCATCACACACTCTTTGGCCGAGGCGATCCGCGTTGCGCTTGCGCCTGTCAACCAACTTTTTTTTGGGTTGGGCAAGCAACTGCAGCAACCACGCCGTTCTCCCGATGGGTCAAGGGATAAAAACCAGAGTCACTCGTATGGGATCTTTCGCTTGTCGCAGTCCGGAGGGCTTCAGCCCCGGTAATCACGTAAGTCGGCTAGCTCTGCTTCGGACCACGTATTAGACTACGGCCTGATGTACAGTGCAAAGAGGTTTTGAGCCGGGTGCCCCAACAGGGTAACACGGGTGATAGGAATTGACTGGACCACAGAGGCGCCGAAGGTCGCCAACATCGGGAAGCGGCAACTATGGTCGAAGACCCTGGAAGTAATCTCGACTGAATGGGTTGCATCTGTGCTCCCGCAGTGCGAGGCACCCGCGCAGGTATAGTGTCTTACGGGGTCTATCGCTCTTGCTGAATCGCCTTACCCGAACGGAGAGATTATGTCTAAGTCTTGACCTCACGAATAAAAGAGTAGCTTACCCTACGAGTTAGTGGGCAGAGTCGCGGGTCCACGAATCGATCGTACCCGCTAGACCGGCTACGGTATACACGACATAGAGCTGGAGCCGTCTTGAATCACAATTAACTGTAATGGGGAGATCATGTCTAAATCTTAGACCAGCAATGGGGAGATAATGTCAAAGTCAGTTTGCTAGCTTGCGACATCGACAGCTGCAGATCTTCTACCTAATATAAGACCGGTAATCCGTCCTAGACCATAACGCTGCCTGGATGTCAAGACCGCCAAGATAGCGTTTGTGCGTGACCACAGGGCGCTTTGGGTCCTACCTGGGTTGTCGACCAAAATCTAGACTCCCGTATCAGACTAAAGGCTCGATTCCCTAATAATTGATGTAATAGGTACAGGTTGTGGCCAGGCAGTAGCACCCCGTTAGGTACACGCGTACCCCTCAACTCTATTGTTATCGTGGTTACCTCAGCACGAAGGGGTGGCTACGTCTTATGATTCCGACTACGCAATGACGGCCAGACTGACAGTAGCAGTGGTATGTGCCCGTCCTCTTATGTGCTAGGTCAGATGGTGGGGCGCCTGGATGCACACCTAATTCACGACTGCACATCGCGACTGACCCTCCAGTTCAGACTACTCGACGGACCGGAGGTGAAACTAGGTTCCGATAGACTAGCAAATCCCTGACAAGAGAAGCCCCTCAGACACGAAACGCCAGGATTGAATTCAATTCGAATTGGTGATCCTGCTTCCCTTTATCGGCTCGTGACGCCGGCGAATGACTGATCACCCGTTGGTCGCACGTCCGTTGTTTGTCTTGCTATGCCTTATACACCCACCGGGAGTGGGGTTTCACTGACCACGACTATATTATTGCGTTTGTCTAAGAAGCAAAATAGGACAGGGGCTCAATGTACTAGTAACTGTACCGACGCCGTTGCTTATGGCGGGGGTTAGCGGTTTTAAATCGACTATTAACGCGCACTGGACATAAAGGTTCAACTCCTGGTTAAAGGACAGCAGCACTTATTAAGCGGGGACGCCAGTGCCGACTTTTAGACATAATTTCCCCATTTGTTCTAGGACTTACTCATTATTTCGCCATTAACCATCCACATGTGCGCAAGTCGGCAACTGCATTGACAGGATTTCAGAACATGACCCGACCCCTTTACCCTACATTATGCTGGCAAAGTATATGTCATATGACAGGAACACGTTGAAAAGCATTCTGTGACACTACTGGTTACGAGAGTCTATCCGACTTCATACGGTGTCCTAAGCGTTCTAAATTTCGCACTGGCGAGGAGCTGTCCGCGGTCCCACCCGGCGTCCACACTTTCAGCTACAATTTGCATACTCCAGAATTCACAGCCATGAACCCGTCATCAACGGCGAGATAATGTCAAAGAGCGTCTAACATGACGCACCTTTCTCCGCGAGAGCTGTCAGATGCCAATAGCCCGTGTAAAACACAATGTTAGCTCTGAGATTGGTCTCTGTATTCATCTATCTTAAGGGCTAGCTCCCGAGACCCTTGTGTTAAAGCTGCGATAGGAAGACTTTTCACGAATACGCAACCACAATGACTCTCCCTTGCGTCTGGAAGAGTTGACCCGACACTTGCGCAGCGCTTTGAACTCTTGAAATTACTTAACCGGCGTAACTTGTGCAGGGTCGTATATTCCTTGAACTATCGGAGCGCAAGGTATCGGGAGTAGGTGTACCGATTTGCGCTCTTCGACATAATTTCACCGTTTATTTGCCTTCAACAAAAATTGCCTCCGCACTGCTAATTTTCAGAGTTTCGGGGTCACCTGACCGTACTATGTGTCAGACAGAATCAATGCCATGCATGTTAGTAACAGCATTGCTCTGGAAACTGCGAACTAGGTAGGAAAGGGCCCCAGTCATCTTGAGAGAAACTAGGGAGCGGGGGCATAGGGGCCTCATCGAGATACTCCATTCCGGATGGTGTTTCGTTCGGTTAAAAGACGTCATACGTCCACTCGCTTCGTAGTGGCCATGAGTATTCCTCAGGGAACCTGAAGTCCATTAGGGTGATAAGCTTACCGAAACTTACAACGACTGCATCGATTGGCAACGTAAAAGGGCCCTGACGACCCCGTACCATTACCAGGCTTATTATAAGGTCGGTACACAGCACAGTGATTGAATGATCGACTGTGTTCAGCTCAACACCCCGGACAGACATCTTATCTTAGTTCCTACTCTGCTTATCGGCTACCCCCATTATAGTAAGCCCCAGCTTCGTAACAACAGAACGGGCTTCTGGGACCACGCTAGAAGACAAGAAATGGGTCGAACGAAACCCGAGAAAACCAAGCGCGGTCCAATACGGCACTTATTACAGCGGGTTCAGCTGAGATGGCGACGGAAATTTGATAGACAAAGGCCCTCCGCTGATGATGAAGGGAAGACAATTATTCTACCCCTCGATCTTCGTATGGTGCCTTACTCGCGTCTCTTTCGTGGAGGCTGATGGGTGACTTGCAGTGTGTCTTCATTTCAATCGACCACGGATTGAAGATGACGGACTCTTCGCGTCGCATGAAAGTAGGCCTTACTATTATCTCACCGCGCCCTACTACCGGGTGGCGCCACGCAAACTAATAAGGGCCTAGCACAAAGGTTTTGCTTGTCAGCGAATAAATGGAAAACTCCAGCAAATCCCGGTCGACAGGTTGCACCATACGAATCTCGAGAATGTCTGGCGCGCATGGAAAGTGGAAAACCTAACGGCGAAATCATGTCCAAGAGCGAAATAGTGATCAAGTACTTCAACGAATACCTGTTGGGTTATCCGGCGTTGTGCCTCGAAACATCGTGCCCGGCCGCGGGGGCGACTTGGACATAATCTCGCCATTTGACCTCACATAGCCACCGCGTTTATAAAACATGTTGCCAGTGCAAGGATACCTTTGTCCACAATACTAGTGCTCAGTGTCGATTCCAACGTCAACTAATACAAATTCCTGGACTTAACGTTCTTAATCAAATTTACAGTTATATCACAGCCATAGCTAGTCGTGTAGTTCTCTCGCTGGAAGTTGCGTCGTGTCCCCGTCATAACAATTAGCGCTAGAAATTGTGGACTAGTTAGAACACCTCAGGCGAGAACAGCATTAACTTACGCGACAACGTAAGCACCCACCAGGGATTTTCGTCAATTCGGTCGGATTGGCCAGCCGACGAAGCCGCGAATGCGCCAAACCATGGGCGGTAGGTTCCCTTGAACATCCAGGTAACGCGAATTGAACTCTTATTGAAATGAGCACTGTGCTGACCTTGAGCAGTCAGGGGTCATCACGGATTAGGCGCACCACATGGGACTGTTCCTCATCCGCAGGCTCCGAATGACGTGAATGTAAGCCCGATTACTCAATTTGCAAACCACCCCCGTGCTTAATCGTGGTAATCTAGCCTCTCAATGGTGAAATAATGAGTAAAAGCCGCGAACGAGTATAGCCCCAGAGCCATCGGTATGAAACTTACATGACGAGGGATCTCCTTCCTATTAGGGTTGACGGCTCGGGCAGCAATGCTCCCAACCGACTCACCAATAATGATAGGATCGGGCACTTTATGGAACGTAAATGATAAAGCAGGTCAGAATGGAGAAATAATGAGCAAGTCGTCGCACCAACCTGCTTGTCCACGTATGCAGGAAGACTCCGTAGCGAAGGTCACGTCCAACCGAACCCCCACCTTCTGCTTATGGTACGCGGTTAATAACAGCACTATTGGGTCAGTGCAAGGTTTTTTTCGCATGCAATACAGGTCATGTACTTAGAGGCTAGCAACTGCTGAATCTGCTCCGGATTTAGACATGATTTCTCCGTTCTTGATTAACTTGTACGGCTCTTACTCATAATCTCCCCATTCCGCACCTAGTAAACGTCACTGCGTGTAAGTATGGATCCGAGGTTATGAACGACCCCACCAGGAATCAAGTCGAATTCTCTTCAGGCATATGCCTCAAGAATTTCCATAAACTTGGCAAATTGATAGCCGCTAAGAAGATTCTAACGTTGCTACAAACACCGCCCAGAGCGAACGACAGCCTCGACCGTATCAAGCTTCAACGGTGAAAGTAACTAAATCAAATGAAAGGGTTCGAGGATTTGCAAACTCACAATCAAGCCACGGATAGCGTTAGTAGCGGTGATTAACTTCTCGCTAGTAGCGAAACGGGAATTTCTGCCCACCCTTGCTCACGAATGTTGGACATTTCATCTGGTGTCCGGAGGAGATAAGTACTGTCTGGGACAGGGGCAAGCTCCTCTACACCTGATGGTTTACTTGGTTAAGCTCGGTCAGCAAGATTTGCGAAACCCAACCGTTTCATACGAACGCATAGAGGAAAAACCAATTAAACTCTTGTCGCACCGCAATCGCTATCTCCTCCCTGGGAGATTTCGATGGTTGGTACCGTGGCAGCGTGACCTCTCTATTACGATAGCTGAAGGCTGGATTGACAGGAGTTACATAATACACGACTGACCGCAGAATGCGAGCGGTTGCGGCACATCTCGAACTAAGAGCTGCTTAGGGTATACACGCTCCGACACGACGACTATAGATATCAGAGAGGGACTCGGCCAAGAGGTGTGGATGTCTGAGCCAACACAATTTTCTAATGGCCAAACAAACTCTGACCTGCAACTCAGAATCCCATCTGCCAGTCCACTCCTCGACTCCTTCGGCCTACAGCTATCACCCTGACACGTCACTCTCTTACTACCTACTTTCTCATGGGCTTGACTCTCATTTAATGCGATACATGATTCCTTAGTTAAGATGTGTTGCTCGCCACCTACTGCTAGATCCGGTGGCACCGACGCATTACCGAGCGGAAGCTAGGATTTGTGCCAGATGTTGCGTCCCCTTACGCGATTATCGAACCATAGATCATCGCACGGCCAGAACCTGCCCGATTATATCATGTAGGACAAGCTCATTGGTGGAGAGGTAGAAAAAGAATCTCAGTGAGTTCGCGGTGTCAGTAAAACGACGGGCCTAATAGCAAGGATCAAGACGATCGGCTATGAGGTGCCTAGGATGCAGGCCTCGAATACTCGGATGTATAGGAACAGACCGGGTCGTCTTGACGGTTATGCTGAGGTCGTCGATGTGGTTAAAGTTTATGAGTCCATGTGCGCGCTCAGGATTCGAAGTGATCAGTGCTATCCGCAAGGAGTGTGCGGGTAAGGGTTCAGCGACGCAAGTAATGCGTGACGTACTTCGTATACGTACCATAACACAACTTTCAACAAATGAACTTTATGTCTGCCTCGTGCTCGGAATGGTTGACAGGGCCATTCAAAATAGAGTTCAACTCGCGGTCGCTTTCGGGCAACAGATTTACAAGTTAATAGTAAATCGTCTGGTGCCGGAGGTCCTGGATGCAGTAATAAGTAGAGACAACCTCTCCTGCGATTTCGTGGCACGCCTCGAGTACATGTTCGCCTTCAGCCGGTGTCACGGCTCGACACAACCTACTATTTGGCTTGACTTAAGAAGCCTATGGAGCAGTGTAACGCAAGTGAAGCCTTTCGATTTGCTTGATAGGATCCTAGTGATGTAGCTAAAGCAGTAGCCGGGTACGTTATTTCACTGCTTACTTAGTTGGCCCAAGTCTGGAGTGGGATGTGCAGTGGCAGACCTTGACATAGTGGCGGGGCGTCGGAGTTGAAACAGTATGTCGCGTTATCTTTCCATCCTGGGATGATTTACCAGATATGGACGCGGCTGTACACCCAAACAAGCTTTAGTGGGTCGGCCTTGCCTAGAA"
tyrocidine = "VKLFPWFNQY"

peptide_string = "FFWKCLPHYVSQDQH"
peptide_string = tyrocidine


def main():
    # pprint(peptide_weights(peptide_string))

    print(count_peptides_with_weight(1024))


if __name__ == '__main__':
    main()