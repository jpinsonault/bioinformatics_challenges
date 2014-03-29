from collections import defaultdict
from pprint import pprint
import re

genome = ""
search_string = "CTTGATCAT"


def main():
    matches = re.finditer(r"(?=({}))".format(search_string), genome)

    positions = [str(match.start()) for match in matches]

    print(" ".join(positions))

if __name__ == '__main__':
    main()