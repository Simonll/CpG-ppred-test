import argparse
import re
from typing import List


def wrapper(input: str, output: str):
    rx = re.compile(r"""(?<!\\)%.+|(\\(?:no)?citep?\{((?!\*)[^{}]+)\})""")

    with open(input, "r") as fh:
        lines = fh.readlines()

        list_of_authors: List[str] = []
        for l in lines:
            list_of_authors += [m.group(2) for m in rx.finditer(l) if m.group(2)]

        print(list_of_authors)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="argument", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
