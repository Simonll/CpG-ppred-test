import argparse
import os

from bintools.align.align import ali
from bintools.align.align import read_fasta
from bintools.align.align import write_phylip


def fasta2phylip(fasta: str, phylip: str) -> bool:
    try:
        os.makedirs(os.path.dirname(phylip), exist_ok=True)
    except Exception as e:
        print(
            "somehting wrong when making dir %s, %s" % (os.path.dirname(phylip), str(e))
        )
        return False
    try:
        with open(fasta, "r") as fh:
            ali_: ali = read_fasta(fh=fh)
            write_phylip(
                filename=phylip,
                align=ali_.get_biopython_align(),
            )
    except Exception as e:
        print("somehting wrong when converting %s to %s, %s" % (fasta, phylip, str(e)))
        return False
    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="argument", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--fasta",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--phylip",
        type=str,
        required=True,
    )

    args = parser.parse_args()
    fasta2phylip(fasta=args.fasta, phylip=args.phylip)
