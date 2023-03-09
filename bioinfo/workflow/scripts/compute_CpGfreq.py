import argparse
import ast
import codecs
import os
from itertools import product
from typing import Any
from typing import Dict
from zipfile import ZipFile

import numpy as np
import pandas as pd
from bintools.align.align import ali
from bintools.align.align import read_fasta
from bintools.align.align import read_phylip
from bintools.align.utils import compute_XpY


def wrapper(input: str, output: str, metadata: str) -> bool:

    os.makedirs(
        os.path.dirname(output) + "/",
        exist_ok=True,
    )

    metadata_: Dict[str, Any] = ast.literal_eval(metadata)
    ss: Dict[str, float] = {}
    if input.endswith(".zip"):
        try:
            dict_of_stats = {}
            k = 0
            with ZipFile(input, "r") as zarch:

                for f in zarch.namelist():
                    drawID: str = f.split("A_")[-1].split(".ali")[0]
                    with zarch.open(f, "r") as fh:
                        align: ali = read_phylip(fh=codecs.getreader("utf-8")(fh))
                        # print(align.get_dict_of_seq())
                        dict_of_XpY = {}
                        for key, v in align.get_dict_of_seq().items():
                            # print(key,v)
                            seq: str = "".join([s for p, s in v.items()])

                            for x, y in list(
                                product(["A", "C", "G", "T"], ["A", "C", "G", "T"])
                            ):
                                XpY: str = x + "p" + y
                                if XpY not in dict_of_XpY:
                                    dict_of_XpY[XpY] = [
                                        compute_XpY(seq=seq, x=x, y=y)[XpY]
                                    ]
                                else:
                                    dict_of_XpY[XpY] += [
                                        compute_XpY(seq=seq, x=x, y=y)[XpY]
                                    ]

                        dict_of_stats[k] = {
                            x + "p" + y: np.mean(dict_of_XpY[x + "p" + y])
                            for x, y in list(
                                product(["A", "C", "G", "T"], ["A", "C", "G", "T"])
                            )
                        }
                        dict_of_stats[k].update({"drawID": drawID})

                        for key, v in metadata_.items():
                            dict_of_stats[k].update({key: v})
                        k += 1

        except Exception as e:
            print("something wrong when reading ziped suffstats %s, %s" % (f, str(e)))
            raise RuntimeError
        try:
            df = pd.DataFrame.from_dict(data=dict_of_stats, orient="index")
        except Exception as e:
            print(
                "something wrong when concatenating files from zip %s, %s" % (f, str(e))
            )
            raise RuntimeError

        try:
            df.to_csv(output, sep="\t")
        except Exception as e:
            print("something wrong writing tsv file %s %s" % (output, str(e)))

    else:
        dict_of_stats = {}
        k = 0
        with open(input, "r") as fh_:
            if input.endswith(".ali") or input.endswith(".phylip"):
                align = read_phylip(fh=fh_)
            else:
                align = read_fasta(fh=fh_)

            dict_of_XpY = {}
            for key, v in align.get_dict_of_seq().items():
                # print(key,v)
                seq = "".join([s for p, s in v.items()])

                for x, y in list(product(["A", "C", "G", "T"], ["A", "C", "G", "T"])):
                    XpY = x + "p" + y
                    if XpY not in dict_of_XpY:
                        dict_of_XpY[XpY] = [compute_XpY(seq=seq, x=x, y=y)[XpY]]
                    else:
                        dict_of_XpY[XpY] += [compute_XpY(seq=seq, x=x, y=y)[XpY]]

            dict_of_stats[k] = {
                x + "p" + y: np.mean(dict_of_XpY[x + "p" + y])
                for x, y in list(product(["A", "C", "G", "T"], ["A", "C", "G", "T"]))
            }
            for key, v in metadata_.items():
                dict_of_stats[k].update({key: v})
            k += 1
        try:
            df = pd.DataFrame.from_dict(data=dict_of_stats, orient="index")
        except Exception as e:
            print(
                "something wrong when concatenating files from zip %s, %s" % (f, str(e))
            )
            raise RuntimeError
        try:
            df.to_csv(output, sep="\t")
        except Exception as e:
            print("something wrong writing tsv file %s %s" % (output, str(e)))

    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="argument", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--input",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--metadata",
        type=str,
        required=False,
    )
    args = parser.parse_args()
    wrapper(input=args.input, output=args.output, metadata=args.metadata)
