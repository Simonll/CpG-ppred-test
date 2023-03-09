import argparse
import glob
import os
from typing import List

import numpy as np
import pandas as pd
from suffstat_parser import suffstat_parser


def wrapper(input_dir: str, pattern: str, output: str) -> bool:

    os.makedirs(
        os.path.dirname(output) + "/",
        exist_ok=True,
    )
    list_of_df: List[pd.DataFrame] = []
    list_of_zipfiles: List[str] = glob.glob(input_dir + pattern)
    for zipf in list_of_zipfiles:
        geneID: str = zipf.split("/")[-1].split("-")[0]
        modelID: str = zipf.split("/")[-1].split("-")[1]
        repID: str = zipf.split("/")[-1].split("-")[3]
        df: pd.DataFrame = suffstat_parser(
            input_zip=zipf,
            metadata={"geneID": geneID, "modelID": modelID, "repID": repID},
        )
        try:
            df["CG>TG|CA"] = df[["CG>TG", "CG>CA"]].sum(axis=1)
            df["GC>GT|AC"] = df[["GC>GT", "GC>AC"]].sum(axis=1)
            sum_: pd.DataFrame = (
                df[
                    [
                        "geneID",
                        "modelID",
                        "repID",
                        "type",
                        "mcmcID",
                        "CG>TG|CA",
                        "CG",
                        "GC>GT|AC",
                        "GC",
                    ]
                ]
                .groupby(["geneID", "modelID", "repID", "type", "mcmcID"])
                .agg([np.sum])
                .droplevel(level=1, axis=1)
            )
            list_of_df += [sum_]
            print(".", end="")
        except Exception as e:
            print("something wrong when computing stats %s" % str(e))
            raise RuntimeError
        try:
            pd.concat(list_of_df, axis=0, ignore_index=False).to_csv(output, sep="\t")
        except Exception as e:
            print("something wrong concatenating %s" % str(e))
    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="argument", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--input_dir",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--pattern",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
    )

    args = parser.parse_args()
    wrapper(
        input_dir=args.input_dir,
        pattern=args.pattern,
        output=args.output,
    )
