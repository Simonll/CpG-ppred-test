import argparse
import os
import re
from typing import Dict

from bintools.run.run import generate_iqtree_cmd


def generate_mapping(local, docker, input) -> str:
    s: str = docker + re.sub(local, "", input)
    return s


def wrapper(simu, length, values, tree, output, local, docker) -> str:

    os.makedirs(
        os.path.dirname(output) + "/",
        exist_ok=True,
    )

    kwargs: Dict[str, str] = {}

    with open(values, "r") as stream:
        gtr4g = stream.read().strip()

    kwargs = {
        "--alisim": generate_mapping(local, docker, simu),
        "--seqtype": "DNA",
        "-af": "fasta",
        "-t": generate_mapping(local, docker, tree),
        "-m": gtr4g,
        "--length": str(length),
    }

    logger: str = " ".join(["2>", output[:-3] + ".log"])

    cmd: str = generate_iqtree_cmd(
        method="iqtree2",
        mapping=local + ":" + docker,
        logger=logger,
        image="ubuntu20.04/iqtree:latest",
        **kwargs
    )

    with open(output, "w") as fh:
        fh.write(cmd)

    return cmd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="argument", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--simu",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--length",
        type=int,
        required=True,
    )
    parser.add_argument(
        "--values",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--tree",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--local",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--docker",
        type=str,
        required=True,
    )

    args = parser.parse_args()
    wrapper(
        simu=args.simu,
        length=args.length,
        values=args.values,
        tree=args.tree,
        output=args.output,
        local=args.local,
        docker=args.docker,
    )
