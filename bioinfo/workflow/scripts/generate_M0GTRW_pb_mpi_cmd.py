import argparse
import os
import re
from typing import Dict

from bintools.run.run import generate_pb_mpi_cmd


def generate_mapping(local, docker, input) -> str:
    s: str = docker + re.sub(local, "", input)
    return s


def wrapper_generate_M0GTR_pb_mpi_cmd(phylip, tree, output, local, docker) -> str:

    os.makedirs(
        os.path.dirname(output) + "/",
        exist_ok=True,
    )

    kwargs: Dict[str, str] = {
        "-np": str(2),
        "-d": generate_mapping(local=local, docker=docker, input=phylip),
        "-T": generate_mapping(local=local, docker=docker, input=tree),
        "-mutsel": "",
        "-catfix": "uniform",
        "-freeomega": "",
        "-s": "",
        "-x": "10 200",
        "-uni": "",
        "-chainname": generate_mapping(local=local, docker=docker, input=output[:-3]),
    }

    logger: str = " ".join(
        [
            "2>",
            output[:-3] + ".log",
        ]
    )

    cmd: str = generate_pb_mpi_cmd(
        method="pb_mpi", mapping=local + ":" + docker, logger=logger, **kwargs
    )

    with open(output, "w") as fh:
        fh.write(cmd)

    return cmd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="argument", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--phylip",
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
    wrapper_generate_M0GTR_pb_mpi_cmd(
        phylip=args.phylip,
        tree=args.tree,
        output=args.output,
        local=args.local,
        docker=args.docker,
    )
