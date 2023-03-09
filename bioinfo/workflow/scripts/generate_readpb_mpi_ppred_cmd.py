import argparse
import os
import re
from typing import Dict

from bintools.run.run import generate_pb_mpi_cmd


def generate_mapping(local: str, docker: str, input: str) -> str:
    s: str = docker + re.sub(local, "", input)
    return s


def wrapper(
    input: str,
    output: str,
    local: str,
    docker: str,
    image: str,
    sampling: str = "100 10 200",
) -> str:

    os.makedirs(
        os.path.dirname(output) + "/",
        exist_ok=True,
    )

    logger: str = " ".join(
        [
            "2>",
            output[:-3] + ".log",
        ]
    )
    cmd: str = ""
    kwargs: Dict[str, str] = {}
    kwargs = {
        "-np": str(2),
        "-x": sampling,
        "-ppred": "",
        "-chainname": generate_mapping(local=local, docker=docker, input=input[:-6]),
    }
    cmd = generate_pb_mpi_cmd(
        method="readpb_mpi",
        mapping=local + ":" + docker,
        logger=logger,
        image=image,
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
        "--local",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--docker",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--image",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--sampling",
        type=str,
        required=True,
    )

    args = parser.parse_args()
    wrapper(
        input=args.input,
        output=args.output,
        local=args.local,
        docker=args.docker,
        image=args.image,
        sampling=args.sampling,
    )
