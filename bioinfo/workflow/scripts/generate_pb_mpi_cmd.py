import argparse
import os
import re
from typing import Dict

from bintools.run.run import generate_pb_mpi_cmd


def generate_mapping(local, docker, input) -> str:
    s: str = docker + re.sub(local, "", input)
    return s


def wrapper(
    phylip,
    tree,
    model,
    sampling,
    np,
    local_root,
    docker_root,
    cluster_root,
    work_dir,
    sh_docker,
    sh_cluster,
) -> str:

    os.makedirs(
        os.path.dirname(local_root + work_dir) + "/",
        exist_ok=True,
    )

    kwargs_docker: Dict[str, str] = {
        "-d": generate_mapping(local=local_root, docker=docker_root, input=phylip),
        "-T": generate_mapping(local=local_root, docker=docker_root, input=tree),
        model: "",
        sampling: "",
        "-np": np,
        "-chainname": generate_mapping(
            local=local_root, docker=docker_root, input=work_dir + sh_docker[:-3]
        ),
    }

    logger_local: str = " ".join(
        [
            "2>",
            local_root + work_dir + sh_docker[:-3] + ".log",
        ]
    )
    cmd: str = ""

    cmd = generate_pb_mpi_cmd(
        method="pb_mpi",
        mapping=local_root + ":" + docker_root,
        logger=logger_local,
        **kwargs_docker
    )

    with open(local_root + work_dir + sh_docker, "w") as fh:
        fh.write(cmd)

    with open(local_root + work_dir + sh_cluster, "w") as fh:
        fh.write(
            "#qsub -cwd -pe mpi %s %s"
            % (np, cluster_root + work_dir + sh_cluster + "\n")
        )
        fh.write("LD_LIBRARY_PATH=/usr/local/lib:/usr/local/lib64:$LD_LIBRARY_PATH\n")
        fh.write("CXX=`which g++`\n")
        fh.write("CC=`which gcc`\n")
        fh.write("export LD_LIBRARY_PATH\n")
        fh.write("export CC\n")
        fh.write("export CXX\n")
        fh.write(
            "mpirun -np %s /home/sll/pbmpi/data/pb_mpi -d %s -T %s %s %s %s"
            % (
                np,
                cluster_root + phylip,
                cluster_root + tree,
                model,
                sampling,
                cluster_root + work_dir + sh_docker[:-3],
            )
        )

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
        "--model",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--sampling",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--np",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--local_root",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--docker_root",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--cluster_root",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--work_dir",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--sh_docker",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--sh_cluster",
        type=str,
        required=True,
    )

    args = parser.parse_args()
    wrapper(
        phylip=args.phylip,
        tree=args.tree,
        model=args.model,
        sampling=args.sampling,
        np=args.np,
        local_root=args.local_root,
        docker_root=args.docker_root,
        cluster_root=args.cluster_root,
        work_dir=args.work_dir,
        sh_cluster=args.sh_cluster,
        sh_docker=args.sh_docker,
    )
