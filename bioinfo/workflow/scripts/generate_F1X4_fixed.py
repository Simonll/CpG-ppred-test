import argparse
import os
from pathlib import Path
from typing import Dict

from bintools.align.align import read_phylip
from bintools.phylobayes.mcmc_parser import posterior_M0_GTR
from bintools.phylobayes.priors import prior_M7M8_fix


def generate_MGF1x4_fixed(phylip, mcmc, burnin, omega, output) -> bool:

    os.makedirs(
        os.path.dirname(output) + "/",
        exist_ok=True,
    )
    Nsite: int = -1
    with open(phylip) as stream:
        Nsite = int(read_phylip(fh=stream).get_n_site() / 3)
        pMGF1x4: posterior_M0_GTR = posterior_M0_GTR.parse_mcmc(
            mcmc_path=Path(mcmc), burnin=int(burnin)
        )
        params: Dict[str, str] = pMGF1x4.sample()
        params.update({"omega_site": prior_M7M8_fix(N=Nsite, value=float(omega))})
        if os.path.exists(output):
            append_write = "a"  # append if already exists
        else:
            append_write = "w"  # make a new file if not

        with open(
            output,
            append_write,
        ) as fl:
            pMGF1x4.write_values(dict_of_params=params, file_handler=fl)
    return True


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
        "--mcmc",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--burnin",
        type=int,
        required=True,
    )

    parser.add_argument(
        "--omega",
        type=float,
        required=True,
    )

    parser.add_argument(
        "--output",
        type=str,
        required=True,
    )

    args = parser.parse_args()

    generate_MGF1x4_fixed(
        phylip=args.phylip,
        mcmc=args.mcmc,
        burnin=args.burnin,
        omega=args.omega,
        output=args.output,
    )
