import argparse
import os
from pathlib import Path
from typing import Dict
from typing import Union

from bintools.phylobayes.gtrg_parser import posterior_GTRG


def wrapper(mcmc, burnin, output_values, output_tree) -> bool:

    os.makedirs(
        os.path.dirname(output_values) + "/",
        exist_ok=True,
    )
    os.makedirs(
        os.path.dirname(output_tree) + "/",
        exist_ok=True,
    )

    pGTR: posterior_GTRG = posterior_GTRG.parse_mcmc(
        mcmc_path=Path(mcmc), burnin=int(burnin)
    )
    params: Dict[str, Union[float, str]] = pGTR.sample()

    with open(
        output_values,
        "w",
    ) as fl:

        rho_sum: float = (
            (float(params["rho_AC"]))
            + (float(params["rho_AG"]))
            + (float(params["rho_AT"]))
            + (float(params["rho_CG"]))
            + (float(params["rho_CT"]))
            + (float(params["rho_TG"]))
        )
        phi_sum: float = (
            float(params["phi_A"])
            + float(params["phi_C"])
            + float(params["phi_G"])
            + float(params["phi_T"])
        )
        fl.write(
            "GTR{"
            + str(float(params["rho_AC"]) / rho_sum)
            + "/"
            + str(float(params["rho_AG"]) / rho_sum)
            + "/"
            + str(float(params["rho_AT"]) / rho_sum)
            + "/"
            + str(float(params["rho_CG"]) / rho_sum)
            + "/"
            + str(float(params["rho_CT"]) / rho_sum)
            + "/"
            + str(float(params["rho_TG"]) / rho_sum)
            + "}+F{"
            + str(float(params["phi_A"]) / phi_sum)
            + "/"
            + str(float(params["phi_C"]) / phi_sum)
            + "/"
            + str(float(params["phi_G"]) / phi_sum)
            + "/"
            + str(float(params["phi_T"]) / phi_sum)
            + "}"
            + "+G4{"
            + str(float(params["alpha"]))
            + "}\n",
        )
    with open(
        output_tree,
        "w",
    ) as fl:
        fl.write(str(params["tree"]) + "\t")

    return True


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="argument", formatter_class=argparse.ArgumentDefaultsHelpFormatter
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
        "--output_values",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--output_tree",
        type=str,
        required=True,
    )

    args = parser.parse_args()

    wrapper(
        mcmc=args.mcmc,
        burnin=args.burnin,
        output_tree=args.output_tree,
        output_values=args.output_values,
    )
