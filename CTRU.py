"""Unified failure-rate test entry point for CTRU variants. version_202607"""

from math import sqrt
from CTRU_failure import (
    Bandwidth,
    ErrorRate_CTRU_3Cyclo_Ring,
    ErrorRate_CTRU_Light_v1,
    ErrorRate_CTRU_Prime_Field,
    build_cbd1_law,
    build_cbd2_law,
    build_cbd3_law,
    build_cbd4_law,
    build_uniform_law,
)
from proba_util import var_of_law
from NTRU_security import NTRU_summarize_attacks, NTRUParameterSet
from MLWE_security import MLWE_summarize_attacks, MLWEParameterSet


class CTRU_ParameterSet:
    def __init__(
        self,
        n,
        q1,
        q2,
        p,
        probability_distribution1,
        probability_distribution2=None,
        probability_distribution3=None,
    ):
        if probability_distribution2 is None:
            probability_distribution2 = probability_distribution1
        if probability_distribution3 is None:
            probability_distribution3 = probability_distribution1

        self.n = n
        self.q1 = q1
        self.q2 = q2
        self.q = q1
        self.p = p
        self.probability_distribution1 = probability_distribution1
        self.probability_distribution2 = probability_distribution2
        self.probability_distribution3 = probability_distribution3
        self.sigma1 = sqrt(var_of_law(probability_distribution1))
        self.sigma2 = sqrt(var_of_law(probability_distribution2))
        self.sigma3 = sqrt(var_of_law(probability_distribution3))
        self.threshold = q1 / p
        self.threshold2 = (q1 - 1) / p
        self.threshold3 = (q1 + 1) / p


def print_parameters(name, ps, distribution_name):
    print(f"\n--- {name} ---")
    print(
        f"Parameters: n={ps.n}, q={ps.q1}, q2={ps.q2}, p={ps.p}, "
        f"dist=({distribution_name}, {distribution_name})"
    )

def CTRU_to_NTRU(ps):
    return NTRUParameterSet(ps.n, ps.q1, ps.sigma1)

def CTRU_to_MLWE(ps):
    return MLWEParameterSet(ps.n, ps.q1, ps.sigma2)

def summarize(ps):
    #NTRU_summarize_attacks(CTRU_to_NTRU(ps))
    MLWE_summarize_attacks(CTRU_to_MLWE(ps))

def run_ctru_parameter(name, ps, distribution_name):
    print_parameters(name, ps, distribution_name)
    Bandwidth(ps)
    summarize(ps)
    for method in ("geometric", "satterthwaite"):
        for use_e8 in (False, True):
            print(
                f"CTRU method={method}, "
                f"e8_volume_factor={use_e8}"
            )
            ErrorRate_CTRU_3Cyclo_Ring(
                ps,
                use_e8_volume_factor=use_e8,
                block_variance_method=method,
            )


def main():
    p = 2

    run_ctru_parameter(
        "CTRU-576",
        CTRU_ParameterSet(576, 3457, 2**10, p, build_cbd4_law()),
        "B4",
    )
    run_ctru_parameter(
        "CTRU-768",
        CTRU_ParameterSet(768, 3457, 2**10, p, build_cbd3_law()),
        "B3",
    )
    run_ctru_parameter(
        "CTRU-1024",
        CTRU_ParameterSet(
            1024, 3457, 2**10, p, build_cbd3_law(), build_cbd2_law()
        ),
        "B3/B2",
    )
    run_ctru_parameter(
        "CTRU-1152",
        CTRU_ParameterSet(
            1152, 3457, 2**10, p, build_cbd3_law(), build_cbd2_law()
        ),
        "B3/B2",
    )
    run_ctru_parameter(
        "CTRU-1536",
        CTRU_ParameterSet(1536, 3457, 2**10, p, build_uniform_law(1)),
        "U1",
    )
    run_ctru_parameter(
        "CTRU-2048",
        CTRU_ParameterSet(2048, 3457, 2**10, p, build_cbd1_law()),
        "B1",
    )

    print("\n======== CTRU-Light Parameter Set ========")
    light = CTRU_ParameterSet(512, 769, 2**8, p, build_cbd1_law())
    print_parameters("CTRU-Light-512", light, "B1")
    Bandwidth(light)
    summarize(light)
    print("CTRU-Light default: v1")
    ErrorRate_CTRU_Light_v1(light)

    print("\n======== CTRU-Prime Parameter Set ========")
    prime = CTRU_ParameterSet(761, 4591, 2**10, p, build_cbd2_law())
    print_parameters("CTRU-Prime-761", prime, "B2")
    Bandwidth(prime)
    summarize(prime)
    print("CTRU-Prime default: geometric + E8 volume factor")
    ErrorRate_CTRU_Prime_Field(
        prime,
        use_e8_volume_factor=True,
        block_variance_method="geometric",
    )


if __name__ == "__main__":
    main()
