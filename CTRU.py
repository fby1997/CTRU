"""Parameter sets and command-line entry point for CTRU variants."""

from math import sqrt

from CTRU_failure import Bandwidth, ErrorRate, build_cbd2_law, build_cbd3_law, build_cbd5_law
from CTRU_Light_failure import ErrorRate_CTRU_Light, build_cbd1_law
from CTRU_Prime_failure import ErrorRate_CTRU_Prime_Field
from MLWE_security import MLWEParameterSet, MLWE_summarize_attacks
from NTRU_security import NTRUParameterSet, NTRU_summarize_attacks
from proba_util import var_of_law


class CTRU_ParameterSet:
    def __init__(
        self,
        n,
        q1,
        q2,
        p,
        probability_distribution1,
        probability_distribution2=None,
    ):
        if probability_distribution2 is None:
            probability_distribution2 = probability_distribution1

        distribution1 = (
            probability_distribution1()
            if callable(probability_distribution1)
            else probability_distribution1
        )
        distribution2 = (
            probability_distribution2()
            if callable(probability_distribution2)
            else probability_distribution2
        )

        self.n = n
        self.q1 = q1
        self.q2 = q2
        self.q = q1
        self.p = p
        self.probability_distribution1 = probability_distribution1
        self.probability_distribution2 = probability_distribution2
        self.sigma1 = sqrt(var_of_law(distribution1))
        self.sigma2 = sqrt(var_of_law(distribution2))
        self.threshold = q1 / p


def CTRU_to_NTRU(ps):
    return NTRUParameterSet(ps.n, ps.q1, ps.sigma1)


def CTRU_to_MLWE(ps):
    return MLWEParameterSet(ps.n, ps.q1, ps.sigma2)


def summarize_security(ps):
    NTRU_summarize_attacks(CTRU_to_NTRU(ps))
    MLWE_summarize_attacks(CTRU_to_MLWE(ps))


def print_parameters(name, ps, distribution_name):
    print(f"\n--- {name} ---")
    print(
        f"Parameters: n={ps.n}, q={ps.q1}, q2={ps.q2}, p={ps.p}, "
        f"dist=({distribution_name}, {distribution_name})"
    )


def main():
    p = 2
    ctru_parameter_sets = [
        ("CTRU-576", CTRU_ParameterSet(576, 3457, 2**10, p, build_cbd5_law), "B5"),
        ("CTRU-768", CTRU_ParameterSet(768, 3457, 2**10, p, build_cbd3_law), "B3"),
        ("CTRU-1024", CTRU_ParameterSet(1024, 3457, 2**10, p, build_cbd2_law), "B2"),
    ]

    print("======== CTRU Parameter Sets ========")
    for name, ps, distribution_name in ctru_parameter_sets:
        print_parameters(name, ps, distribution_name)
        Bandwidth(ps)
        ErrorRate(ps)

    print("\n======== CTRU-Light Parameter Set ========")
    light = CTRU_ParameterSet(512, 769, 2**8, p, build_cbd1_law)
    print_parameters("CTRU-Light", light, "B1")
    Bandwidth(light)
    ErrorRate_CTRU_Light(light)

    print("\n======== CTRU-Prime Parameter Set ========")
    prime = CTRU_ParameterSet(761, 4591, 2**10, p, build_cbd2_law)
    print_parameters("CTRU-Prime-761", prime, "B2")
    Bandwidth(prime)
    ErrorRate_CTRU_Prime_Field(prime, use_e8_volume_factor=True)


if __name__ == "__main__":
    main()
