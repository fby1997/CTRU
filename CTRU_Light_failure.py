"""CTRU-Light decryption-failure estimator."""

from math import log, sqrt

from scipy.stats import chi2

from CTRU_failure import build_rounding_law_rlwr
from proba_util import build_centered_binomial_law, var_of_law


def ErrorRate_CTRU_Light(ps):
    """Return log2(delta) for CTRU-Light (n=512, q=769, B1)."""
    sigma1 = ps.sigma1
    sigma2 = ps.sigma2
    sigma_epsilon = sqrt(var_of_law(build_rounding_law_rlwr(ps)))

    # Power-of-two cyclotomic ring: the reduced expected term count is 2.75n/6.
    term_count = int(ps.n * 2.75 / 6)
    s1 = term_count * (
        sigma1**2 * sigma2**2
        + sigma2**2 * (sigma1**2 + sigma1**2)
    )

    scaled_epsilon_variance = 0.5 * (
        ps.q1 / ps.q2 * sigma_epsilon
    ) ** 2
    s2 = term_count * ps.p**2 * (
        scaled_epsilon_variance * sigma1**2
        + scaled_epsilon_variance * (sigma1**2 + sigma1**2)
    )
    s3 = scaled_epsilon_variance
    standard_deviation = sqrt(s1 + s2 + s3)

    result = chi2.logsf(
        (ps.threshold / standard_deviation) ** 2, 8
    ) / log(2) + 2 * log(ps.n / 8, 2)
    print(f"err rate: 2^({result:.2f})")
    return result


def build_cbd1_law():
    return build_centered_binomial_law(1)


if __name__ == "__main__":
    from CTRU import CTRU_ParameterSet

    p = 2
    ps = CTRU_ParameterSet(512, 769, 2**8, p, build_cbd1_law)
    ErrorRate_CTRU_Light(ps)
