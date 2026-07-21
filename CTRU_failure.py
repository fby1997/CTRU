"""CTRU decryption-failure estimator.
"""

from math import exp, log, sqrt
from scipy.stats import chi2
from proba_util import *


def _distribution(value):
    """Consume either a distribution dictionary or its constructor."""
    return value() if callable(value) else value


def build_rounding_law_rlwr(ps):
    distribution = {}
    for u in range(ps.q1):
        temp = u - ps.q1 / ps.q2 * round(ps.q2 / ps.q1 * u)
        epsilon = -ps.q2 / ps.q1 * temp
        distribution[epsilon] = distribution.get(epsilon, 0) + 1.0 / ps.q1
    return distribution


def geometric_mean(values):
    return exp(sum(log(value) for value in values) / len(values))


def log2_sum(values):
    maximum = max(values)
    return maximum + log(sum(2 ** (value - maximum) for value in values), 2)


def _ring_half_term_counts(n):
    return (
        [7 * n / 4] * (int(n / 2) - 1)
        + [n] * 2
        + [3 * n / 2] * (int(n / 2) - 1)
    )


def ErrorRate(ps):
    """Return log2(delta) for CTRU-576, CTRU-768, or CTRU-1024."""
    sigma1 = sqrt(var_of_law(_distribution(ps.probability_distribution1)))
    sigma2 = sqrt(var_of_law(_distribution(ps.probability_distribution2)))
    sigma_epsilon2 = var_of_law(build_rounding_law_rlwr(ps))

    counts = _ring_half_term_counts(ps.n)
    variance_gr = [
        count * sigma1**2 * sigma2**2 for count in counts
    ]
    variance_epsilon_f = [
        ps.p**2 * count * sigma1**2 * sigma_epsilon2 for count in counts
    ]

    block_logs = []
    for start in range(0, ps.n, 8):
        block_variance = (
            geometric_mean(variance_gr[start : start + 8])
            + (ps.q1 / ps.q2) ** 2
            * geometric_mean(variance_epsilon_f[start : start + 8])
        )
        block_logs.append(
            chi2.logsf(ps.threshold**2 / block_variance, 8) / log(2)
        )

    # The quadratic term is the finite-dimension correction to the geometric
    # block-variance approximation.
    rho = log(ps.n / 16, 2) - log(768 / 16, 2)
    result = log2_sum(block_logs) + log(16, 2) - 38 * rho - 275 * rho**2
    print(f"err rate: 2^({result:.2f})")
    return result


def build_cbd1_law():
    return build_centered_binomial_law(1)


def build_cbd2_law():
    return build_centered_binomial_law(2)


def build_cbd3_law():
    return build_centered_binomial_law(3)


def build_cbd4_law():
    return build_centered_binomial_law(4)


def build_cbd5_law():
    return build_centered_binomial_law(5)


def Bandwidth(ps):
    from math import ceil

    public_key = ceil(ps.n * ceil(log(ps.q1, 2)) / 8)
    ciphertext = ceil(ps.n * ceil(log(ps.q2, 2)) / 8)
    print(
        f"|pk| = {public_key}, |ct| = {ciphertext}, "
        f"bandwidth = {public_key + ciphertext}"
    )
