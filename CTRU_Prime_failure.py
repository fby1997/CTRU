"""CTRU-Prime decryption-failure estimator."""

from math import log, pi, sqrt

from scipy.stats import chi2

from CTRU_failure import build_rounding_law_rlwr, geometric_mean, log2_sum
from proba_util import build_centered_binomial_law, var_of_law


def e8_volume_threshold_factor():
    return (384 / pi**4) ** (1 / 64)


def ErrorRate_CTRU_Prime_Field(ps, use_e8_volume_factor=True):
    """Return log2(delta) for the CTRU-Prime parameter set."""
    sigma1 = ps.sigma1
    sigma2 = ps.sigma2
    sigma_epsilon2 = var_of_law(build_rounding_law_rlwr(ps))

    # Z[x]/(x^n-x-1): n, then 2n-k (1<=k<=n-2), then n+1.
    counts = [ps.n]
    counts.extend(2 * ps.n - k for k in range(1, ps.n - 1))
    counts.append(ps.n + 1)

    variance_list = []
    for count in counts:
        variance_gr = count * sigma1**2 * sigma2**2
        variance_epsilon_f = 8 * count * sigma1**2 * sigma_epsilon2
        variance_list.append(
            variance_gr + (ps.q1 / ps.q2) ** 2 * variance_epsilon_f
        )

    threshold = ps.threshold
    if use_e8_volume_factor:
        threshold *= e8_volume_threshold_factor()

    block_logs = []
    for start in range(0, ps.n, 8):
        effective_variance = geometric_mean(variance_list[start : start + 8])
        block_logs.append(
            chi2.logsf(
                threshold**2 / effective_variance, 8
            ) / log(2)
        )

    result = log2_sum(block_logs) + log(8, 2)
    print(f"err rate: 2^({result:.2f})")
    return result


def build_cbd1_law():
    return build_centered_binomial_law(1)


def build_cbd2_law():
    return build_centered_binomial_law(2)


if __name__ == "__main__":
    from CTRU import CTRU_ParameterSet

    p = 2
    ps = CTRU_ParameterSet(761, 4591, 2**10, p, build_cbd2_law)
    ErrorRate_CTRU_Prime_Field(ps)
