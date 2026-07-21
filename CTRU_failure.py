"""Failure-probability estimators for CTRU, CTRU-Light, and CTRU-Prime."""

from math import ceil, exp, log, pi, sqrt
from scipy.stats import chi2
from proba_util import build_centered_binomial_law, var_of_law


def build_uniform_law(k):
    return {value: 1.0 / (2 * k + 1) for value in range(-k, k + 1)}


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


def geometric_mean(values):
    return exp(sum(log(value) for value in values) / len(values))


def log2_sum(log_values):
    maximum = max(log_values)
    return maximum + log(
        sum(2 ** (value - maximum) for value in log_values), 2
    )


def e8_volume_threshold_factor():
    return (384 / pi**4) ** (1 / 64)


def block_tail_log2(variances, threshold, method):
    """Approximate an eight-dimensional weighted Gaussian tail."""
    if method == "geometric":
        effective_variance = geometric_mean(variances)
        return chi2.logsf(threshold**2 / effective_variance, 8) / log(2)

    if method == "satterthwaite":
        variance_sum = sum(variances)
        variance_square_sum = sum(value**2 for value in variances)
        scale = variance_square_sum / variance_sum
        degrees_of_freedom = variance_sum**2 / variance_square_sum
        return chi2.logsf(
            threshold**2 / scale, degrees_of_freedom
        ) / log(2)

    raise ValueError(
        "block_variance_method must be 'geometric' or 'satterthwaite'"
    )


def build_rounding_law_rlwr(ps):
    distribution = {}
    for u in range(ps.q1):
        temp = u - ps.q1 / ps.q2 * round(ps.q2 / ps.q1 * u)
        epsilon = -ps.q2 / ps.q1 * temp
        distribution[epsilon] = distribution.get(epsilon, 0) + 1.0 / ps.q1
    return distribution


def ErrorRate_CTRU_3Cyclo_Ring(
    ps,
    use_e8_volume_factor=False,
    block_variance_method="geometric",
):
    """Error rate for CTRU over Z[x]/(x^n-x^(n/2)+1)."""
    sigma_epsilon = sqrt(var_of_law(build_rounding_law_rlwr(ps)))

    s1_list = []
    for i in range(int(ps.n / 2) - 1):
        count = 3 * ps.n / 2 - i - 1
        s1_list.append(count * ps.sigma1**2 * ps.sigma2**2)
    for _ in range(2):
        s1_list.append(ps.n * ps.sigma1**2 * ps.sigma2**2)
    for _ in range(int(ps.n / 2) + 1, ps.n):
        s1_list.append(3 * ps.n / 2 * ps.sigma1**2 * ps.sigma2**2)

    s2_list = []
    for i in range(int(ps.n / 2) - 1):
        count = 3 * ps.n / 2 - i - 1
        s2_list.append(4 * count * ps.sigma1**2 * sigma_epsilon**2)
    for _ in range(2):
        s2_list.append(4 * ps.n * ps.sigma1**2 * sigma_epsilon**2)
    for _ in range(int(ps.n / 2) + 1, ps.n):
        s2_list.append(
            4 * (3 * ps.n / 2) * ps.sigma1**2 * sigma_epsilon**2
        )

    variance_list = [
        variance_gr
        + (ps.q1 / ps.q2) ** 2 * variance_epsilon_f
        for variance_gr, variance_epsilon_f in zip(s1_list, s2_list)
    ]

    threshold = ps.threshold
    if use_e8_volume_factor:
        threshold *= e8_volume_threshold_factor()

    block_logs = []
    for start in range(0, ps.n, 8):
        block_logs.append(
            block_tail_log2(
                variance_list[start : start + 8],
                threshold,
                block_variance_method,
            )
        )

    result = log2_sum(block_logs)
    print(f"err rate: 2^({result:.2f})")
    return result


# Compatibility with the name used in the supplied formula.
ErrorRate_CTRU_3Cyclo_Ring_2 = ErrorRate_CTRU_3Cyclo_Ring


def build_rounding_law_rlwr_1(ps):
    distribution = {}
    for u in range(ps.q1):
        base_term = ps.q2 / ps.q1 * u
        epsilon = round(base_term) - base_term
        distribution[epsilon] = distribution.get(epsilon, 0) + 1.0 / ps.q1
    return distribution


def build_rounding_law_rlwr_2(ps):
    lower = ceil(-ps.q1 / (2 * ps.q2))
    upper = ceil(ps.q1 / (2 * ps.q2))
    length = upper - lower
    distribution = {}
    for u in range(lower, upper):
        epsilon = -ps.q2 / ps.q1 * u
        distribution[epsilon] = distribution.get(epsilon, 0) + 1.0 / length
    return distribution


def _ctru_light_standard_deviation(ps, include_message_term=False):
    sigma_epsilon_1 = sqrt(var_of_law(build_rounding_law_rlwr_1(ps)))
    sigma_epsilon_2 = sqrt(var_of_law(build_rounding_law_rlwr_2(ps)))

    variance_gr = ps.n * ps.sigma1**2 * ps.sigma2**2
    variance_epsilon_1_f = 2 * ps.n * ps.sigma1**2 * sigma_epsilon_1**2
    variance = (
        variance_gr
        + (ps.q1 / ps.q2) ** 2
        * (variance_epsilon_1_f + sigma_epsilon_1**2)
        + 4 * ps.n * ps.sigma1**2 * sigma_epsilon_2**2
        + sigma_epsilon_2**2
    )
    if include_message_term:
        variance += 1.0 / 8
    return sqrt(variance)


def _ctru_light_error_rate(ps, threshold, include_message_term):
    standard_deviation = _ctru_light_standard_deviation(
        ps, include_message_term=include_message_term
    )
    block_log = chi2.logsf(
        (threshold / standard_deviation) ** 2, 8
    ) / log(2)
    result = log(ps.n / 16, 2) + block_log
    print(f"err rate: 2^({result:.2f})")
    return result


def ErrorRate_CTRU_Light_v1(ps):
    return _ctru_light_error_rate(
        ps, threshold=ps.threshold, include_message_term=False
    )


def ErrorRate_CTRU_Light_v2(ps):
    return _ctru_light_error_rate(
        ps, threshold=ps.threshold3, include_message_term=True
    )


def ErrorRate_CTRU_Light_v3(ps):
    return _ctru_light_error_rate(
        ps, threshold=ps.threshold2, include_message_term=True
    )


# Default CTRU-Light estimator requested for this package.
ErrorRate_CTRU_Light = ErrorRate_CTRU_Light_v1


def _prime_term_counts(n):
    counts = [n]
    counts.extend(2 * n - k for k in range(1, n - 1))
    counts.append(n + 1)
    return counts


def ErrorRate_CTRU_Prime_Field(
    ps,
    use_e8_volume_factor=True,
    block_variance_method="geometric",
):
    """CTRU-Prime estimator"""
    sigma_epsilon2 = var_of_law(build_rounding_law_rlwr(ps))
    variance_list = []
    for count in _prime_term_counts(ps.n):
        variance_gr = count * ps.sigma1**2 * ps.sigma2**2
        variance_epsilon_f = 8 * count * ps.sigma1**2 * sigma_epsilon2
        variance_list.append(
            variance_gr + (ps.q1 / ps.q2) ** 2 * variance_epsilon_f
        )

    threshold = ps.threshold
    if use_e8_volume_factor:
        threshold *= e8_volume_threshold_factor()

    block_logs = []
    for start in range(0, ps.n, 8):
        block = variance_list[start : start + 8]
        block_logs.append(
            block_tail_log2(block, threshold, block_variance_method)
        )

    result = log2_sum(block_logs) + log(8, 2)
    print(f"err rate: 2^({result:.2f})")
    return result


def Bandwidth(ps):
    public_key = ceil(ps.n * ceil(log(ps.q1, 2)) / 8)
    ciphertext = ceil(ps.n * ceil(log(ps.q2, 2)) / 8)
    print(
        f"|pk| = {public_key}, |ct| = {ciphertext}, "
        f"bandwidth = {public_key + ciphertext}"
    )
