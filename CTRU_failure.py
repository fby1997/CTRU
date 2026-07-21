from scipy.stats import chi2
from math import sqrt, log, ceil, erf
from math import factorial as fac
from proba_util import *
import numpy as np

def build_rounding_law_rlwr(ps):
    D = {}
    for u in range(0, ps.q1):    
        temp = u - 1.*ps.q1/ps.q2*round(ps.q2/ps.q1*u)
        epsilon = -ps.q2/ps.q1*temp
        D[epsilon] = D.get(epsilon,0)+1./ps.q1 
    return D  

def build_rounding_law_message(ps):
    D = {}
    for k1 in range(0, ps.p):    
        temp = ps.q1/ps.p*k1
        epsilon = round( temp ) - temp
        D[epsilon] = D.get(epsilon,0)+1./ps.p    
    return D   
        

def build_rounding_law_ciphertext(ps):
    D = {}
    for sigma1 in range(0, ps.q1):
        for k1 in range(0, ps.p):    
            temp = 1.*ps.q2/ps.q1*( sigma1 + round(ps.q1/ps.p*k1))
            if temp >= ps.q2:
                temp -= ps.q2
            if temp <= -ps.q2:
                temp += ps.q2                
            epsilon = round( temp ) - temp
            if epsilon >= ps.q2/2:
                epsilon -= ps.q2
            if epsilon <= -ps.q2/2:
                epsilon += ps.q2              
            D[epsilon] = D.get(epsilon,0)+1./ps.q1 * 1./ps.p    
    return D  
def build_rounding_law_rlwr_non_power_of_two_2(ps,k):
    H = {}
    for h in range(0, ps.q1):
        H[h] = 1./ps.q1
    R = ps.probability_distribution2()
    HR_each = law_product(H, R)
    HR = iter_law_convolution_modulo_q(HR_each,k,ps.q1)
    #HR = iter_law_convolution(HR_each, int(ps.n/2))
    #HR = iter_law_convolution_modulo_q(HR_each, int(ps.n/2), ps.q1)
    C = {}
    for i in HR:
        c = i % ps.q1
        C[c] = C.get(c, 0) + HR[i]
    D = {}
    for c in C:    
        temp = c - 1.*ps.q1/ps.q2*round(ps.q2/ps.q1*c)
        epsilon = -ps.q2/ps.q1*temp
        D[epsilon] = D.get(epsilon,0) + C[c] 
    return D  
def build_rounding_law_rlwr_non_power_of_two_3(ps):
    H = {}
    for h in range(0, ps.q1):
        H[h] = 1./ps.q1
    R = ps.probability_distribution2()
    HR_each = law_product_over_non_power_of_2(H, R);
    #HR = iter_law_convolution(HR_each, int(ps.n/2))
    HR = iter_law_convolution_modulo_q(HR_each, int(ps.n*2.75/6), ps.q1)
    C = {}
    for i in HR:
        c = i % ps.q1
        C[c] = C.get(c, 0) + HR[i]
    D = {}
    for c in C:    
        temp = c - 1.*ps.q1/ps.q2*round(ps.q2/ps.q1*c)
        epsilon = -ps.q2/ps.q1*temp
        D[epsilon] = D.get(epsilon,0) + C[c] 
    return D 
# CTRU-Prime
def ErrorRate_Prime(ps):
    sigma_epsilon = sqrt( var_of_law( build_rounding_law_rlwr_non_power_of_two_3(ps) ) )    
    s1 = int(ps.n*2.75/6)*(ps.sigma1**2 * ps.sigma2**2 + ps.sigma2**2 * (ps.sigma1**2+ps.sigma1**2) ) 
    temp = (ps.q1/ps.q2*sigma_epsilon)**2
    s2 = int(ps.n*2.75/6)*(ps.p**2)*(temp * ps.sigma1**2 + temp * (ps.sigma1**2+ps.sigma1**2) ) 
    s3 = (ps.q1/ps.q2*sigma_epsilon)**2  
    s = sqrt(s1 + s2 + s3)
    pr = chi2.logsf( (ps.threshold/s)**2, 8 ) / log(2) + log(ps.n/8, 2)
    print("err:")
    print("    = 2^%.2f"% pr)
# CTRU
from math import exp, log, sqrt
from scipy.stats import chi2
from proba_util import *


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
        self.n = n
        self.q1 = q1
        self.q2 = q2
        self.q = q1
        self.p = p
        self.probability_distribution1 = probability_distribution1
        self.probability_distribution2 = probability_distribution2
        self.sigma1 = sqrt(var_of_law(probability_distribution1))
        self.sigma2 = sqrt(var_of_law(probability_distribution2))
        self.threshold = 1.0 * q1 / p


def build_rounding_law_rlwr(ps):
    D = {}
    for u in range(ps.q1):
        temp = u - 1.0 * ps.q1 / ps.q2 * round(ps.q2 / ps.q1 * u)
        epsilon = -ps.q2 / ps.q1 * temp
        D[epsilon] = D.get(epsilon, 0) + 1.0 / ps.q1
    return D


def geometric_mean(data):
    return exp(sum(log(x) for x in data) / len(data))


def log2_sum(log_values):
    maximum = max(log_values)
    return maximum + log(sum(2 ** (x - maximum) for x in log_values), 2)


def ErrorRate(ps):
    sigma_epsilon = sqrt(var_of_law(build_rounding_law_rlwr(ps)))

    # Coefficient term counts for the common 3-cyclotomic Ring-Half model:
    #   7n/4 for the first n/2-1 coefficients,
    #   n    for the two central coefficients,
    #   3n/2 for the remaining n/2-1 coefficients.
    term_count_list = (
        [7 * ps.n / 4] * (int(ps.n / 2) - 1)
        + [ps.n] * 2
        + [3 * ps.n / 2] * (int(ps.n / 2) - 1)
    )

    # For a coefficient containing t product terms:
    #   Var(g*r)       = t*sigma_1^2*sigma_2^2,
    #   Var(epsilon*f) = p^2*t*sigma_1^2*sigma_epsilon^2.
    s1_list = [
        count * ps.sigma1**2 * ps.sigma2**2
        for count in term_count_list
    ]
    s2_list = [
        ps.p**2 * count * ps.sigma1**2 * sigma_epsilon**2
        for count in term_count_list
    ]

    block_log_probabilities = []
    for i in range(0, ps.n, 8):
        block_variance = (
            geometric_mean(s1_list[i : i + 8])
            + (ps.q1 / ps.q2) ** 2
            * geometric_mean(s2_list[i : i + 8])
        )
        block_log_probabilities.append(
            chi2.logsf(ps.threshold**2 / block_variance, 8) / log(2)
        )
    rho = log(ps.n / 16, 2) - log(768 / 16, 2)
    result = (
        log2_sum(block_log_probabilities)
        + log(16, 2)
        - 38 * rho
        - 275 * rho**2
    )
    print(f"err rate: 2^({result:.2f})")
    return result

# Helper functions for probability distributions
def build_cbd1_law(): return build_centered_binomial_law(1)
def build_cbd2_law(): return build_centered_binomial_law(2)
def build_cbd3_law(): return build_centered_binomial_law(3)
def build_cbd4_law(): return build_centered_binomial_law(4)
def build_cbd5_law(): return build_centered_binomial_law(5)

if __name__ == "__main__":
    p = 2
    parameter_sets = [
        ("CTRU-576", 576, build_cbd5_law()),
        ("CTRU-768", 768, build_cbd3_law()),
        ("CTRU-1024", 1024, build_cbd2_law()),
    ]

    print("======== CTRU error-rate model version_202505 ========")
    for name, n, distribution in parameter_sets:
        print(f"\n--- {name} ---")
        ps = CTRU_ParameterSet(n, 3457, 2**10, p, distribution)
        ErrorRate(ps)

def Bandwidth(ps):
    pk = ceil(ps.n*ceil(log(ps.q1,2))/8)
    ct = ceil(ps.n*ceil(log(ps.q2,2))/8)
    print('|pk| = %d, |ct| = %d, bandwidth = %d\n'%(pk, ct, (pk+ct)))

def geometric_mean(data):
    return exp(sum(log(x) for x in data) / len(data))

def satterthwaite_effective_degrees_of_freedom(variances, degrees_of_freedom):
    """
    计算Satterthwaite近似下的有效自由度。
    
    参数:
    variances (list of float): 各个分布的方差列表。
    degrees_of_freedom (list of int): 各个分布的自由度列表。
    
    返回:
    float: 近似卡方分布的有效自由度。
    """
    numerator = (np.sum(variances))**2
    denominator = np.sum([(var**2) / df for var, df in zip(variances, degrees_of_freedom)])
    nu_eff = numerator / denominator
    
    return nu_eff
    
