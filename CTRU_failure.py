from scipy.stats import chi2
from math import sqrt, log, ceil, erf
from math import factorial as fac
from proba_util import *

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
    R = ps.probability_distribution2
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
    R = ps.probability_distribution2
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
def ErrorRate(ps):
    sigma_epsilon_1 = sqrt( var_of_law( build_rounding_law_rlwr_non_power_of_two_2(ps,int(ps.n*5/4) ) ))    
    s1_1 = ps.n/2*(ps.sigma1**2 * ps.sigma2**2 + ps.sigma2**2 * (ps.sigma1**2+ps.sigma1**2) ) 
    temp_1 = (ps.q1/ps.q2*sigma_epsilon_1)**2
    s2_1 = ps.n/2*(ps.p**2)*(temp_1 * ps.sigma1**2 + temp_1 * (ps.sigma1**2+ps.sigma1**2) ) 
    s3_1 = (ps.q1/ps.q2*sigma_epsilon_1)**2  
    s_1 = sqrt(s1_1 + s2_1 + s3_1)
    pr_1 = chi2.logsf( (ps.threshold/s_1)**2, 8 ) / log(2) + log(ps.n/16, 2)

    sigma_epsilon_2 = sqrt( var_of_law( build_rounding_law_rlwr_non_power_of_two_2(ps,int(ps.n*3/2) ) ))    
    s1_2 = ps.n/2*(ps.sigma1**2 * ps.sigma2**2 + ps.sigma2**2 * (ps.sigma1**2+ps.sigma1**2) ) 
    temp_2 = (ps.q1/ps.q2*sigma_epsilon_2)**2
    s2_2 = ps.n/2*(ps.p**2)*(temp_2 * ps.sigma1**2 + temp_2 * (ps.sigma1**2+ps.sigma1**2) ) 
    s3_2 = (ps.q1/ps.q2*sigma_epsilon_2)**2  
    s_2 = sqrt(s1_2 + s2_2 + s3_2)
    pr_2 = chi2.logsf( (ps.threshold/s_2)**2, 8 ) / log(2) + log(ps.n/16, 2)
    
    pr = log(2**(pr_1)+2**(pr_2),2)
    
    print("err:")
    print("    = 2^%.2f"% pr)


def Bandwidth(ps):
    pk = ceil(ps.n*ceil(log(ps.q1,2))/8)
    ct = ceil(ps.n*ceil(log(ps.q2,2))/8)
    print('|pk| = %d, |ct| = %d, bandwidth = %d\n'%(pk, ct, (pk+ct)))
    
