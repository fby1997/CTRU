from math import log
import operator as op
from math import factorial as fac
from math import sqrt, log, ceil
from proba_util import *
from CTRU_failure import *
from NTRU_security import NTRU_summarize_attacks, NTRUParameterSet
from MLWE_security import MLWE_summarize_attacks, MLWEParameterSet


class CTRU_ParameterSet:
    def __init__(self, n, q1, q2, p, probability_distribution1, probability_distribution2=None):
        if probability_distribution2 is None:
            probability_distribution2 = probability_distribution1
        self.n = n
        self.q1 = q1
        self.q2 = q2
        self.q = q1  
        self.p = p
        self.probability_distribution1 = probability_distribution1
        self.probability_distribution2 = probability_distribution2
        self.sigma1 = sqrt (var_of_law(probability_distribution1()))
        self.sigma2 = sqrt (var_of_law(probability_distribution2()))
        self.threshold = 1.0*q1/p

def CTRU_to_NTRU(ps):
    return NTRUParameterSet(ps.n, ps.q1, ps.sigma1)

def CTRU_to_MLWE(ps):
    return MLWEParameterSet(ps.n, ps.q1, ps.sigma2)

def summarize(ps):
    #NTRU_summarize_attacks(CTRU_to_NTRU(ps))
    MLWE_summarize_attacks(CTRU_to_MLWE(ps))
    ErrorRate(ps)
    

if __name__ == "__main__":

    # ParameterSet(n, q1, q2, p, probability_distribution1, probability_distribution2)

    ps_ctru512_cbd2 = CTRU_ParameterSet(653, 4621, 2**10, 2, build_cbd3_law)
    ps_ctru512_cbd3 = CTRU_ParameterSet(653, 4621, 2**11, 2, build_cbd3_law)
    ps_ctru512_cbd4 = CTRU_ParameterSet(653, 4621, 2**12, 2, build_cbd3_law)       

    ps_ctru768_1024_cbd2 = CTRU_ParameterSet(761, 4591, 2**10, 2, build_cbd2_law)
    ps_ctru768_3457_cbd3 = CTRU_ParameterSet(761, 4591, 2**11, 2, build_cbd2_law)      
    ps_ctru768_2048_cbd3 = CTRU_ParameterSet(761, 4591, 2**12, 2, build_cbd2_law)      

    # ps_ctru512_cbd2  
    ps_temp = ps_ctru512_cbd2
    print("**************************************")
    print("ps_ctru512_cbd2:")
    print("CTRU: n = %d, q = %d, q2 = %d, p = %d\n"%(ps_temp.n, ps_temp.q1, ps_temp.q2, ps_temp.p))  
    Bandwidth(ps_temp)                               
    summarize(ps_temp)
    print()
    
    # ps_ctru512_cbd3  
    ps_temp = ps_ctru512_cbd3
    print("**************************************")
    print("ps_ctru512_cbd3:")
    print("CTRU: n = %d, q = %d, q2 = %d, p = %d\n"%(ps_temp.n, ps_temp.q1, ps_temp.q2, ps_temp.p))  
    Bandwidth(ps_temp)                               
    summarize(ps_temp)
    print()
    
        # ps_ctru512_cbd4
    ps_temp = ps_ctru512_cbd4
    print("**************************************")
    print("ps_ctru512_cbd3:")
    print("CTRU: n = %d, q = %d, q2 = %d, p = %d\n"%(ps_temp.n, ps_temp.q1, ps_temp.q2, ps_temp.p))  
    Bandwidth(ps_temp)                               
    summarize(ps_temp)
    print()
     
    # ps_ctru768_1024_cbd2  
    ps_temp = ps_ctru768_1024_cbd2
    print("**************************************")
    print("ps_ctru768_1024_cbd2:")
    print("CTRU: n = %d, q = %d, q2 = %d, p = %d\n"%(ps_temp.n, ps_temp.q1, ps_temp.q2, ps_temp.p))  
    Bandwidth(ps_temp)                               
    summarize(ps_temp)
    print()
    
    # ps_ctru768_3457_cbd3  
    ps_temp = ps_ctru768_3457_cbd3
    print("**************************************")
    print("ps_ctru768_3457_cbd3:")
    print("CTRU: n = %d, q = %d, q2 = %d, p = %d\n"%(ps_temp.n, ps_temp.q1, ps_temp.q2, ps_temp.p))  
    Bandwidth(ps_temp)                               
    summarize(ps_temp)
    print()
    
    # ps_ctru768_2048_cbd3  
    ps_temp = ps_ctru768_2048_cbd3
    print("**************************************")
    print("ps_ctru768_2048_cbd3:")
    print("CTRU: n = %d, q = %d, q2 = %d, p = %d\n"%(ps_temp.n, ps_temp.q1, ps_temp.q2, ps_temp.p))  
    Bandwidth(ps_temp)                               
    summarize(ps_temp)
    print()