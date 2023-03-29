#def sig_noerr(bk_xs_pb, bk_eff, alp_xs_pb, alp_eff):
#    import numpy as np
#    import sympy as sy
#
#    Lum = sy.Symbol('Lum', positive=True, real=True)
#
#    sn = alp_xs_pb*1000*alp_eff
#    bk = bk_xs_pb*1000*bk_eff
#    eq = sy.sqrt(2*Lum(sn-bk*sy.log(1+sn/bk ) ) ) 
#
#    L = sy.nsolve(eq,Lum,1,verify=False)
#    
#    return L
#

def sig_noerr(bk_xs_pb, bk_eff, alp_xs_pb, alp_eff):
    import math

    return 4.0/2.0/(sn - bk * math.log(1 + sn/bk))

L = sig_noerr(replace_1, replace_2, replace_3, replace_4)

f = open("significance_output.dat", "w")
f.write("lum_fbinv\n")
f.write(str(L))
f.close()
