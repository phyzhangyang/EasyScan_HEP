#########################################################################
#    Statistical funtions                                               #
#########################################################################
from math import log10

def prior(x,par):
    prior_type = par[1].upper()
    if prior_type == 'FLAT':
        max_ = float(par[3])
        min_ = float(par[2])
        return x * (max_ - min_) + min_ 
    elif prior_type == 'LOG':
        max_ = log10(float(par[3]))
        min_ = log10(float(par[2]))
        return 10.0** ( x * (max_ - min_) + min_ )
    elif prior_type == 'FIXED':
        return float(par[2])
    # TODO add log_flat, gauss
    else:
        sf.ErrorStop( 'Not ready. Only "flat" and "log" prior can be used.' )
