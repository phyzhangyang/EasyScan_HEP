#########################################################################
#    Statistical funtions                                               #
#########################################################################
from math import log10

def prior(x,par):
    prior_type = par[1].lower()
    if prior_type == 'flat':
        max_ = float(par[3])
        min_ = float(par[2])
        return x * (max_ - min_) + min_ 
    elif prior_type == 'log':
        max_ = log10(float(par[3]))
        min_ = log10(float(par[2]))
        return 10.0** ( x * (max_ - min_) + min_ )
    elif prior_type == 'fixed':
        return float(par[2])
    else:
        sf.ErrorStop( 'Not ready. Only "flat" and "log" prior can be used.' )
