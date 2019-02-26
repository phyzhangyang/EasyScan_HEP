#########################################################################
#    Statistical funtions                                               #
#########################################################################
from math import log10
import init as sf


def prior(x,par):
    if par[1].lower() == 'flat':
        max_ = float(par[3])
        min_ = float(par[2])
        return x * (max_ - min_) + min_ 
    elif par[1].lower() == 'log':
        max_ = log10(float(par[3]))
        min_ = log10(float(par[2]))
        return 10.0** ( x * (max_ - min_) + min_ )
    else:
        sf.ErrorStop( 'Not ready. Only "flat" and "log" prior can be used.' )
