#########################################################################
#    Statistical funtions                                               #
#########################################################################
from math import log10


def prior(x,par):
    if par[1] == 'flat':
        return x * (par[3] - par[2]) + par[2]
    elif par[1] == 'log':
        return 10.0** ( x * (log10(par[3]) - log10(par[2])) + log10(par[2]) )
        
        
        
        

