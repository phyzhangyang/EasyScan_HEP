#!/usr/bin/env python
####################################################################
#    Some usefull functions
####################################################################

# External modules.
import time, os, sys
from collections import OrderedDict
# Internal modules
from initialize import logger
from initialize import resume

log_zero = -1e+100
NaN = float('NaN')
CurrentPath = os.getcwd()
# Name of scanner
_onepoint = "ONEPOINT"
_onepointbatch = "ONEPOINTBATCH"
_random = "RANDOM"
_grid = "GRID"
_mcmc = "MCMC"
_multinest = "MULTINEST"
_plot = "PLOT"
_postprocess = "POSTPROCESS"
_read = "READ"

_all = [_onepoint, _onepointbatch, _random, _grid, _mcmc, _multinest, _postprocess, _plot, _read]
_no_random = [_onepoint, _onepointbatch, _grid, _postprocess, _plot, _read]
_no_like   = [_onepoint, _onepointbatch, _random, _grid, _postprocess, _plot, _read]
_post = [_postprocess, _plot]
# Define name of result data file
ResultFile = 'ScanResult.txt'
ResultFile_MCMC = 'All_ScanResult.txt'
ResultFile_MultiNest = 'MultiNestData/.txt'
ResultFile_post = 'Previous_ScanResult.txt'


# Define screen print functions
def ColorText(i,text,j=1):
    return '\033[%i;3%i;2m %s\033[0m' %(j,i,text)
def GotoWeb():
    print(ColorText(1,'# Goto ') + ColorText(1,'https://github.com/phyzhangyang/EasyScan_HEP',4) + ColorText(1,' for detail.'))
def WarningWait(warinfo):
    logger.warning(ColorText(1,warinfo))
    print(ColorText(1,'# Waiting 3 seconds for WARNING.'))
    time.sleep(3)
def WarningNoWait(warinfo):
    logger.warning(ColorText(1,warinfo))
def ErrorStop(errinfo):
    logger.error( ColorText(1,errinfo) )
    print(ColorText(1,'# Exiting with ERROR.'))
    sys.exit(1)
def Info(debinfo):
    logger.info( ColorText(2,debinfo) )
def Debug(debinfo,debvalue=''):
    if debvalue=='':
        logger.debug( ColorText(5, str(debinfo) ) )
    else:
        logger.debug( ColorText(5, str(debinfo)+': '+str(debvalue) ) )

# Transform str into int or float if possible
def autotype(s):
    if type(s) is not str:
        return s
    else:
        if s.isdigit():
            return int(s)
        else:
            try:
                return float(s)
            except ValueError:
                return s

# Transform str to int or float forcely 
def forcetype(s):
    if type(s) is not str:
        ErrorStop("%s is not str"%s)
    if s.isdigit():
        return int(s)
    try:
        return float(s)
    except:
        ErrorStop("%s is not int or float. Please check settings for Output variable in [program]."%s)

# Check if input is a number
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


# Transform str into list
def string2list(s):
    s = [ autotype(ss.strip()) for ss in s.split('\n') ]
    return s

# Parse string of input variable and output variable in configure file to list of items.
def string2nestlist(s):
    s = [x.split(',') for x in s.split('\n')]
    s = [[autotype(x.strip()) for x in ss] for ss in s]
    return s

# Write information of result into file
def WriteResultInf(InPar, FixedPar, OutPar, Constraint, Path, ScanMethod):
    if ScanMethod == _plot or resume: 
      return
    inf = ''
    if ScanMethod == _multinest:
      inf += 'probability,-2lnlike,'
    inf += ','.join(list(InPar.keys())+list(FixedPar.keys())+list(OutPar.keys())+list(Constraint.keys()))

    if ScanMethod == _mcmc:
      inf += ",dwell"
    outfile = open(os.path.join(Path, ResultFile),'w')
    outfile.write(inf+'\n')
    outfile.close()
    if ScanMethod == _mcmc:
      outfile = open(os.path.join(Path, ResultFile_MCMC),'w')
      outfile.write(inf+'\n')
      outfile.close()
    
# Evaluate a math string
# http://lybniz2.sourceforge.net/safeeval.html
import math 
# Make a list of safe functions, add any needed builtins back in.
safe_dict = {'abs':abs, 'float':float, 'int':int,
             'acos':math.acos, 'asin':math.asin, 'atan':math.atan, 'atan2':math.atan2, 
             'ceil':math.ceil, 'cos':math.cos, 'cosh':math.cosh, 'degrees':math.degrees, 
             'e':math.e, 'exp':math.exp, 'fabs':math.fabs, 'floor':math.floor, 'fmod':math.fmod, 'frexp':math.frexp, 
             'hypot':math.hypot, 'ldexp':math.ldexp, 'log':math.log, 'log10':math.log10, 
             'modf':math.modf, 'pi':math.pi, 'pow':math.pow, 'radians':math.radians, 
             'sin':math.sin, 'sinh':math.sinh, 'sqrt':math.sqrt, 
             'tan':math.tan, 'tanh':math.tanh,
             '__builtins__': None}
# names that can not be used for variable
forbidden_names = [*safe_dict] + _all + ['dwell', 'probability', '-2lnlike'] 
def parseMath(par):
    safe_dict.update(par)
    for key,value in par.items():
        expr = ','.join(key.split(';'))
        try:
            cal = eval(expr, safe_dict)
        except SyntaxError:
            for keyAlt,valueAlt in par.items():
                if type(valueAlt) == str:
                    expr=expr.replace(keyAlt, valueAlt)
                else:
                    if isnan(valueAlt): 
                        expr=expr.replace(keyAlt, "float('NaN')")
                    else:
                        expr=expr.replace(keyAlt, str(valueAlt))
            cal = eval(expr, safe_dict)

        #print key, expr, cal; raw_input("math") 
        par[key] = cal

# Sort parameters so that we can output them in right order
def sortDic(Dic):
    return OrderedDict(sorted(list(Dic.items()), key = lambda t: t[0]))
       
# Divide samples for multi-processes
def divide_jobs(num_jobs, num_works, i_work):
    return int(num_jobs/num_works) + 1 if i_work < num_jobs%num_works else int(num_jobs/num_works)
