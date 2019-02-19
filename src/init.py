#!/usr/bin/env python
####################################################################
#    Initialization and some usefull functions.
####################################################################

# External modules.
import os
import sys
import optparse
import logging
import logging.config
import time

print '\033[1;36;2m'
print '''
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
____ ____ ____ _   _ ____ ____ ____ _  _ _  _ ____ ___
|___ |__| [__   \_/  [__  |    |__| |\ | |__| |___ |__]
|___ |  | ___]   |   ___] |___ |  | | \| |  | |___ |

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
'''
print '\033[0m'

## Add debug option
usage = "usage: %prog [options] [FILE] "
parser = optparse.OptionParser(usage=usage)
parser.add_option("-l", "--logging", default='INFO',
        help="logging level (DEBUG|INFO|WARNING|ERROR|CRITICAL) [%default]")
parser.add_option("-d","--debug", action="store_true", default=False,
        dest='debug', help='force to launch debug mode')
(options, args) = parser.parse_args()
if len(args) == 0:
    args = ''

## Configure logging info
if options.debug and options.logging == 'INFO':
    options.logging = 'DEBUG'
logging.config.fileConfig(os.path.join(os.path.split(os.path.split(os.path.realpath(__file__))[0])[0], "src/easyscan_logging.conf"),
        defaults={'logfilename': 'EASYSCAN.LOG'})
logging.root.setLevel(vars(logging)[options.logging])
logging.getLogger('easyscan').setLevel(vars(logging)[options.logging])
                              
logger=logging.getLogger('easyscan.main')

##check whether the configure file exists
try:
    open(sys.argv[1],'r')
except IndexError:
    logger.error('No configfile for the script easyscan')
    print 'Usage: ./easyscan configfile'
    sys.exit(1)
except IOError:
    logger.error('Configfile not exist')
    sys.exit(1)


## define 'negetive infinity'
log_zero = -1e+100
## define 'not a number'
NaN = float('NaN')
## Current path
CurrentPath = os.getcwd()

def ColorText(i,text,j=1):
    return '\033[%i;3%i;2m %s\033[0m' %(j,i,text)

def GotoWeb():
    print ColorText(1,'# Goto ') + ColorText(1,'http://easyscanhep.hepforge.org',4) + ColorText(1,' for detail.')

def WarningWait(warinfo):
    logger.warning(ColorText(1,warinfo))
    print ColorText(1,'# Waiting 3 seconds for WARNING.')
    time.sleep(3)

def WarningNoWait(warinfo):
    logger.warning(ColorText(1,warinfo))

def ErrorStop(errinfo):
    logger.error( ColorText(1,errinfo) )
    GotoWeb()
    print ColorText(1,'# Exiting with ERROR.')
    sys.exit(1)

def Info(debinfo):
    logger.info( ColorText(2,debinfo) )

def Debug(debinfo,debvalue=''):
    if debvalue=='':
        logger.debug( ColorText(5, str(debinfo) ) )
    else:
        logger.debug( ColorText(5, str(debinfo)+': '+str(debvalue) ) )


def autotype(s):
    if type(s) is not str:
        return s
    if s.isdigit():
        return int(s)
    try:
        return float(s)
    except ValueError:
        return s

def string2list(s):
    s = [ autotype(ss.strip()) for ss in s.split('\n') ]
    return s

## used for parsing string of input variable and output variable in configure file to list of list of items.
def string2nestlist(s):
    s = map( lambda x: x.split(','), s.split('\n') )
    s = [[autotype(x.strip()) for x in ss] for ss in s]
    return s

## "File" parameter readd 20180416 liang
def WriteResultInf(InPar,OutPar,Chi2,Path, ScanMethod,File):
    if ScanMethod == 'PLOT': return
    #if ScanMethod == 'READ': os.rename(os.path.join(Path,'ScanInf.txt'),os.path.join(Path,'ScanInf_old.txt'))
    file_inf = open(os.path.join(Path,'ScanInf.txt'),'w')
    #unmark 20180416 liang
    file_inf.write('\t'.join([Path, File])+'\n')
    i   = 0
    if ScanMethod == 'MULTINEST':
        file_inf.write('probability\t0\n')
        file_inf.write('-2*loglikehood\t1\n')
        i = 2
    for name in InPar:
        file_inf.write('\t'.join([name,str(i)])+'\n')
        i += 1
    for name in OutPar :
        file_inf.write('\t'.join([name,str(i)])+'\n')
        i += 1
    ## new 20180428 liang
    for name in Chi2:
        file_inf.write('\t'.join([name,str(i)])+'\n')
        i += 1
    if ScanMethod == 'MCMC':
        file_inf.write('\t'.join(['mult',str(i+0)])+'\n')
    file_inf.close()

def parseMath(par):
    ## Thanks to authors at the web page http://lybniz2.sourceforge.net/safeeval.html

    from math import *
    ## make a list of safe functions
    safe_list = ['math','acos', 'asin', 'atan', 'atan2', 'ceil', 'cos', 'cosh',
                 'degrees', 'e', 'exp', 'fabs', 'floor', 'fmod', 'frexp', 'hypot',
                 'ldexp', 'log', 'log10', 'modf', 'pi', 'pow', 'radians', 'sin',
                 'sinh', 'sqrt', 'tan', 'tanh']

    ## use the list to filter the local namespace
    safe_dict = dict([ (k, locals().get(k, None)) for k in safe_list ])
    ## add any needed builtins back in.
    safe_dict['abs'] = abs
    safe_dict['float'] = float 

    safe_dict.update(par)
    safe_dict.update({"__builtins__": None})

    for key,value in par.items():
        expr = ','.join(key.split(':'))

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

#new 20180419 liang
def checkItemInList(List):
    for item in List:
        counter = List.count(item)
        if counter>1:
            ErrorStop('Figure name / output variable name "%s" duplicating %i times! Please correct in [plot] / [programX] in your input file!!'%(item, counter))
        
#new 20180420 liang
def checkFileInList(List):
    newList=[]
    files=[]
    for item in List:
        try:
            float(item)
            newList.append(item)
        except ValueError:
            if item.find('/') >= 0:
                files.append(item)
                newList.append(item.split('/')[-1])
            else:
                newList.append(item)

    return newList, files
 
## new 20180428 liang
def sortDic(Dic):
    from collections import OrderedDict
    return OrderedDict(sorted(Dic.items(), key = lambda t: t[0]))

       
