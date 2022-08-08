#!/usr/bin/env python
####################################################################
#    Initialization 
####################################################################

# External modules.
import os
import sys
import optparse
import logging
logging.getLogger("matplotlib").setLevel(logging.WARNING)
import logging.config
    
print('\033[1;36;2m')
print('''
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
____ ____ ____ _   _ ____ ____ ____ _  _ _  _ ____ ___
|___ |__| [__   \_/  [__  |    |__| |\ | |__| |___ |__]
|___ |  | ___]   |   ___] |___ |  | | \| |  | |___ |

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
''')
print('\033[0m')

## Set usage options
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

flag_resume = False

## Check whether the configure file exists
try:
    open(sys.argv[1],'r')
except IndexError:
    logger.error('No configfile for the script easyscan')
    print('Usage: ./easyscan configfile')
    sys.exit(1)
except IOError:
    logger.error('Configfile not exist')
    sys.exit(1)

       
