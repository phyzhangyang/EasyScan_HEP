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
from pathlib import Path
    
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
parser.add_option("-r","--resume", action="store_true", default=False,
        dest='resume', help='launch resume mode')
        
(options, args) = parser.parse_args()
if len(args) == 0:
    args = ''

def logging_config_path():
    candidates = []
    env_root = os.environ.get("EASYSCAN_ROOT", "")
    if env_root:
        candidates.append(Path(env_root) / "src" / "easyscan_logging.conf")
    candidates.append(Path(__file__).resolve().with_name("easyscan_logging.conf"))
    candidates.append(Path(__file__).resolve().parents[1] / "src" / "easyscan_logging.conf")
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return None


## Configure logging info
if options.debug and options.logging == 'INFO':
    options.logging = 'DEBUG'
log_config = logging_config_path()
if log_config:
    logging.config.fileConfig(str(log_config), defaults={'logfilename': 'EASYSCAN.LOG'})
else:
    logging.basicConfig(
        level=getattr(logging, options.logging),
        format='%(levelname)s: %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler('EASYSCAN.LOG', mode='w'),
        ],
    )
logging.root.setLevel(vars(logging)[options.logging])
logging.getLogger('easyscan').setLevel(vars(logging)[options.logging])
                              
logger=logging.getLogger('easyscan.main')

resume = options.resume

## Check whether the configure file exists
try:
    open(sys.argv[1],'r')
except IndexError:
    logger.error('No configfile for the script easyscan')
    print('Usage: easyscan configfile')
    sys.exit(1)
except IOError:
    logger.error('Configfile not exist')
    sys.exit(1)

       
