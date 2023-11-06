=======
EasyScan_HEP
=======

:EasyScan_HEP: A tool for connecting programs to scan the parameter space of physics models
:Author: Yang Zhang， Liangliang Shang
:Version: 1.0
:GitHub: https://github.com/phyzhangyang/EasyScan_HEP
:Website: https://easyscanhep.hepforge.org
:Documentation: https://arxiv.org/pdf/2304.03636.pdf


Installation instructions
-------------------------

EasyScan_HEP is a Python3 code with dependencies on *numpy*, *scipy* and *ConfigParser* libraries. The optional plot functions and MultiNest sampler further require *matplotlib*, *pandas* and *pymultinest* libraries. The dependencies can be installed via pip:
:: 
    sudo apt install python3-pip python3-tk 
    
    pip3 install numpy scipy matplotlib ConfigParser pandas pymultinest

The "easyscan.py" in folder "bin" is the main program, which is executed with configuration file through the command line,
::
    ./bin/easyscan.py templates/example_random.ini

Here *example_random.ini* is an example configuration file provided in EasyScan_HEP. It performs a scan on a simplified model,
:：
    f(x,y) = sin^2 x + cos^2 y,
    
using random sampler, where *x* and *y* are input parameters in range *[0,\pi]* and *[-\pi,\pi]*, respectively, and *f* is output parameter. 

Three other example configuration files in *templates* folder (*example_grid.ini*, *example_mcmc.ini* and *example_multinest.ini*) exhibit briefly usages of other samplers in EasyScan_HEP.

Configuration file *templates/scan_MSSM_for_mW.ini* is an simply physical examples. Relevant programs need to be installed beforehand, using
::
    cd utils/MSSM_mW
    make
    
and then it can be executed with 
::
    ./bin/easyscan.py templates/scan_MSSM_for_mW.ini

Package content:

	- bin
		- easyscan.py
	- src
		- easyscan\_logging.conf
		- program.py
		- initialize.py
		- readin\_config.py
		- scanner.py
		- scan\_controller.py
		- constraint.py
		- statfun.py
		- ploter.py
		- auxfun.py
	- utils
		- TestFunction.py
		- TestFunction\_input.dat
		- OnePointBatch.in
		- MSSM\_mW
	- tools
		- TestFunction.py
		- TestFunction\_input.dat
	- templates
		- example\_random.ini
		- example\_random\_parallel.ini
		- example\_grid.ini
		- example\_mcmc.ini
		- example\_mcmc\_bound.ini
		- example\_multinest.ini
		- example\_onepoint.ini
		- example\_onepointbatch.ini
		- example\_plot.ini
		- scan\_MSSM\_for\_mW.ini
		- bound.txt
	- README.rst 
	- LICENSE 
	
