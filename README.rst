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
    
    sudo pip3 install numpy scipy matplotlib ConfigParser pandas pymultinest

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

The directory structure information
- \es.
	- bin/\DTcomment{Executable file}.
		- easyscan.py.
	- src/\DTcomment{Internal functions}.
		- easyscan\_logging.conf.
		- program.py.
		- initialize.py.
		- readin\_config.py.
		- scanner.py.
		- scan\_controller.py.
		- constraint.py.
		- statfun.py.
		- ploter.py.
		- auxfun.py.
	- utils/\DTcomment{Auxiliary functions}.
		- TestFunction.py.
		- TestFunction\_input.dat.
		- OnePointBatch.in.
		- MSSM\_mW.
		- ALP\_simulation.
	- templates/\DTcomment{Example configuration files}.
		- example\_random.ini.
		- example\_grid.ini.
		- example\_mcmc.ini.
		- example\_mcmc\_bound.ini.
		- example\_multinest.ini.
		- example\_onepoint.ini.
		- example\_onepointbatch.ini.
		- example\_plot.ini.
		- scan\_MSSM\_for\_mW.ini.
		- bound.txt.
	- README.rst \DTcomment{Readme}.
	- LICENSE \DTcomment{Apache license}.
	- EasyScan\_HEP.pdf \DTcomment{Document}.

.. raw:: latex
\dirtree{%
		.1 \es.
		.2 bin/\DTcomment{Executable file}.
		.3 easyscan.py.
		.2 src/\DTcomment{Internal functions}.
		.3 easyscan\_logging.conf.
		.3 program.py.
		.3 initialize.py.
		.3 readin\_config.py.
        .3 scanner.py.
        .3 scan\_controller.py.
        .3 constraint.py.
	    .3 statfun.py.
	   	.3 ploter.py.
	    .3 auxfun.py.
		.2 utils/\DTcomment{Auxiliary functions}.
		.3 TestFunction.py.
		.3 TestFunction\_input.dat.
        .3 OnePointBatch.in.
		.3 MSSM\_mW.
		.3 ALP\_simulation.
		.2 templates/\DTcomment{Example configuration files}.
	    .3 example\_random.ini.
	    .3 example\_grid.ini.
	    .3 example\_mcmc.ini.
	    .3 example\_mcmc\_bound.ini.
	    .3 example\_multinest.ini.
	    .3 example\_onepoint.ini.
	    .3 example\_onepointbatch.ini.
	    .3 example\_plot.ini.
	    .3 scan\_MSSM\_for\_mW.ini.
	    .3 bound.txt.
		.2 README.rst \DTcomment{Readme}.
		.2 LICENSE \DTcomment{Apache license}.
		.2 EasyScan\_HEP.pdf \DTcomment{Document}.
	}
