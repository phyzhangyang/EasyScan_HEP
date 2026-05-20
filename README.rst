============
EasyScan_HEP
============

:EasyScan_HEP: A tool for connecting programs to scan the parameter space of physics models
:Author: Yang Zhang, Liangliang Shang, Yang Xiao, Yuanfang Yue
:Version: 2.0
:GitHub: https://github.com/phyzhangyang/EasyScan_HEP
:Website: https://easyscanhep.hepforge.org
:Documentation: https://arxiv.org/pdf/2304.03636.pdf


Installation instructions
-------------------------

EasyScan_HEP is a Python3 code. The core scan workflow depends on
*numpy*, *scipy* and *ConfigParser*. The plot and result inspection
functions further require *matplotlib* and *pandas*. The MultiNest sampler
requires the optional *pymultinest* library.

On Ubuntu or Debian systems, install the basic system packages first::

    sudo apt install python3-pip python3-venv python3-tk

Install the Python dependencies into the Python environment you will use
to run EasyScan_HEP::

    python3 -m pip install numpy scipy matplotlib ConfigParser pandas

If you prefer an isolated environment, you can create a local virtual
environment first and then install the same dependencies there::

    python3 -m venv .venv
    .venv/bin/python -m pip install numpy scipy matplotlib ConfigParser pandas

Install *pymultinest* only if the MultiNest sampler is needed::

    pip install pymultinest

Install *dynesty* only if the Dynesty nested sampler is needed::

    pip install dynesty

Install *emcee* only if the EMCEE ensemble MCMC sampler is needed::

    pip install emcee

The local Web UI has additional dependencies: *fastapi*, *uvicorn*,
*jinja2* and *python-multipart*. They can be installed with::

    python3 -m pip install fastapi uvicorn jinja2 python-multipart

The "easyscan.py" in folder "bin" is the main program, which is executed with configuration file through the command line,
::

    ./bin/easyscan.py templates/example_random.ini

Check a configuration file without running a scan with::

    ./bin/easyscan.py --check templates/example_random.ini

Here *example_random.ini* is an example configuration file provided in EasyScan_HEP. It performs a scan on a simplified model,
::

    f(x,y) = sin^2 x + cos^2 y,
    
using random sampler, where *x* and *y* are input parameters in range *[0,\pi]* and *[-\pi,\pi]*, respectively, and *f* is output parameter. 

Other example configuration files in *templates* folder, including *example_grid.ini*, *example_bestfit.ini*, *example_mcmc.ini*, *example_emcee.ini*, *example_multinest.ini* and *example_dynesty.ini*, exhibit briefly usages of other samplers in EasyScan_HEP.

Configuration file *templates/scan_MSSM_for_mW.ini* is an simply physical examples. Relevant programs need to be installed beforehand, using
::

    cd utils/MSSM_mW
    make
    
and then it can be executed with 
::

    ./bin/easyscan.py templates/scan_MSSM_for_mW.ini

Local Web UI
------------

EasyScan_HEP also provides a local single-user Web UI. Start it with the
same main entrypoint::

    ./bin/easyscan.py -ui

The UI runs at ``http://127.0.0.1:8000/`` and opens the browser
automatically when possible. Keep the terminal open while using the UI.
Press ``Control-C`` in that terminal to stop the server.

The Web UI can build ``.ini`` configuration files visually, choose local
paths from the browser UI, run EasyScan in the background, show live logs,
stop running jobs, inspect result files and generated plots, and check
configuration consistency before a run. Generated UI ``.ini`` and ``.log``
artifacts are stored in the directory where the UI was opened. The
Configuration File panel can also call an OpenAI-compatible large language
model API to turn a natural-language request into a checked ``.ini`` setup.
