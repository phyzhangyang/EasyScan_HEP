============
EasyScan_HEP
============

:EasyScan_HEP: A tool for connecting programs to scan the parameter space of physics models
:Author: Yang Zhang, Liangliang Shang, Yang Xiao
:Version: 2.0
:GitHub: https://github.com/phyzhangyang/EasyScan_HEP
:Website: https://easyscanhep.hepforge.org
:Documentation: https://arxiv.org/pdf/2304.03636.pdf


Installation instructions
-------------------------

EasyScan_HEP is a Python3 code. The core scan workflow depends on
*numpy*, *scipy*, *matplotlib* and *pandas*. These core Python
dependencies are installed automatically when EasyScan_HEP is installed
with pip. The MultiNest sampler additionally requires the optional
*pymultinest* Python package and a working MultiNest library.

On Ubuntu or Debian systems, install the basic system packages first::

    sudo apt install python3-pip python3-venv python3-tk

Install EasyScan_HEP itself from the source directory with::

    python3 -m pip install .

After installation, the command ``easyscan`` is available from any
directory. Use ``easyscan`` rather than ``./easyscan``; the ``./`` prefix
means "run a file in the current directory" and is only appropriate for
local scripts. Relative paths inside a configuration file are resolved
from the directory where ``easyscan`` is launched, so run it from your
project directory or use absolute paths in the configuration file.

If you prefer an isolated environment, you can create a local virtual
environment first and then install EasyScan_HEP there::

    python3 -m venv .venv
    .venv/bin/python -m pip install .

Install the local Web UI dependencies with::

    python3 -m pip install ".[ui]"

Install *dynesty* only if the Dynesty nested sampler is needed::

    python3 -m pip install ".[dynesty]"

Install *emcee* only if the EMCEE ensemble MCMC sampler is needed::

    python3 -m pip install ".[emcee]"

Install *pymultinest* only if the MultiNest sampler is needed::

    python3 -m pip install ".[multinest]"

The Python-only optional features can be installed together with::

    python3 -m pip install ".[all]"

The ``all`` extra includes the Web UI, Dynesty and EMCEE dependencies.
MultiNest is kept separate because it may require additional native
library setup outside pip.

The installed ``easyscan`` command is executed with a configuration file through the command line,
::

    easyscan templates/example_random.ini

The source-tree entry point remains available for development use::

    ./bin/easyscan.py templates/example_random.ini

Check a configuration file without running a scan with::

    easyscan --check templates/example_random.ini

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

    easyscan templates/scan_MSSM_for_mW.ini

Local Web UI
------------

EasyScan_HEP also provides a local single-user Web UI. Start it with the
same main entrypoint::

    easyscan -ui

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

Package release check
---------------------

Before publishing a release package, install the packaging tools once::

    python3 -m pip install build twine

Then build and smoke-test the source distribution and wheel with::

    python3 scripts/check_package.py

The check builds EasyScan_HEP in a temporary directory, verifies the
distribution metadata with ``twine check``, installs the wheel into an
isolated target directory, and runs ``easyscan --check`` on the built-in
test configuration.
