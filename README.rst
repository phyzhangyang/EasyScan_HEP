=======
EasyScan_HEP
=======

:EasyScan_HEP: A tool for easily connecting programs to scan parameter space of high energy physics models
:Author: EasyScan_HEP collaboration
:Organiser: Yang Zhang, Liangliang Shang
:Version: 0.1.0
:GitHub: https://github.com/phyzhangyang/EasyScan_HEP
:Website: https://easyscanhep.hepforge.org
:Documentation: https://arxiv.org/pdf/190X.XXXXX.pdf


Installation instructions
-------------------------

EasyScan_HEP is a Python3 codes with dependencies on numpy, scipy and ConfigParser libraries. The optional plot functions and MultiNest sampler further require matplotlib, pandas and pymultinest libraries. The dependencies can be installed via pip:
.. code:: bash
    sudo apt install python3-pip python3-tk
    sudo pip3 install numpy scipy matplotlib ConfigParser pandas pymultinest

The "easyscan.py" in folder "bin" is the main program, which is executed with
configuration file through the command line,
.. code:: bash
    ./bin/easyscan.py templates/example_random.ini


