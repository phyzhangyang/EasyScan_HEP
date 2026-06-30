============
EasyScan_HEP
============

:EasyScan_HEP 2: Agent-Ready Parameter Scans for High-Energy Physics
:Author: Yang Zhang, Liangliang Shang, Yang Xiao, Yuanfang Yue
:Version: 2.0
:GitHub: https://github.com/phyzhangyang/EasyScan_HEP
:Website: https://easyscanhep.hepforge.org
:Documentation: https://arxiv.org/pdf/2304.03636.pdf  https://arxiv.org/pdf/2607.xxxxx.pdf


Overview
--------

EasyScan_HEP connects external physics programs to parameter scans through an
explicit ``.ini`` configuration file. EasyScan_HEP 2 keeps this configuration
workflow and adds an installable ``easyscan`` command, a dry-run checker,
structured run/result interfaces, a local Web UI, an agent skill, and the
``BESTFIT``, ``EMCEE`` and ``DYNESTY`` scan methods.


Installation
------------

Install from GitHub::

    python3 -m pip install git+https://github.com/phyzhangyang/EasyScan_HEP.git

or from a local source tree::

    python3 -m pip install .

Optional features can be installed from a local source tree with::

    python3 -m pip install ".[ui]"
    python3 -m pip install ".[emcee]"
    python3 -m pip install ".[dynesty]"
    python3 -m pip install ".[all]"

``MULTINEST`` additionally needs ``pymultinest`` and a working MultiNest
library.


Command-line use
----------------

Run a configuration file::

    easyscan templates/example_random.ini

Check a configuration without running scan points::

    easyscan --check templates/example_random.ini
    easyscan --check templates/example_random.ini --json

Run non-interactively and choose the overwrite policy explicitly::

    easyscan --run templates/example_random.ini --overwrite stop --json

Summarize an existing result directory::

    easyscan --results example_random --json

Relative paths in a configuration file are resolved from the directory where
``easyscan`` is launched.


Local Web UI
------------

Start the local single-user UI with::

    easyscan --ui

or preload a configuration file::

    easyscan --ui templates/example_random.ini

The UI edits the same ``.ini`` configuration used by the command line. It can
load, edit, preview, check, save and run configurations, show logs and output
files, and enable or disable method-dependent controls according to the selected
scan method.


Agent skill
-----------

The agent skill is maintained separately from the EasyScan_HEP package:

    https://github.com/Contract-Mediated-Agent/easyscan-skill.git

For Codex-style skill folders, install it with::

    mkdir -p ~/.codex/skills
    git clone https://github.com/Contract-Mediated-Agent/easyscan-skill.git ~/.codex/skills/easyscan-hep

The skill is ``.ini``-first: it should prepare the configuration, run
``easyscan --check --json`` when possible, and wait for user review before
launching a scan. The checker validates syntax and operational consistency; it
does not validate the physics model or constraints.


Examples
--------

Basic examples are provided in ``templates/``, including ``example_random.ini``,
``example_grid.ini``, ``example_bestfit.ini``, ``example_mcmc.ini``,
``example_emcee.ini``, ``example_multinest.ini`` and ``example_dynesty.ini``.

The MSSM ``m_W`` example uses external programs prepared in ``utils/MSSM_mW``::

    cd utils/MSSM_mW
    make
    cd ../..
    easyscan templates/scan_MSSM_for_mW.ini

The scalar-singlet dark-matter/electroweak-phase-transition example prepares
fixed external-program versions, ``micrOMEGAs 7.1`` and ``PhaseTracer 2.2.0``::

    python3 utils/SSM_DM_EWPT/bootstrap.py
    easyscan templates/scan_SSM_DM_EWPT.ini


Development check
-----------------

Build and smoke-test the package with::

    python3 scripts/check_package.py
