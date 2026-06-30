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
explicit ``.ini`` configuration file. EasyScan_HEP 2 adds an installable
``easyscan`` command, a configuration checker, structured run/result interfaces,
a local Web UI, an agent skill, and the ``BESTFIT``, ``EMCEE`` and ``DYNESTY``
scan methods.


Install
-------

Install from GitHub or from a local source tree::

    python3 -m pip install git+https://github.com/phyzhangyang/EasyScan_HEP.git
    python3 -m pip install .

Optional local extras::

    python3 -m pip install ".[ui]"
    python3 -m pip install ".[all]"

``MULTINEST`` additionally requires ``pymultinest`` and a working MultiNest
library.


Basic Commands
--------------

::

    easyscan templates/example_random.ini
    easyscan --check templates/example_random.ini --json
    easyscan --run templates/example_random.ini --overwrite stop --json
    easyscan --results example_random --json
    easyscan --ui templates/example_random.ini

Relative paths in ``.ini`` files are resolved from the directory where
``easyscan`` is launched.


Local Web UI
------------

The Web UI edits the same ``.ini`` files used by the command line. It can load,
edit, preview, check, save and run configurations, show logs and result files,
and mark method-dependent controls according to the selected scan method.


Agent Skill
-----------

The agent skill is maintained separately:

    https://github.com/Contract-Mediated-Agent/easyscan-skill.git

Install it for Codex-style skill folders with::

    mkdir -p ~/.codex/skills
    git clone https://github.com/Contract-Mediated-Agent/easyscan-skill.git ~/.codex/skills/easyscan-hep

The skill is ``.ini``-first: it prepares a configuration, runs
``easyscan --check --json`` when possible, and waits for user review before
launching a scan. The checker validates syntax and operational consistency, not
the underlying physics model.


Examples
--------

The ``templates/`` directory contains basic examples for ``RANDOM``, ``GRID``,
``BESTFIT``, ``MCMC``, ``EMCEE``, ``MULTINEST`` and ``DYNESTY``.

Physical examples include::

    python3 utils/SSM_DM_EWPT/bootstrap.py
    easyscan templates/scan_SSM_DM_EWPT.ini


