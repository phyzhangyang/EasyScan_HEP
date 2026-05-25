---
name: easyscan-hep
description: >
  Use when the user mentions HEP parameter scans, scanning, sampling, samplers,
  MCMC, EMCEE, Dynesty, MultiNest, BESTFIT, nested sampling, likelihood fitting,
  exclusion scans, constraint scans, or the Chinese terms 扫描, 采样, 参数空间.
  Generate EasyScan_HEP .ini configuration files first, check them when possible,
  and do not run scans until the user explicitly confirms.
---

# EasyScan_HEP Skill

EasyScan_HEP connects external physics programs to scan high-energy-physics
parameter spaces. Your default job is to generate the `.ini` configuration file
that EasyScan_HEP needs. Stop after producing the `.ini` and any non-running
config check; only run a scan when the user later asks clearly to run it.

## Core Workflow

1. **Detect or install EasyScan_HEP**
   - First check whether the current project or a parent directory is an EasyScan_HEP source tree. Look for `easyscan_hep/cli.py`, `bin/easyscan.py`, `pyproject.toml`, or `setup.py`.
   - If a source tree is present, prefer running checks with `python3 -m easyscan_hep.cli --check <config.ini> --json` from that tree; `./bin/easyscan.py --check <config.ini> --json` is acceptable for source-tree compatibility.
   - If no source tree is present, check system installation with `command -v easyscan`, `python3 -m pip show easyscan-hep`, and `python3 -c "import easyscan_hep"`.
   - If EasyScan_HEP is missing, install it with `python3 -m pip install easyscan-hep`.
   - If the package imports but `easyscan` is not on `PATH`, use `python3 -m easyscan_hep.cli` as the fallback command.

2. **Collect scan requirements**
   - Required information: scan purpose/model, scan method, input parameters with ranges or fixed values, scan size (`Number of points`, GRID point/interval counts, optimizer iterations, or live points), result folder, connected program command, input/output file mappings when EasyScan must write/read files, output variables, constraints, and requested plots.
   - If important information is missing, do not invent physics details. Ask the user for the missing fields in a compact checklist.
   - If the user only wants a quick EasyScan_HEP test, generate a minimal built-in `TestFunction` style config.

3. **Generate the EasyScan_HEP `.ini`**
   - Load `references/ini_config_reference.md` when writing config syntax.
   - Load `references/scan_methods.md` when selecting or explaining scan methods.
   - Load `references/external_programs.md` only when the request involves external physics programs such as SUSY-HIT, HiggsBounds, FeynHiggs, GM2Calc, micrOMEGAs, or a user-provided executable.
   - Put the generated `.ini` in the directory where `easyscan` should be launched unless the user asks for another location.
   - Preserve user-provided paths and names. Prefer relative paths from the directory where `easyscan` will be launched, or `~` paths. Absolute paths are acceptable when the user needs them, but parallel mode still requires relative paths for `Parallel folder`, `Command path`, input files, and output files.
   - For external Python programs, prefer an explicit interpreter that matches the intended EasyScan runtime, such as `/usr/bin/python3 TestFunction.py` or `.venv/bin/python run.py`, instead of relying on `./script.py` and its shebang.

4. **Check, then stop**
   - Run `--check --json` if EasyScan_HEP is available and checking is feasible. This is allowed because it validates the config without running the scan.
   - Run `--check` from the intended launch directory, because EasyScan_HEP resolves config paths relative to the current working directory.
   - Show the generated `.ini` path, summarize the scan method, parameters, program mappings, constraints, plots, and report the check result.
   - After generating the `.ini` and reporting the check result, ask the user whether they want to use the EasyScan_HEP UI `Check Config` workflow to inspect the same file.
   - Do not execute the scan in the same turn unless the user explicitly asked to run after reviewing the `.ini`.

## Missing Information Checklist

Ask for only the fields needed for the requested scan:

- Scan method: `RANDOM`, `GRID`, `BESTFIT`, `MCMC`, `EMCEE`, `DYNESTY`, `MULTINEST`, `ONEPOINT`, or `ONEPOINTBATCH`.
- Input parameters: names and ranges for scanned parameters; fixed values for fixed parameters; `GRID` point counts or interval counts; interval divisor and initial value for `MCMC`/`EMCEE`.
- Scan size: total points for `RANDOM`, saved/live/surviving target as appropriate for `MCMC`, `EMCEE`, `DYNESTY`, `MULTINEST`, or max iterations for `BESTFIT`.
- Program connection: execute command, command directory, input file(s), how variables are written, output file(s), and how outputs are read. For Python programs, ask which Python interpreter should run the external code when dependencies may differ. For parallel runs, ask for an existing relative `Parallel folder` that prefixes command paths, often `utils`.
- Constraints: Gaussian values, bounds, or free-form chi-square definitions.
- Output: result folder and requested plot types.

## Scan Method Defaults

- Use `RANDOM` for quick exploration when the user gives ranges but no fitting goal.
- Use `GRID` only for low-dimensional regular grids. EasyScan_HEP's GRID column is the number of intervals (`BinNum`), so if the user asks for `N` grid points including endpoints, write `N-1`.
- Use `BESTFIT` for finding a best-fit point with SciPy differential evolution.
- Use `MCMC` or `EMCEE` for likelihood/posterior sampling; require constraints plus initial values and interval divisors for each scanned parameter. If the user gives a desired proposal step size, convert it to `interval = (max - min) / step`.
- Use `DYNESTY` or `MULTINEST` for nested sampling; mention optional dependencies (`dynesty` or `pymultinest` plus MultiNest system library).
- Use `ONEPOINT`/`ONEPOINTBATCH` for fixed point evaluation, not scans over priors.

## Hard Rules

- Generate EasyScan_HEP `.ini` files; do not write new scan-driver code when an `.ini` can express the task.
- Do not run a scan before the user has had a chance to inspect the `.ini`.
- Do not delete or overwrite result directories without explicit user approval.
- Prefer the installed `easyscan` command; fallback to `python3 -m easyscan_hep.cli` when PATH is not configured.
- Use `--check --json` before suggesting a config is ready whenever the checker is available.
- Keep generated paths relative to the launch directory or use `~` paths unless there is a strong reason otherwise.
- For Python external programs, avoid ambiguous shebang execution; write the interpreter explicitly when environment consistency matters.

## References

- `references/ini_config_reference.md`: current `.ini` syntax and field rules.
- `references/scan_methods.md`: sampler selection and method-specific requirements.
- `references/external_programs.md`: external program connection examples.
