---
name: easyscan-hep
description: >
  Use when the user mentions HEP parameter scans, scanning, sampling, samplers,
  MCMC, EMCEE, Dynesty, MultiNest, BESTFIT, nested sampling, likelihood fitting,
  exclusion scans, constraint scans, or the Chinese terms 扫描, 采样, 参数空间.
  Generate EasyScan_HEP .ini configuration files first, check them when possible,
  and do not run scans until the user explicitly confirms. If scan details are
  missing, ask for the minimum required fields before generating a config.
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
   - If the user only invokes this skill, says only "帮我做一个扫描", "do a scan", "sampling", or gives similarly vague intent with no concrete scan details, stop after discovery/installation checks and ask for the minimum fields listed below. Do not choose a default scan method, physics model, parameter range, external program, output variable, constraint, plot, result folder, or point count.
   - If the user provides partial details, create an `.ini` only when enough information exists to express a valid EasyScan_HEP config. Otherwise ask for the missing fields. Do not read external program source, binaries, input cards, output files, existing `.ini` files, templates, examples, or prior generated configs to infer input/output mappings unless the user explicitly asks you to inspect those files for mappings.
   - If the request names a user-provided executable or script but does not provide `Input file`, `Input variable`, `Output file`, and `Output variable` mappings, stop and ask for those mappings before writing the `.ini`. A program path or command alone is not enough to infer file mappings, even for `utils/TestFunction.py`.
   - Treat every executable or script path in the user's request as user-provided unless the user explicitly says it is the built-in demo. Names such as `TestFunction.py` do not automatically imply the bundled demo mapping.
   - If important information is missing, do not invent physics details. Ask the user for the missing fields in a compact checklist.
   - Generate a minimal built-in `TestFunction` style config only when the user explicitly asks for a quick EasyScan_HEP test or demo and accepts the bundled demo I/O mapping.

3. **Generate the EasyScan_HEP `.ini`**
   - Before writing any `.ini`, perform this gate:
     - If the config contains `[programN]` with file mappings, each `Input file`, `Input variable`, `Output file`, and `Output variable` must be either directly supplied by the user in the current request/conversation or explicitly authorized by the user from a named example.
     - If mappings are missing or only found in examples/templates, do not write the `.ini`; ask for the missing mappings.
     - If the user asks only for the `[scan]` section or a skeleton config, write no `[programN]` mappings and clearly label the file as incomplete/uncheckable.
   - Load `references/ini_config_reference.md` when writing config syntax. Use examples only for syntax and field formatting, never as evidence for the current task's input/output mappings.
   - Load `references/scan_methods.md` when selecting or explaining scan methods.
   - Load `references/external_programs.md` only when the request involves external physics programs such as SUSY-HIT, HiggsBounds, FeynHiggs, GM2Calc, micrOMEGAs, or a user-provided executable.
   - Put the generated `.ini` in the directory where `easyscan` should be launched unless the user asks for another location.
   - Preserve user-provided paths and names. Prefer relative paths from the directory where `easyscan` will be launched, or `~` paths. Absolute paths are acceptable when the user needs them, but parallel mode still requires relative paths for `Parallel folder`, `Command path`, input files, and output files.
   - For external Python programs, prefer an explicit interpreter that matches the intended EasyScan runtime, such as `/usr/bin/python3 TestFunction.py` or `.venv/bin/python run.py`, instead of relying on `./script.py` and its shebang.

4. **Check, then stop**
   - Run `--check --json` if EasyScan_HEP is available and checking is feasible. This is allowed because it validates the config without running the scan.
   - Do not add guessed mappings, constraints, plots, or other fields solely to make `--check` pass. If checking is not feasible because required mappings were not supplied, report that the config was not generated and ask for the missing mappings.
   - Run `--check` from the intended launch directory, because EasyScan_HEP resolves config paths relative to the current working directory.
   - Show the generated `.ini` path, summarize the scan method, parameters, program mappings, constraints, plots, and report the check result.
   - After generating the `.ini` and reporting the check result, ask the user whether they want to use the EasyScan_HEP UI `Check Config` workflow to inspect the same file.
   - Do not execute the scan in the same turn unless the user explicitly asked to run after reviewing the `.ini`.

5. **Open the UI with the generated config when requested**
   - When the user agrees to inspect the generated file in the UI, start EasyScan_HEP UI with the config path argument so the UI loads that `.ini` immediately: `./bin/easyscan.py -ui <config.ini>` from the intended launch directory, or `/usr/bin/python3 bin/easyscan.py -ui <config.ini>` when using an explicit Python interpreter.
   - Prefer a path relative to the UI launch directory, such as `mytest.ini`, when the `.ini` is in that directory. Absolute paths are acceptable if needed.
   - Do not start the UI as plain `./bin/easyscan.py -ui` after generating a specific config; that opens the default/example state and can mislead the user.
   - If the UI server is already running, open or navigate the browser to `http://127.0.0.1:8000/?config=<URL-encoded-config-path>` for the generated `.ini`, or restart the UI with the config argument if the running instance cannot load the requested file.
   - After opening the UI, stop and ask the user to inspect the loaded configuration in the browser. Point them to verify the result folder, scan method, input ranges, program command/path, input/output mappings, constraints, and plots as applicable.
   - Do not automatically click UI actions such as `Check Config`, `Export to INI`, or any run-related button after the UI is opened unless the user explicitly asks for that specific UI action. If the user asks for the UI `Check Config` workflow, click only `Check Config`, report the result, and then stop for user review.
   - Opening the UI and using its `Check Config` workflow is still not permission to run the scan. Do not click `Generate and Run` or start a run unless the user explicitly asks.

## Missing Information Checklist

Ask for only the fields needed for the requested scan:

- Scan method: `RANDOM`, `GRID`, `BESTFIT`, `MCMC`, `EMCEE`, `DYNESTY`, `MULTINEST`, `ONEPOINT`, or `ONEPOINTBATCH`.
- Input parameters: names and ranges for scanned parameters; fixed values for fixed parameters; `GRID` point counts or interval counts; interval divisor and initial value for `MCMC`/`EMCEE`.
- Scan size: total points for `RANDOM`, saved/live/surviving target as appropriate for `MCMC`, `EMCEE`, `DYNESTY`, `MULTINEST`, or max iterations for `BESTFIT`.
- Program connection: execute command, command directory, input file(s), how variables are written, output file(s), and how outputs are read. For Python programs, ask which Python interpreter should run the external code when dependencies may differ. For parallel runs, ask for an existing relative `Parallel folder` that prefixes command paths, often `utils`.
- Constraints: Gaussian values, bounds, or free-form chi-square definitions.
- Output: result folder and requested plot types.

For a vague request such as "帮我做一个扫描", ask for this minimum set before writing any `.ini`:

- Scan method and scan size.
- Input parameters with ranges or fixed values.
- Result folder.
- External program command and command directory, unless the user explicitly wants the built-in `TestFunction` example.
- Input-file mapping for each scanned parameter when EasyScan_HEP must write values.
- Output file and output variable definitions when EasyScan_HEP must read results.
- Constraints or objective definition for likelihood, nested-sampling, or best-fit methods.
- Requested plots, or confirmation that no plots are needed.

For a request that provides scan ranges and a program command but omits I/O
mappings, ask only for the missing program contract before writing the config:

- Input file path(s) that EasyScan_HEP should edit.
- How each scanned parameter is written: `Position`, `Label`, `Replace`,
  `SLHA`, or `File`, with the needed row/column/label/block details.
- Output file path(s) that EasyScan_HEP should read.
- Output variable name(s) and how each is read.

## Scan Method Defaults

- Use `RANDOM` for quick exploration when the user gives ranges but no fitting goal.
- Use `GRID` only for low-dimensional regular grids. EasyScan_HEP's GRID column is the number of intervals (`BinNum`), so if the user asks for `N` grid points including endpoints, write `N-1`.
- Use `BESTFIT` for finding a best-fit point with SciPy differential evolution.
- Use `MCMC` or `EMCEE` for likelihood/posterior sampling; require constraints plus initial values and interval divisors for each scanned parameter. If the user gives a desired proposal step size, convert it to `interval = (max - min) / step`.
- Use `DYNESTY` or `MULTINEST` for nested sampling; mention optional dependencies (`dynesty` or `pymultinest` plus MultiNest system library).
- Use `ONEPOINT`/`ONEPOINTBATCH` for fixed point evaluation, not scans over priors.

## Hard Rules

- Generate EasyScan_HEP `.ini` files; do not write new scan-driver code when an `.ini` can express the task.
- Do not create a default `.ini` from an empty or vague scan request. Ask for the minimum fields first.
- Do not silently invent output variables, constraints, plots, or program mappings just to satisfy the checker.
- Do not proactively inspect external program files, existing configs, templates, examples, or prior generated configs to guess `Input file`, `Input variable`, `Output file`, or `Output variable`. Ask the user for these mappings. File existence checks such as `ls` or `test -e` are acceptable; content inspection is allowed only after an explicit user request.
- Existing examples may be used to learn syntax, but must not be copied as the current task's program mapping unless the user explicitly says to use that example mapping.
- When assumptions come from user-provided mappings, report them explicitly before writing the `.ini`.
- Do not run a scan before the user has had a chance to inspect the `.ini`.
- Do not delete or overwrite result directories without explicit user approval.
- Prefer the installed `easyscan` command; fallback to `python3 -m easyscan_hep.cli` when PATH is not configured.
- Use `--check --json` before suggesting a config is ready whenever the checker is available.
- When opening the UI for a generated config, pass the config path to `-ui`; do not open the default UI state for config inspection.
- Keep generated paths relative to the launch directory or use `~` paths unless there is a strong reason otherwise.
- For Python external programs, avoid ambiguous shebang execution; write the interpreter explicitly when environment consistency matters.

## References

- `references/ini_config_reference.md`: current `.ini` syntax and field rules.
- `references/scan_methods.md`: sampler selection and method-specific requirements.
- `references/external_programs.md`: external program connection examples.
