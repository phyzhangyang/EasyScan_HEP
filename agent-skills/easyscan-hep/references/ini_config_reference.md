# EasyScan_HEP INI Configuration Reference

Use this reference when generating the `.ini` file that EasyScan_HEP needs.
Prefer concise, explicit configs. Do not run a scan before the user has checked
the generated file.

## File Structure

An EasyScan_HEP config usually contains:

```ini
[scan]
[program1]
[program2]
[constraint]
[plot]
```

`[scan]` and at least one `[programN]` section are required for ordinary scans.
`[constraint]` is required for likelihood-driven methods such as `BESTFIT`,
`MCMC`, `EMCEE`, `DYNESTY`, and `MULTINEST`; use `Gaussian` or
`FreeFormChi2`. `[plot]` is optional.

## `[scan]` Section

Required keys:

| Key | Meaning |
|---|---|
| `Result folder name` | Output directory. Do not include spaces. |
| `Scan method` | `RANDOM`, `GRID`, `BESTFIT`, `MCMC`, `EMCEE`, `DYNESTY`, `MULTINEST`, `ONEPOINT`, or a batch form such as `ONEPOINT/path/to/file` or `ONEPOINTBATCH/path/to/file`. |
| `Input parameters` | Multi-line parameter definitions. Format depends on scan method. |

Common optional keys:

| Key | Used by | Meaning |
|---|---|---|
| `Number of points` | `RANDOM`, `BESTFIT`, `MCMC`, `EMCEE`, `DYNESTY`, `MULTINEST` | Total points, optimizer iterations, saved samples, or live points depending on method. EasyScan defaults to `10` if omitted, but generated configs should set it explicitly. |
| `Random seed` | `RANDOM`, `MCMC`, `EMCEE`, `MULTINEST` | Reproducibility seed. |
| `Parallel threads` | `RANDOM`, `GRID`, `MCMC`, `ONEPOINTBATCH`, `MULTINEST` | Number of processes. Avoid adding it for other methods unless the implementation/checker changes. |
| `Parallel folder` | when `Parallel threads > 1` | Existing working directory used for parallel runs. It must prefix each `Command path`; prefer a top-level relative folder such as `utils`. EasyScan copies it to `p0_<folder>`, `p1_<folder>`, etc. |
| `Interval of print` | all methods | Progress print interval. |
| `Acceptance rate` | `MCMC` | Target acceptance rate, often `0.25`. |
| `Resume mode` | resume workflows | Continue from an existing result folder. |
| `Debug mode` | debugging | More verbose execution. |

### Input Parameter Formats

`RANDOM`, `BESTFIT`, `DYNESTY`, `MULTINEST`:

```ini
Input parameters:  x, flat, 0, 3.14
                   y, log,  1e-3, 1e3
                   z, fixed, 1.0
```

`GRID`:

```ini
Input parameters:  x, flat, 0, 3.14, 20
                   y, flat, -3.14, 3.14, 30
```

The fifth column is EasyScan_HEP's `BinNum`, meaning interval count. The actual
number of grid points is `BinNum + 1` because both endpoints are included. If a
user asks for `N` grid points in a range, write `N-1` in the config. Example:
10 points from `0` to `pi` should be written as `0, 3.141592653589793, 9`.

`MCMC` and `EMCEE`:

```ini
Input parameters:  x, flat, 0, 3.14, 10, 1.5
                   y, flat, -3.14, 3.14, 20, 0.0
```

The fifth column is an interval divisor, not the proposal step size itself.
EasyScan uses `step = (max - min) / interval`. If the user asks for a step
size `s`, write `(max - min) / s` as the interval. The sixth column is the
initial value. For `log` priors, minimum, maximum, and initial values must be
positive.

`ONEPOINT`:

```ini
Input parameters:  x, fixed, 1.5
                   y, fixed, 0.0
```

`ONEPOINTBATCH`:

```ini
Scan method:       ONEPOINTBATCH/path/to/points.dat
Input parameters:  x, fixed, 0
                   y, fixed, 0
```

The batch file supplies the actual point values.

Prior types:

| Prior | Meaning |
|---|---|
| `flat` | Uniform scan between minimum and maximum. |
| `log` | Logarithmic scan between positive minimum and maximum values. |
| `fixed` or `Fixed` | Fixed value. |

## `[programN]` Sections

Program sections connect external executables. Number them in order:
`[program1]`, `[program2]`, etc.

Recommended keys for most connected programs:

| Key | Meaning |
|---|---|
| `Program name` | Human-readable program label. Optional; defaults to the section name such as `program1`. |
| `Execute command` | Command run inside `Command path`, e.g. `./TestFunction.py` or `python3 run.py`. |
| `Command path` | Directory from which the command is executed. Optional; defaults to the launch directory. |
| `Input file` | File ID and input path, e.g. `1, utils/input.dat`. Needed only when EasyScan writes inputs. |
| `Input variable` | How scan variables are written to input files. Must be paired with `Input file`. |
| `Output file` | File ID and output path. Needed when EasyScan reads outputs. |
| `Output variable` | How EasyScan_HEP reads outputs. Must be paired with `Output file`. |

`Execute command` is required. In ordinary scans, include all file and variable
mappings for clarity, even though the parser only loads input/output mappings
when both members of each pair are present.

Optional keys:

| Key | Meaning |
|---|---|
| `Bound` | Program-level acceptance or boundary conditions. |
| `Time limit in minute` | Timeout; when set, EasyScan_HEP uses `subprocess.popen`. |
| `Clean output file` | Whether to clean output before execution. |
| `Command executor` | `os.system` or `subprocess.popen`; not needed when time limit is set. |

Use `Command executor: os.system` or `Command executor: subprocess.popen`.
Invalid values fall back to `os.system`. If `Time limit in minute` is set,
EasyScan forces `subprocess.popen`.

For external Python scripts, prefer an explicit interpreter in `Execute command`
that matches the intended environment, such as `/usr/bin/python3 TestFunction.py`
or `.venv/bin/python run.py`. `./script.py` follows the script shebang and may
use a different Python with different installed packages.

Prefer paths relative to the directory where `easyscan` will be launched, or
`~` paths. Absolute paths are supported for ordinary serial runs. In parallel
mode, keep `Parallel folder`, `Command path`, input files, and output files
relative to the launch directory, because EasyScan copies the parallel folder
to worker directories.

For agent workflows, check the file with machine-readable output:

```bash
easyscan --check config.ini --json
```

After the user confirms the configuration, non-interactive runs should choose
an explicit result-folder policy:

```bash
easyscan --run config.ini --overwrite stop --json
```

### File Declarations

```ini
Input file:      1, utils/TestFunction_input.dat
Output file:     1, utils/TestFunction_output.dat
```

The numeric file ID is referenced by `Input variable` and `Output variable`.

### Input/Output Variable Mapping

Position:

```ini
Input variable:  x, 1, Position, 1, 1
Output variable: f, 1, Position, 1, 2
```

Label:

```ini
Output variable: mh, 1, Label, "mh =", 2
```

SLHA:

```ini
Input variable:  mu,  1, SLHA, BLOCK, EXTPAR, 31
Output variable: mh1, 1, SLHA, BLOCK, MASS, 25
```

Replace:

```ini
Input variable:  tanb, 1, Replace, "tanb"
```

File:

```ini
Input variable:  input_card, 1, File, SAVE
```

### Bounds

Simple inequality:

```ini
Bound: x, >=, 0
       y, <=, 10
```

Range:

```ini
Bound: x, 0, 3.14
```

Boundary curve:

```ini
Bound: y, x, MIN, templates/bound.txt
```

## `[constraint]` Section

Use constraints when a method needs likelihood or chi-square information.

Gaussian:

```ini
[constraint]
Gaussian: f, 1.0, 0.2
          mW, 80.4, 0.01, symm, MW_constraint
```

Format:

```ini
Gaussian: <variable_or_expression>, <mean>, <uncertainty>, [symm|upper|lower], [name]
```

Free-form chi-square from an output variable:

```ini
[constraint]
FreeFormChi2: chi2_total
```

Mathematical expressions can be used when variable names are already available:

```ini
Gaussian: sin(2*alpha), 1.0, 0.1
```

## `[plot]` Section

Supported plot rows:

```ini
[plot]
Histogram: f, hist_f
Scatter:   x, y, scatter_xy
Color:     x, y, f, color_xy_f
Contour:   x, y, f, contour_xy_f
```

Formats:

| Plot | Format |
|---|---|
| `Histogram` | `x, figure_name` |
| `Scatter` | `x, y, figure_name` |
| `Color` | `x, y, value, figure_name` |
| `Contour` | `x, y, value, figure_name` |

## Minimal TestFunction Example

Use this when the user wants a simple EasyScan_HEP test and has not supplied an
external program.

```ini
[scan]
Result folder name:  example_random
Scan method:       RANDOM
Input parameters:  x, flat, 0, 3.14
                   y, flat, -3.14, 3.14
Number of points:  100
Parallel threads:  1
Interval of print: 1

[program1]
Program name:    TestFunction
Execute command: /usr/bin/python3 TestFunction.py
Command path:    utils/
Input file:      1, utils/TestFunction_input.dat
Input variable:  x, 1, Position, 1, 1
                 y, 1, Position, 1, 2
Output file:     1, utils/TestFunction_output.dat
Output variable: f, 1, Position, 1, 2

[constraint]
Gaussian: f, 1, 0.2

[plot]
Color: x, y, f, test_run
```

## Generation Rules

- Use uppercase scan method names in generated configs.
- Include only fields used by the selected scan method unless the user asks for future-edit placeholders.
- Preserve disabled/not-used values only when editing an existing user config; do not add unused fields to a new config.
- Put generated configs in the intended launch directory when possible.
- Keep relative paths relative to the directory where `easyscan` will be launched.
- Avoid macOS-style `/Users/...` absolute paths in generated runtime config unless the user explicitly asks; use relative paths or `~`.
- For GRID, if the user asks for `N` endpoint-including points, write `N-1` as the fifth column.
- For MCMC/EMCEE, if the user gives a proposal step size `s`, write `(max - min) / s` as the fifth column.
- Avoid reserved variable names such as `pi`, `sin`, `cos`, scan method names, `dwell`, `probability`, and `-2lnlike`.
- If `Parallel threads` is larger than `1`, include an existing `Parallel folder`.
- Before telling the user a config is ready, run `easyscan --check <file>` or `python3 -m easyscan_hep.cli --check <file>` from the intended launch directory when available.
