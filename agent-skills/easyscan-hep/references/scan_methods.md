# EasyScan_HEP Scan Methods

Use this file when choosing or explaining a sampler for an EasyScan_HEP
configuration. Method names are case-insensitive in user input, but write them
in uppercase in generated `.ini` files.

## Method Comparison

| Method | Main use | Required scan-size field | Input parameter columns |
|---|---|---|---|
| `RANDOM` | quick exploration of parameter ranges | `Number of points` = total trial points | `name, prior, min, max` |
| `GRID` | low-dimensional regular grids | intervals per input parameter | `name, prior, min, max, intervals` |
| `BESTFIT` | best-fit search with SciPy differential evolution | `Number of points` = max iterations | `name, prior, min, max` |
| `MCMC` | single-chain likelihood/posterior sampling | `Number of points` = target accepted/surviving points | `name, prior, min, max, interval, initial` |
| `EMCEE` | ensemble MCMC sampling | `Number of points` = target saved samples | `name, prior, min, max, interval, initial` |
| `DYNESTY` | nested sampling with Dynesty | `Number of points` = live points | `name, prior, min, max` |
| `MULTINEST` | nested sampling with MultiNest | `Number of points` = live points | `name, prior, min, max` |
| `ONEPOINT` | evaluate one fixed point | not used | `name, fixed, value` |
| `ONEPOINTBATCH` | evaluate a file of fixed points | not used | point file supplies values |

## Selection Rules

- Use `RANDOM` when the user gives parameter ranges and wants a first scan.
- Use `GRID` only when dimensionality is small and the user wants regular grid points; convert requested point counts to interval counts with `BinNum = N - 1`.
- Use `BESTFIT` when the goal is a best-fit point or chi-square minimization.
- Use `MCMC` when the user wants posterior-like sampling and can provide an initial point and interval divisor per parameter.
- Use `EMCEE` when the user asks for ensemble MCMC, walkers, or a more modern MCMC sampler.
- Use `DYNESTY` when the user asks for nested sampling without requiring the MultiNest library.
- Use `MULTINEST` only when the user explicitly wants MultiNest or Bayesian evidence and can install `pymultinest` plus the native MultiNest library.
- Use `ONEPOINT` or `ONEPOINTBATCH` when all parameters are fixed and the task is evaluation, not sampling.

## Method Notes

### RANDOM

```ini
[scan]
Scan method:       RANDOM
Input parameters:  x, flat, 0, 3.14
                   y, flat, -3.14, 3.14
Number of points:  1000
```

`Number of points` is the total number of trial points. `Random seed` is useful
for reproducible tests.

### GRID

```ini
[scan]
Scan method:       GRID
Input parameters:  x, flat, 0, 3.14, 20
                   y, flat, -3.14, 3.14, 20
```

The last input-parameter column is EasyScan_HEP's `BinNum`, meaning the number
of intervals, not the final number of grid points. EasyScan_HEP evaluates both
endpoints, so one parameter produces `BinNum + 1` points. If the user asks for
`N` grid points in a range, write `N-1` in the `.ini`. For example, "10 points
from 0 to pi" should be `0, 3.141592653589793, 9`. Total grid points are the
product of `(BinNum + 1)` over all scanned parameters.

### BESTFIT

```ini
[scan]
Scan method:       BESTFIT
Input parameters:  x, flat, 0, 3.14
                   y, flat, -3.14, 3.14
Number of points:  100
```

`BESTFIT` uses `scipy.optimize.differential_evolution`. It does not need
derivatives. `Number of points` is passed as the maximum optimizer iterations.
It is appropriate only when a likelihood or chi-square constraint is defined.

### MCMC

```ini
[scan]
Scan method:       MCMC
Input parameters:  x, flat, 0, 3.14, 10, 1.5
                   y, flat, -3.14, 3.14, 20, 0.0
Number of points:  10000
Acceptance rate:   0.25
```

The fifth column is an interval divisor. EasyScan uses `step = (max - min) /
interval`, so convert a user-provided proposal step `s` to `(max - min) / s`.
The sixth column is the initial value. Use `Acceptance rate` only for `MCMC`; a
typical value is `0.25`.

### EMCEE

```ini
[scan]
Scan method:       EMCEE
Input parameters:  x, flat, 0, 3.14, 10, 1.5
                   y, flat, -3.14, 3.14, 20, 0.0
Number of points:  10000
```

`EMCEE` requires the optional `emcee` package. It uses the same parameter-column
format as `MCMC`: the fifth column is an interval divisor and the sixth is the
initial value. EasyScan_HEP derives the number of walkers from the number of
dimensions and parallel threads, but the current config checker does not list
EMCEE as a parallel method; do not add parallel settings unless the user asks
or the checker has been updated.

### DYNESTY

```ini
[scan]
Scan method:       DYNESTY
Input parameters:  x, flat, 0, 3.14
                   y, flat, -3.14, 3.14
Number of points:  500
```

`DYNESTY` requires the optional `dynesty` package. `Number of points` is the
number of live points.

### MULTINEST

```ini
[scan]
Scan method:       MULTINEST
Input parameters:  x, flat, 0, 3.14
                   y, flat, -3.14, 3.14
Number of points:  500
Random seed:       1234
```

`MULTINEST` requires `pymultinest` and the native MultiNest library. Prefer
`DYNESTY` if the user wants a pip-friendlier nested sampler.

### ONEPOINT

```ini
[scan]
Scan method:       ONEPOINT
Input parameters:  x, fixed, 1.5
                   y, fixed, 0.0
```

Only fixed values are meaningful. `Number of points`, priors, bins, proposal
intervals, and random seeds are not used.

### ONEPOINTBATCH

```ini
[scan]
Scan method:       ONEPOINTBATCH/path/to/points.dat
Parallel threads:  4
```

The path after `ONEPOINTBATCH/` points to a tabular file of fixed parameter
values. Do not invent the batch file content unless the user supplies the
parameter names and rows.

The source templates also use `ONEPOINT/path/to/points.dat`; both forms are
parsed as `ONEPOINTBATCH` when a slash and file path are present.

## Parallel Notes

- Recommend parallel settings only for methods currently recognized by the
  checker: `RANDOM`, `GRID`, `MCMC`, `ONEPOINTBATCH`, and `MULTINEST`.
- When `Parallel threads > 1`, `Parallel folder` is required, must already
  exist, and will be copied to `p0_<folder>`, `p1_<folder>`, etc. It must be a
  prefix of each `Command path`; a top-level relative folder such as `utils` is
  safer than a nested path.
- Do not recommend parallel execution for `BESTFIT`, `DYNESTY`, or `EMCEE`
  unless the user explicitly accepts checker warnings or the code/checker is
  updated.

## Dependency Reminders

- Core EasyScan_HEP needs `numpy`, `scipy`, `matplotlib`, and `pandas`.
- `BESTFIT` uses SciPy and does not calculate derivatives.
- `DYNESTY` needs `python3 -m pip install dynesty` or `easyscan-hep[dynesty]`.
- `EMCEE` needs `python3 -m pip install emcee` or `easyscan-hep[emcee]`.
- `MULTINEST` needs `pymultinest` plus the native MultiNest library.
