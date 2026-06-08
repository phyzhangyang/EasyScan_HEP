# SSM_DM_EWPT EasyScan Example

This example scans the Z2 real-singlet scalar dark-matter model with EasyScan.
EasyScan writes the external-program input files, runs micrOMEGAs and PhaseTracer
directly, and reads their output files.
The default GRID setup uses 10 intervals in each direction, i.e. 11 x 11 scan
points in the `mS` and `lambdaHS` plane. Both directions use EasyScan's `flat`
prior, so the scan points are uniformly spaced in the plotted parameter plane.

Fixed external versions:

- micrOMEGAs 7.1
- PhaseTracer 2.2.1

Prepare the fixed-version external programs:

```bash
cd EasyScan_HEP
python3 utils/SSM_DM_EWPT/bootstrap.py
```

If the micrOMEGAs Zenodo download is unavailable, provide a local archive:

```bash
python3 utils/SSM_DM_EWPT/bootstrap.py --micromegas-tar /path/to/micromegas_7.1.tgz
```

PhaseTracer also needs its normal CMake dependencies, including Boost, Eigen,
ALGLIB and NLopt. Install those dependencies before running `bootstrap.py`; if
one is missing, the script stops during the PhaseTracer CMake step and asks the
user to install it. On macOS with Homebrew, install the commonly available
packages with:

```bash
brew install boost eigen nlopt pkgconf
```

ALGLIB must also be installed in a location visible to CMake; follow the
PhaseTracer/ALGLIB installation instructions for your system.

Then check and run the example:

```bash
python3 bin/easyscan.py --check templates/scan_SSM_DM_EWPT.ini --json
python3 bin/easyscan.py templates/scan_SSM_DM_EWPT.ini
```

The scanned `lambdaHS` uses the paper/micrOMEGAs convention. The PhaseTracer
input file receives `lambda_hs = 2 * lambdaHS`.
The micrOMEGAs patch only adds two EasyScan-readable labels to the native
`SingletDM/main` output: `ES_Omega_h2` and `ES_sigmaSIp_pb`.

EasyScan automatically writes color plots under `SSM_DM_EWPT_grid/Figures/`.
To make a compact dark-matter/EWPT summary plot from the same `ScanResult.txt`,
run:

```bash
python3 utils/SSM_DM_EWPT/plot_dm_ewpt.py
```
