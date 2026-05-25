# Common External Program Configuration Templates

This document provides EasyScan_HEP configuration templates for common high-energy physics programs.

---

## Critical Syntax Rules From EasyScan_HEP Runtime

- SLHA mappings must include `BLOCK` or `DECAY`: use
  `var, file_id, SLHA, BLOCK, BLOCKNAME, key...`, not
  `var, file_id, SLHA, BLOCKNAME, key`.
- Prefer relative paths from the directory where `easyscan` is launched, or `~`
  paths. Avoid macOS-style `/Users/...` absolute paths in generated configs
  unless the user explicitly asks.
- For Python external programs, write the intended interpreter explicitly in
  `Execute command`, for example `/usr/bin/python3 script.py` or
  `.venv/bin/python script.py`.
- `Program name` and `Command path` are optional in the parser, but include them
  in generated configs for readability and reproducibility.

---

## Program Installation Detection

Run the following commands to detect installed programs:

```bash
# Check MSSM_mW directory
ls -la utils/MSSM_mW/

# Should see the following subdirectories (if compiled):
# susyhit/      - SUSY-HIT
# feynhiggs/    - FeynHiggs
# higgsbounds/  - HiggsBounds
# gm2calc/      - GM2Calc
```

**Criteria**:
- Directory exists and contains executables → Installed and ready
- Directory exists but only has Makefile/source → Needs compilation
- Directory does not exist → Not installed

### Installation Commands (MSSM Programs)

```bash
cd utils/MSSM_mW
make susyhit     # SUSY-HIT 1.6
make feynhiggs   # FeynHiggs 2.19.0
make higgsbounds # HiggsBounds 5.10.2
make gm2calc     # GM2Calc 2.3.1
```

---

## 1. SUSY-HIT

### Program Info
- **Version**: 1.6
- **Purpose**: Calculate MSSM particle spectrum and decays (SuSpect + HDECAY + SDECAY)
- **Input Format**: SLHA (SUSY Les Houches Accord)
- **Output Format**: SLHA

### Configuration Template

```ini
[program1]
Program name:    SUSY-HIT
Execute command: ./run
Command path:    utils/MSSM_mW/susyhit/
Input file:      1, utils/MSSM_mW/for_MSSM_mW.susyhit

# SLHA input variable mapping
Input variable:  mu,        1,      SLHA,     BLOCK,  EXTPAR,  31
                 tanb,      1,      SLHA,     BLOCK,  MINPAR,  5
                 mA,        1,      SLHA,     BLOCK,  EXTPAR,  24
                 MQ3,       1,      SLHA,     BLOCK,  EXTPAR,  7
                 MU3,       1,      SLHA,     BLOCK,  EXTPAR,  8
                 MD3,       1,      SLHA,     BLOCK,  EXTPAR,  9
                 At,        1,      SLHA,     BLOCK,  EXTPAR,  10
                 Ab,        1,      SLHA,     BLOCK,  EXTPAR,  11

# SLHA output variable mapping
Output file:     1, utils/MSSM_mW/susyhit/output.susyhit
Output variable: mh1,       1,      SLHA,     BLOCK,  HMIX,    2
                 mA_out,    1,      SLHA,     BLOCK,  HMIX,    3
                 mH2,       1,      SLHA,     BLOCK,  HMIX,    4
                 mHpm,      1,      SLHA,     BLOCK,  HMIX,    5
```

### SLHA Block Description

| Block | Description | Common Entries |
|-------|-------------|----------------|
| `MINPAR` | Minimal parameters | 3 (m0), 4 (m1/2), 5 (tanb) |
| `EXTPAR` | Extended parameters | 24 (mA), 31 (mu), 34-39 (soft-breaking masses) |
| `HMIX` | Higgs mixing | 2 (mh1), 3 (mA), 4 (mH2) |

### Notes
- Input file must be in complete SLHA format
- Output file contains complete particle spectrum and decay information
- Runtime: ~1-10 seconds/point

---

## 2. HiggsBounds

### Program Info
- **Version**: 5.10.2
- **Purpose**: Test Higgs spectrum against LEP/Tevatron/LHC exclusion limits
- **Input Format**: SLHA (from SUSY-HIT)
- **Output Format**: Text (includes observed values, expected values, etc.)

### Configuration Template

```ini
[program2]
Program name:    HiggsBounds
Execute command: ./HiggsBounds
Command path:    utils/MSSM_mW/higgsbounds/
Input file:      1, utils/MSSM_mW/susyhit/output.susyhit

# Read Higgs exclusion information from output file
Output file:     1, utils/MSSM_mW/higgsbounds/HiggsBoundsResults.dat
Output variable: obs,       1,      Label,    "obs",  2
                 exp,        1,      Label,    "exp",  2
```

### Output Description

Typical `HiggsBoundsResults.dat` format:
```
obs  0.523  # observed signal strength
exp  1.000  # expected exclusion limit
```

### Constraint Application

```ini
[constraint]
# Higgs exclusion constraint (free-form chi-square)
FreeFormChi2:   obs
```

### Notes
- Depends on FeynHiggs (optional, for more precise Higgs mass calculation)
- Runtime: ~1 second/point

---

## 3. GM2Calc

### Program Info
- **Version**: 2.3.1
- **Purpose**: Calculate anomalous magnetic moment (g-2), weak mixing angle, W mass
- **Input Format**: SLHA
- **Output Format**: Text/SLHA

### Configuration Template

```ini
[program3]
Program name:    GM2Calc
Execute command: ./gm2calc
Command path:    utils/MSSM_mW/gm2calc/
Input file:      1, utils/MSSM_mW/susyhit/output.susyhit

# Output variables (read by position)
Output file:     1, utils/MSSM_mW/gm2calc/gm2_output.dat
Output variable: gm2,      1,      Position, 5,  2
                 mW,       1,      Position, 8,  2
                 sw2eff,   1,      Position, 10, 2
```

### Typical Output Format

```
# GM2Calc output
# a_muon (g-2): 2.51e-9
# mW: 80.452 GeV
# sin^2(theta_W^eff): 0.23121
```

### Constraint Application

```ini
[constraint]
# W mass constraint
Gaussian:   mW,      80.4,    0.01

# Weak mixing angle constraint
Gaussian:   sw2eff,  0.23121, 0.000108

# muon g-2 constraint
Gaussian:   gm2,     25.1e-10, 6.2e-10
```

### Notes
- Depends on SUSY-HIT output for complete SLHA file
- Runtime: <1 second/point

---

## 4. FeynHiggs

### Program Info
- **Version**: 2.19.0
- **Purpose**: Precise calculation of MSSM Higgs masses and properties
- **Input Format**: SLHA
- **Output Format**: Text/SLHA

### Configuration Template

```ini
[programX]
Program name:    FeynHiggs
Execute command: ./FeynHiggs
Command path:    utils/MSSM_mW/feynhiggs/
Input file:      1, utils/MSSM_mW/susyhit/output.susyhit

Output file:     1, utils/MSSM_mW/feynhiggs/fh.out
Output variable: mh1_fh,   1,      Label,    "mh0",  2
                 mH2_fh,   1,      Label,    "mH0",  2
```

### Coordination with HiggsBounds

FeynHiggs can be used together with HiggsBounds to provide more precise Higgs mass input.

### SLHA Block Reference

Common SLHA blocks for output variable extraction:

| Block | Entry | Description | Example Value |
|-------|-------|-------------|---------------|
| **HMIX** | 1 | h1 mass | 125.1 |
| | 2 | h2 mass | 125.1 |
| | 3 | A mass | 300.0 |
| | 4 | H2 mass | 300.0 |
| | 5 | H+ mass | 310.0 |
| **MASS** | 6 | Top quark mass | 173.0 |
| | 24 | W boson mass | 80.4 |
| | 25 | SM-like Higgs mass | 125.1 |
| **PRECOBS** | 1 | W mass prediction | 80.45 |
| | 2 | sin^2(theta_W) | 0.231 |
| **GM2CalcOutput** | 5 | a_muon (g-2) | 2.51e-9 |
| | 8 | mW prediction | 80.452 |
| | 10 | sin^2(theta_W^eff) | 0.23121 |

---

## 5. Generic Program Templates

### Position Method (Simple Text)

For programs with simple format output.

```ini
[programX]
Program name:    MyProgram
Execute command: ./run
Command path:    utils/myprogram/
Input file:      1, utils/myprogram/input.dat

# Write input by position
Input variable:  param1,  1,      Position, 1,  1
                 param2,  1,      Position, 2,  1

Output file:     1, utils/myprogram/output.dat
# Read output by position
Output variable: result1, 1,      Position, 1,  2
                 result2, 1,      Position, 2,  2
```

### Label Method (Labeled Text)

For programs with clear label markers.

```ini
[programX]
Program name:    MyProgram
Execute command: ./run
Command path:    utils/myprogram/
Input file:      1, utils/myprogram/input.dat

# Write input by label
Input variable:  mass,    1,      Label,    "MASS=",  2
                 width,   1,      Label,    "WIDTH=", 2

Output file:     1, utils/myprogram/output.dat
# Read output by label
Output variable: xsec,    1,      Label,    "XSEC=",  2
```

### Replace Method (Parameter Replacement)

For scenarios requiring replacement of specific parameters in files.

```ini
[programX]
Program name:    MyProgram
Execute command: ./run
Command path:    utils/myprogram/
Input file:      1, utils/myprogram/card.dat

# Replace parameters in file
Input variable:  tanb,    1,      Replace,  "tanb"
                 mA,      1,      Replace,  "mA"
                 mu,      1,      Replace,  "mu"
```

### SLHA Method (HEP Standard Format)

For SUSY model calculations.

```ini
[programX]
Program name:    SUSYProgram
Execute command: ./run
Command path:    utils/susyprog/
Input file:      1, utils/susyprog/input.slha

# SLHA BLOCK input
Input variable:  mA,      1,      SLHA,     BLOCK,  MODSEL,  3
                 mu,      1,      SLHA,     BLOCK,  EXTPAR,  31

# SLHA BLOCK output
Output file:     1, utils/susyprog/output.slha
Output variable: mh1,     1,      SLHA,     BLOCK,  HMIX,    2
```

---

## Multi-Program Chain Configuration

Typical MSSM scan involves chained execution of multiple programs:

```ini
# Program 1: SUSY-HIT calculates particle spectrum
[program1]
Program name:    SUSY-HIT
Execute command: ./run
Command path:    utils/MSSM_mW/susyhit/
Input file:      1, utils/MSSM_mW/for_MSSM_mW.susyhit
Input variable:  mu,      1,      SLHA,     BLOCK,  EXTPAR,  31
Output file:     1, utils/MSSM_mW/susyhit/output.susyhit
Output variable: mh1,     1,      SLHA,     BLOCK,  HMIX,    2

# Program 2: HiggsBounds tests exclusion
[program2]
Program name:    HiggsBounds
Execute command: ./HiggsBounds
Command path:    utils/MSSM_mW/higgsbounds/
Input file:      1, utils/MSSM_mW/susyhit/output.susyhit
Output file:     1, utils/MSSM_mW/higgsbounds/HiggsBoundsResults.dat
Output variable: obs,     1,      Label,    "obs",  2

# Program 3: GM2Calc calculates g-2 and W mass
[program3]
Program name:    GM2Calc
Execute command: ./gm2calc
Command path:    utils/MSSM_mW/gm2calc/
Input file:      1, utils/MSSM_mW/susyhit/output.susyhit
Output file:     1, utils/MSSM_mW/gm2calc/gm2_output.dat
Output variable: gm2,     1,      Position, 5,  2
                 mW,      1,      Position, 8,  2
```

**Execution Order**: program1 → program2 → program3
**Parameter Passing**: Scan parameters → SUSY-HIT → (output SLHA) → HiggsBounds + GM2Calc

---

## Time Limit Setting

For programs that may hang, set a timeout:

```ini
[programX]
Program name:    SlowProgram
Execute command: ./run
Command path:    utils/slowprog/
Time limit in minute: 10  # 10 minute timeout
```

Note: After setting time limit, program will use `subprocess.popen` instead of `os.system`.

---

## Boundary Constraint Setting

For scenarios requiring boundary reading from files:

```ini
[programX]
Program name:    MyProgram
Execute command: ./run
Command path:    utils/myprogram/
Input file:      1, utils/myprogram/input.dat
Input variable:  x,       1,      Position, 1,  1
                 y,       1,      Position, 1,  2
Output file:     1, utils/myprogram/output.dat
Output variable: f,       1,      Position, 1,  2

# Boundary constraint: y <= f(x), boundary read from file
Bound:           f,       x,      MIN,     templates/bound.txt
```

**bound.txt format**:
```
# x  y_boundary
0    1.0
1    2.0
2    3.0
...
```

---

## Schematic MSSM Scan Example

This is a compact chained-program example. Treat the exact parameter list as a
template and align names/SLHA entries with the real input card before generating
a user config.

```ini
[scan]
Result folder name:  MSSM_mW_scan
Scan method:       MCMC
Input parameters:  mu,      flat,   100,   1000,   50,   500
                   tanb,    flat,   2,     60,     2,    10
                   mA,      flat,   100,   2000,   100,  500
                   MQ3,     flat,   100,   5000,   200,  2000
                   MD3,     flat,   100,   5000,   200,  2000
                   ML3,     flat,   100,   5000,   200,  2000
                   MU3,     flat,   100,   5000,   200,  2000
                   ME3,     flat,   100,   5000,   200,  2000
                   At,      flat,   -5000, 5000,   500,  0
                   Ab,      flat,   -5000, 5000,   500,  0
                   Atau,    flat,   -5000, 5000,   500,  0
                   M1,      flat,   100,   2000,   100,  500
                   M2,      flat,   100,   2000,   100,  500
                   M3,      flat,   100,   5000,   200,  2000
Number of points:  10000
Parallel threads:  4
Parallel folder:   utils

[program1]
Program name:    SUSY-HIT
Execute command: ./run
Command path:    utils/MSSM_mW/susyhit/
Input file:      1, utils/MSSM_mW/for_MSSM_mW.susyhit
Input variable:  mu,      1,      SLHA,     BLOCK,  EXTPAR,  31
                 tanb,    1,      SLHA,     BLOCK,  MINPAR,  5
                 mA,      1,      SLHA,     BLOCK,  EXTPAR,  24
                 MQ3,     1,      SLHA,     BLOCK,  EXTPAR,  7
                 MU3,     1,      SLHA,     BLOCK,  EXTPAR,  8
                 MD3,     1,      SLHA,     BLOCK,  EXTPAR,  9
                 ML3,     1,      SLHA,     BLOCK,  EXTPAR,  12
                 ME3,     1,      SLHA,     BLOCK,  EXTPAR,  14
                 At,      1,      SLHA,     BLOCK,  EXTPAR,  10
                 Ab,      1,      SLHA,     BLOCK,  EXTPAR,  11
                 Atau,    1,      SLHA,     BLOCK,  EXTPAR,  13
                 M1,      1,      SLHA,     BLOCK,  EXTPAR,  41
                 M2,      1,      SLHA,     BLOCK,  EXTPAR,  42
                 M3,      1,      SLHA,     BLOCK,  EXTPAR,  43
Output file:     1, utils/MSSM_mW/susyhit/output.susyhit
Output variable: mh1,     1,      SLHA,     BLOCK,  HMIX,    2
                 mA_out,  1,      SLHA,     BLOCK,  HMIX,    3

[program2]
Program name:    HiggsBounds
Execute command: ./HiggsBounds
Command path:    utils/MSSM_mW/higgsbounds/
Input file:      1, utils/MSSM_mW/susyhit/output.susyhit
Output file:     1, utils/MSSM_mW/higgsbounds/HiggsBoundsResults.dat
Output variable: obs,     1,      Label,    "obs",  2

[program3]
Program name:    GM2Calc
Execute command: ./gm2calc
Command path:    utils/MSSM_mW/gm2calc/
Input file:      1, utils/MSSM_mW/susyhit/output.susyhit
Output file:     1, utils/MSSM_mW/gm2calc/gm2_output.dat
Output variable: gm2,     1,      Position, 5,  2
                 mW,      1,      Position, 8,  2

[constraint]
# W mass constraint
Gaussian:   mW,      80.4,    0.01

# Muon g-2 constraint
Gaussian:   gm2,     25.1e-10, 6.2e-10

# Higgs mass constraint
Gaussian:   mh1,     125.4,   3.0

# Higgs exclusion (free-form chi-square from HiggsBounds)
FreeFormChi2:   obs

[plot]
Color: mu,    tanb,      mW,      mW_mu_tanb
       mu,    mA,        Chi2
```

**Execution Flow**:
1. MCMC generates 17-dimensional parameter point
2. SUSY-HIT calculates particle spectrum (writes SLHA output)
3. HiggsBounds reads SLHA, tests against LEP/Tevatron/LHC limits
4. GM2Calc reads SLHA, calculates mW and g-2
5. All constraints combined into total chi-square
6. MCMC accepts/rejects point based on likelihood

**Runtime**: ~4 hours with 4 parallel threads (10000 points)
