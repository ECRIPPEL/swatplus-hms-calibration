# swatplus-hms-calibration

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19793607.svg)](https://doi.org/10.5281/zenodo.19793607)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Hierarchical Multi-Scale Calibration of SWAT+ guided by Adaptive Sensitivity Analysis.**

A reproducible R-based framework for calibrating SWAT+ hydrological models across three temporal scales (annual water balance → monthly streamflow → daily streamflow), using Morris sensitivity screening re-executed at each phase, σ-informed parameter-range narrowing, and cross-phase robustness filtering to manage equifinality. Built on top of [SWATrunR](https://chrisschuerz.github.io/SWATrunR/).

If you use this software in research, please see [Citation](#citation) below.

---

## Table of Contents

1. [Prerequisites](#1-prerequisites)
2. [Project Structure](#2-project-structure)
3. [Before You Start](#3-before-you-start)
4. [Understanding the Workflow](#4-understanding-the-workflow)
5. [Configuration Files (YAML)](#5-configuration-files-yaml)
6. [Running the Full Pipeline](#6-running-the-full-pipeline)
7. [Running Individual Steps](#7-running-individual-steps)
8. [Resuming After Interruption](#8-resuming-after-interruption)
9. [Understanding the Outputs](#9-understanding-the-outputs)
10. [Adapting to Your Catchment](#10-adapting-to-your-catchment)
11. [Troubleshooting](#11-troubleshooting)

---

## 1. Prerequisites

### Software

| Software   | Version  | Notes                                                                        |
|------------|----------|------------------------------------------------------------------------------|
| R          | >= 4.3.3                                 | Required for SWATrunR compatibility          |
| SWAT+      | >= swatplus-61.0.2.61-ifx-win_amd64-Rel  | Model executable must be inside project_path |
| RTools     | >= 4.3                                   | Windows only — needed to compile R packages  |

### R Packages

Install all required packages before running:

```r
install.packages(c(
  "SWATrunR",     # SWAT+ simulation runner
  "yaml",         # YAML config parser
  "dplyr",        # Data manipulation
  "tibble",       # Modern data frames
  "tidyr",        # Data reshaping
  "purrr",        # Functional programming
  "hydroGOF",     # NSE, KGE, PBIAS functions
  "lubridate",    # Date handling
  "ggplot2",      # Plotting
  "ggrepel",      # Label placement in plots
  "lhs"           #Latin Hypercube Sampling
))
```

> **Note:** If `SWATrunR` is not available on CRAN, install from GitHub:
> ```r
> remotes::install_github("chrisschuerz/SWATrunR")
> ```

### SWAT+ Project

Your SWAT+ TxtInOut folder (referred to as `project_path` in configs) must contain:
- A working SWAT+ model that runs successfully
- `cal_parms.cal` (used for parameter range validation in Phase 1)
- `soils.sol`, `hydrology.hyd`, `topography.hyd` (read during range checks)

### Observed Data

You need CSV files with observed streamflow placed inside your `project_path`:

| File                   | Columns       | Format              | Used in        |
|------------------------|---------------|---------------------|----------------|
| `obs_flow_monthly.csv` | Date, Flow    | YYYY-MM-DD, m^3/s   | Phases 1, 2    |
| `obs_flow_daily.csv`   | Date, Flow    | YYYY-MM-DD, m^3/s   | Phase 3        |

- `Date` must be in `YYYY-MM-DD` format (e.g., `1980-01-01`)
- `Flow` is streamflow in m^3/s
- Monthly file: one row per month (use the 1st day of each month as date)
- Daily file: one row per day

---

## 2. Project Structure

```
swatplus-hms-calibration/
├── config/
│   ├── phase1.yaml         # Phase 1: annual water balance (soft calibration)
│   ├── phase2.yaml         # Phase 2: monthly streamflow (hard calibration)
│   └── phase3.yaml         # Phase 3: daily streamflow (hard calibration)
├── R/
│   ├── check_range.R       # Parameter range validation vs cal_parms.cal
│   ├── lhoat_engine.R      # LH-OAT trajectory generation + Elementary Effects
│   ├── metrics.R           # NSE, KGE, PBIAS, water balance metrics
│   ├── morris_classify.R   # Morris screening classification + range narrowing
│   ├── param_inherit.R     # Cross-phase parameter inheritance
│   ├── plots.R             # Publication-ready scientific plots (TIFF)
│   ├── run_filter.R        # Invalid run detection and filtering
│   └── versioning.R        # Auto-versioned output folder management
├── master.R                # Orchestrator: runs all 6 steps sequentially
├── run_gsa.R               # Generic GSA script (works for any phase)
└── run_cal.R               # Generic calibration script (works for any phase)
```

**You should NOT need to edit any `.R` file.** All customization is done through the YAML config files.

---

## 3. Before You Start

### Set Your Working Directory

All scripts expect to run from the project root (`swatplus-hms-calibration/`).

**In RStudio:**
```r
setwd("C:/path/to/your/swatplus/txtinout/swatplus-hms-calibration")
```

**In VSCode (R terminal):**
```r
setwd("C:/path/to/your/swatplus/txtinout/swatplus-hms-calibration")
```

**From command line:**
```bash
cd C:/path/to/your/swatplus/txtinout/swatplus-hms-calibration
```

### Verify Your SWAT+ Model Runs

Before calibration, ensure your SWAT+ model runs successfully:

```r
library(SWATrunR)
sim_test <- run_swatplus(
  project_path = "C:/path/to/your/swatplus/txtinout",
  output = define_output(file = "channel_sd", variable = "flo_out", unit = 33),
  start_date = "1979-01-01",
  end_date   = "1980-12-31",
  years_skip = 1,
  n_thread   = 4
)
```

If this fails, fix your SWAT+ model before proceeding.

---

## 4. Understanding the Workflow

The calibration follows a **hierarchical multi-scale** approach with 3 phases, each consisting of a **Global Sensitivity Analysis (GSA)** followed by a **Calibration (CAL)**:

```
GSA Phase 1  →  CAL Phase 1  →  GSA Phase 2  →  CAL Phase 2  →  GSA Phase 3  →  CAL Phase 3
 (step 1)        (step 2)        (step 3)        (step 4)        (step 5)        (step 6)
```

### Phase 1 — Annual Water Balance (Soft Calibration)
- **Scale:** Annual
- **Goal:** Filter parameter sets that reproduce realistic water balance (ET/P, WYLD/P ratios)
- **Method:** Behavioural filtering — runs within ±15% of target ratios are kept
- **Parameters:** 12 hydrological parameters (cn2, awc, esco, epco, perco, etc.)
- **GSA:** LH-OAT Morris screening on annual water balance ratios
- **Output:** Narrowed parameter ranges passed to Phase 2

### Phase 2 — Monthly Streamflow (Hard Calibration)
- **Scale:** Monthly
- **Goal:** Calibrate to monthly observed streamflow
- **Method:** Hard thresholds — NSE >= 0.70, KGE >= 0.70, |PBIAS| <= 10%
- **Parameters:** 13 parameters (inherited from Phase 1 + refined ranges)
- **GSA:** LH-OAT Morris screening on monthly NSE/KGE/PBIAS
- **Output:** Further narrowed ranges passed to Phase 3

### Phase 3 — Daily Streamflow (Hard Calibration)
- **Scale:** Daily
- **Goal:** Calibrate to daily observed streamflow
- **Method:** Hard thresholds — NSE >= 0.50, KGE >= 0.50, |PBIAS| <= 15%
- **Parameters:** 17 parameters (inherited from Phase 2 + 4 new: surlag, lattime, mann, ovn)
- **GSA:** LH-OAT Morris screening on daily NSE/KGE/PBIAS
- **Output:** Final calibrated parameter sets

### How Phases Connect

Each phase inherits parameter ranges from the previous one:

1. **Phase 1 CAL** produces `novos_ranges_phase2_from_morris_behavioral.csv`
2. **Phase 2 CAL** reads those ranges and produces `novos_ranges_phase3_from_hard_phase2.csv`
3. **Phase 3 CAL** reads those ranges for the final calibration

Parameter ranges are **narrowed progressively** using quantile rules based on Morris sensitivity groups:
- **High importance + interaction** (mu* > mean, sigma > mean): Q5–Q95 (wide exploration)
- **High importance** (mu* > mean, sigma <= mean): Q10–Q90 (moderate narrowing)
- **Low importance** (mu* <= mean): Q25–Q75 (aggressive narrowing)

---

## 5. Configuration Files (YAML)

Each YAML config controls one phase. Here are the key sections you may want to customize:

### General Settings

```yaml
phase: 1                                              # Phase number (1, 2, or 3)
project_path: "C:/path/to/your/swatplus/txtinout"     # Path to your SWAT+ TxtInOut folder
seed: 123                                             # Random seed for reproducibility
threads: 16                                           # Number of parallel SWAT+ threads
script_version: "1.0.0"                               # Version tag saved with results
```

### GSA Settings

```yaml
gsa:
  prefix: "GSA_phase1_v"            # Output folder prefix (auto-versioned)
  m_traj: 1                         # Number of LH-OAT trajectories
  delta: 0.20                       # Perturbation step size (fraction of range)
  start_date: "1979-01-01"          # GSA simulation period
  end_date: "1980-12-31"
  warmup: 1                         # Warmup years (skipped in analysis)
```

> **Tip:** For Phase 1 screening, `m_traj: 1` is often enough. For Phases 2-3, use
> `m_traj: 15` or more for robust Morris indices. More trajectories = more SWAT+ runs
> but more reliable sensitivity ranking. Total runs = m_traj × (p + 1), where p is
> the number of parameters.

### Calibration Settings

**Soft calibration (Phase 1):**
```yaml
cal:
  type: "soft"                      # Behavioural filtering
  prefix: "SOFT_phase1_v"           # Output folder prefix
  n_simulations: 250                # Number of random parameter sets
  threshold_pct: 15                 # Maximum relative error (%) for water balance
  targets:
    et_rto: 0.48147                 # Target ET/P ratio
    wyld_rto: 0.51853               # Target WYLD/P ratio (= 1 - ET/P)
```

**Hard calibration (Phases 2-3):**
```yaml
cal:
  type: "hard"                      # Performance thresholds
  prefix: "HARD_phase2_v"
  n_simulations: 250
  threshold_nse: 0.70               # Minimum NSE
  threshold_kge: 0.70               # Minimum KGE
  threshold_abs_pbias: 10           # Maximum |PBIAS| (%)
```

### Parameter Definitions (Phase 1 only)

Phase 1 defines parameters explicitly. Each parameter needs:

```yaml
parameters:
  - short: "cn2"                                        # Short name (used in plots)
    par_name: "cn2::cn2.hru | change = pctchg"          # SWATrunR notation
    min: -15                                             # Lower bound
    max: 15                                              # Upper bound
```

**SWATrunR parameter notation:**
- `variable::file | change = type`
- Change types: `absval` (absolute value), `pctchg` (% change), `abschg` (absolute change)
- Example: `awc::soils.sol | change = pctchg` means modify AWC in soils.sol by a percentage

### Parameter Inheritance (Phases 2-3)

Phases 2-3 inherit ranges from the previous phase:

```yaml
inherit:
  from_prefix: "SOFT_phase1_v"                                          # Folder prefix to search
  ranges_file: "novos_ranges_phase2_from_morris_behavioral.csv"         # CSV with new ranges

param_map:
  cn2:    "cn2::cn2.hru | change = pctchg"
  awc:    "awc::soils.sol | change = pctchg"
  # ... maps short names to SWATrunR notation
```

### Adding Extra Parameters in Later Phases

Phase 3 introduces new parameters not present in Phase 2:

```yaml
extra_params:
  - short: "surlag"
    par_name: "surlag::bsn_prm.bsn | change = absval"
    min: 0.5
    max: 24
```

---

## 6. Running the Full Pipeline

### Option A: Run Everything from Step 1

```r
setwd("C:/path/to/your/swatplus/txtinout/swatplus-hms-calibration")
source("master.R")
```

This will execute all 6 steps sequentially:
```
Step 1/6 - GSA Phase 1
Step 2/6 - CAL Phase 1
Step 3/6 - GSA Phase 2
Step 4/6 - CAL Phase 2
Step 5/6 - GSA Phase 3
Step 6/6 - CAL Phase 3
```

Each step prints timing information when completed.

### Option B: From the Command Line

```bash
cd C:/path/to/your/swatplus/txtinout/swatplus-hms-calibration
Rscript master.R
```

---

## 7. Running Individual Steps

You can run any step independently.

### GSA (any phase)

**From R console / RStudio / VSCode:**
```r
PHASE_CONFIG <- "config/phase1.yaml"
source("run_gsa.R")
```

**From command line:**
```bash
Rscript run_gsa.R config/phase1.yaml
Rscript run_gsa.R config/phase2.yaml
Rscript run_gsa.R config/phase3.yaml
```

### Calibration (any phase)

**From R console / RStudio / VSCode:**
```r
PHASE_CONFIG <- "config/phase2.yaml"
source("run_cal.R")
```

**From command line:**
```bash
Rscript run_cal.R config/phase1.yaml
Rscript run_cal.R config/phase2.yaml
Rscript run_cal.R config/phase3.yaml
```

> **Important:** Calibration steps depend on GSA results from the same phase.
> Always run GSA before CAL for each phase. Phase 2+ also depends on calibration
> results from the previous phase (for parameter range inheritance).

### Correct Execution Order

If running steps individually, follow this order:

```
1. Rscript run_gsa.R config/phase1.yaml    # GSA Phase 1
2. Rscript run_cal.R config/phase1.yaml    # CAL Phase 1 (needs GSA Phase 1 results)
3. Rscript run_gsa.R config/phase2.yaml    # GSA Phase 2 (needs CAL Phase 1 ranges)
4. Rscript run_cal.R config/phase2.yaml    # CAL Phase 2 (needs GSA Phase 2 results)
5. Rscript run_gsa.R config/phase3.yaml    # GSA Phase 3 (needs CAL Phase 2 ranges)
6. Rscript run_cal.R config/phase3.yaml    # CAL Phase 3 (needs GSA Phase 3 results)
```

---

## 8. Resuming After Interruption

If the pipeline is interrupted (crash, power loss, manual stop), you can resume from any step:

```r
START_FROM <- 3        # Resume from GSA Phase 2 (step 3)
source("master.R")
```

Steps before `START_FROM` will be skipped with `[SKIP]` in the console output.

### Step Reference

| Step | Label        | What It Does                                    |
|------|--------------|-------------------------------------------------|
| 1    | GSA Phase 1  | Sensitivity analysis on annual water balance     |
| 2    | CAL Phase 1  | Soft calibration (behavioural filtering)         |
| 3    | GSA Phase 2  | Sensitivity analysis on monthly streamflow       |
| 4    | CAL Phase 2  | Hard calibration (NSE/KGE/PBIAS thresholds)      |
| 5    | GSA Phase 3  | Sensitivity analysis on daily streamflow          |
| 6    | CAL Phase 3  | Hard calibration on daily streamflow              |

> **Tip:** If a step fails, fix the issue and resume from that step. You don't need
> to re-run previous steps — their results are already saved in versioned folders.

---

## 9. Understanding the Outputs

### Output Folder Structure

Each step creates an auto-versioned folder inside your `project_path`:

```
C:/path/to/your/swatplus/txtinout/
├── GSA_phase1_v1/          # GSA Phase 1 results (1st run)
├── GSA_phase1_v2/          # GSA Phase 1 results (2nd run, if re-run)
├── SOFT_phase1_v1/         # CAL Phase 1 results
├── GSA_phase2_v1/          # GSA Phase 2 results
├── HARD_phase2_v1/         # CAL Phase 2 results
├── GSA_phase3_v1/          # GSA Phase 3 results
└── HARD_phase3_v1/         # CAL Phase 3 results
```

Version numbers auto-increment. If you re-run a step, a new folder (v2, v3...) is created — previous results are never overwritten.

### GSA Output Files

| File                                    | Description                                        |
|-----------------------------------------|----------------------------------------------------|
| `resultado_lhoat_phaseX.rds`            | Full R object with all GSA results                 |
| `ranking_global_lhoat_phaseX.csv`       | Morris sensitivity ranking (mu*, sigma per param)  |
| `metricas_lhoat_phaseX.csv`             | Performance metrics for all LH-OAT runs            |
| `ee_all_lhoat_phaseX.csv`              | Elementary Effects for all trajectories             |
| `param_info_lhoat_phaseX.csv`          | Parameter info (names, ranges) used in GSA         |
| `Morris_sensitivity_screening_phaseX.tif` | Morris plot (mu* vs sigma) — publication-ready   |
| `check_param_range_phase1.csv`          | Parameter range validation (Phase 1 only)          |

### Calibration Output Files

| File                                            | Description                                      |
|-------------------------------------------------|--------------------------------------------------|
| `resultado_cal_phaseX.rds`                      | Full R object with all calibration results       |
| `results_cal_phaseX.csv`                        | All runs with metrics + parameter values         |
| `satisfactory_results_phaseX.csv`               | Only runs that passed the thresholds             |
| `novos_ranges_phaseY_from_*.csv`                | Narrowed parameter ranges for next phase         |
| `scatter_performance_phaseX.tif`                | NSE vs KGE scatter (hard calibration)            |
| `balance_scatter_phase1.tif`                    | ET/P vs WYLD/P scatter (soft calibration)        |
| `hydrograph_phaseX.tif`                         | Observed vs simulated streamflow                 |
| `boxplot_params_phaseX.tif`                     | Normalized parameter distributions (behavioural) |

### Key Columns in Results CSV

**Soft calibration (Phase 1):**
- `rto_et`, `rto_wyld` — simulated ET/P and WYLD/P ratios
- `err_rto_et`, `err_rto_wyld` — relative errors vs targets
- `satisfactory_run` — TRUE if both errors <= threshold

**Hard calibration (Phases 2-3):**
- `NSE_mon` or `NSE_day` — Nash-Sutcliffe Efficiency
- `KGE_mon` or `KGE_day` — Kling-Gupta Efficiency
- `PBIAS_mon` or `PBIAS_day` — Percent Bias
- `obj_total` — Combined objective function
- `satisfactory_run` — TRUE if NSE, KGE, and PBIAS pass thresholds

### Loading Saved Results in R

```r
# Load full results object
res <- readRDS("C:/path/to/your/swatplus/txtinout/HARD_phase2_v1/resultado_cal_phase2.rds")

# Access components
res$config              # YAML config used
res$param_info           # Parameters and ranges
res$results              # All runs with metrics
res$satisfactory_results # Behavioural runs only
res$sim                  # Raw SWATrunR simulation output

# Quick summary
cat("Total runs:", nrow(res$results), "\n")
cat("Satisfactory:", nrow(res$satisfactory_results), "\n")
```

---

## 10. Adapting to Your Catchment

### Step-by-Step Customization

#### 1. Set your project path

Edit all three YAML configs:
```yaml
project_path: "C:/path/to/your/swatplus/txtinout"
```

#### 2. Set your outlet channel unit

In the YAML config files, find the output definitions and update the `unit` to match your outlet channel number:
```yaml
gsa_outputs:
  - file: "channel_sd"
    variable: "flo_out"
    unit: 33                  # <-- Change to YOUR outlet channel number
```

> **How to find your outlet channel:** Open `channel-lte.cha` in your TxtInOut folder.
> The outlet is typically the channel with the largest drainage area or the one
> corresponding to your gauging station location.

#### 3. Update observed data file names

```yaml
obs_files:
  monthly: "obs_flow_monthly.csv"    # Your monthly observed streamflow file
  daily:   "obs_flow_daily.csv"      # Your daily observed streamflow file
```

#### 4. Set simulation period and warmup

```yaml
start_date: "1976-01-01"     # Simulation start (including warmup)
end_date:   "1989-12-31"     # Simulation end
warmup: 4                    # Years to skip (warmup period)
obs_start_date: "1980-01-01" # First date of observed data to use
```

#### 5. Set water balance targets (Phase 1)

Calculate from observed data or literature for your catchment:
```yaml
cal:
  targets:
    et_rto: 0.48147          # Long-term ET / Precipitation ratio
    wyld_rto: 0.51853        # Long-term Water Yield / Precipitation ratio (= 1 - ET/P)
```

> **Tip:** `ET/P + WYLD/P ≈ 1.0`. Use long-term annual averages from your catchment.
> ET can be estimated from remote sensing (MODIS, GLEAM) or water balance closure.

#### 6. Adjust parameters for your region

Phase 1 parameters are defined directly in `config/phase1.yaml`. You can:
- **Add parameters:** Add new entries to the `parameters` list
- **Remove parameters:** Delete entries you don't need
- **Change ranges:** Adjust `min` and `max` values

For Phases 2-3, the `param_map` section controls which parameters are carried forward.

#### 7. Adjust calibration thresholds

Based on your data quality and expectations:
```yaml
# Phase 2 (monthly) — stricter thresholds
cal:
  threshold_nse: 0.70
  threshold_kge: 0.70
  threshold_abs_pbias: 10

# Phase 3 (daily) — relaxed thresholds
cal:
  threshold_nse: 0.50
  threshold_kge: 0.50
  threshold_abs_pbias: 15
```

> **Reference:** Thresholds follow Moriasi et al. (2015) guidelines.
> NSE > 0.50 is "satisfactory", > 0.65 is "good", > 0.75 is "very good".

#### 8. Adjust computational resources

```yaml
threads: 16                  # Number of CPU cores for parallel SWAT+ runs
cal:
  n_simulations: 250         # More simulations = better exploration but slower
gsa:
  m_traj: 15                 # More trajectories = more robust Morris indices
```

---

## 11. Troubleshooting

### Common Errors

**"Error: Usage: Rscript run_gsa.R config/phaseX.yaml"**
- You forgot to set `PHASE_CONFIG` before `source()`, or didn't pass the YAML path as argument.

**"Error in yaml::read_yaml(...) : file not found"**
- Check your working directory: `getwd()`. It must be the `swatplus-hms-calibration/` folder.

**"Error: ABORTING: some parameters exceed the cal_parms.cal bounds"**
- Phase 1 range check found parameters that would produce values outside SWAT+ limits.
- Review the printed table and adjust `min`/`max` in `config/phase1.yaml`.

**"Error: No folder found with prefix '...'"**
- The inheritance chain is broken. For Phase 2, you need Phase 1 CAL results first.
- Run the previous phase before proceeding.

**"Error: Column 'Flow' not found" or empty observed data**
- Check your observed CSV file: must have columns named exactly `Date` and `Flow`.
- Check date format: must be `YYYY-MM-DD`.
- Check `obs_start_date` falls within your observed data range.

**"Satisfactory runs: 0"**
- No runs passed the thresholds. Options:
  - Relax thresholds (lower NSE/KGE, increase PBIAS tolerance)
  - Widen parameter ranges
  - Increase `n_simulations` for better coverage
  - Check if observed data is correctly aligned with simulated period

**"Error in run_swatplus(...)" or SWAT+ crashes**
- Test your SWAT+ model independently first (see Section 3).
- Check that `project_path` points to a valid TxtInOut folder.
- Reduce `n_thread` if you're running out of memory.

### Performance Tips

- **Phase 1 GSA** is fast (few runs, annual output). Start here to verify the setup works.
- **Phase 2-3 daily simulations** can be slow. Start with fewer simulations (e.g., 100) to test, then increase.
- **Disk space:** SWATrunR creates temporary run folders. The scripts clean them before each run, but check disk space if runs fail mid-execution.
- **Memory:** Large numbers of simulations can consume significant RAM. Monitor with Task Manager.

### Re-running a Step

Simply run the same step again — a new versioned folder (v2, v3...) will be created automatically. Previous results are preserved.

```r
# Re-run GSA Phase 2 (will create GSA_phase2_v2/)
PHASE_CONFIG <- "config/phase2.yaml"
source("run_gsa.R")
```

> **Note:** When Phase 2+ inherits ranges, it always picks the **latest version** folder
> from the previous phase. If you re-run Phase 1 CAL, Phase 2 will automatically use the
> newest results.

---

## Quick Reference

```r
# ---- Full pipeline ----
source("master.R")

# ---- Resume from step 3 ----
START_FROM <- 3
source("master.R")

# ---- Run single GSA ----
PHASE_CONFIG <- "config/phase1.yaml"
source("run_gsa.R")

# ---- Run single calibration ----
PHASE_CONFIG <- "config/phase2.yaml"
source("run_cal.R")

# ---- From command line ----
# Rscript run_gsa.R config/phase1.yaml
# Rscript run_cal.R config/phase1.yaml
# Rscript master.R

# ---- Load saved results ----
res <- readRDS("C:/path/to/your/swatplus/txtinout/HARD_phase2_v1/resultado_cal_phase2.rds")
```

---

## Citation

If you use this software, please cite it via the metadata in [`CITATION.cff`](CITATION.cff). After a tagged release, a versioned DOI is automatically minted by Zenodo and shown as a badge at the top of this README.

GitHub natively renders `CITATION.cff` as a "Cite this repository" button on the repository landing page (top right), with ready-to-paste BibTeX/APA strings.

---

## Contact

For technical questions, bug reports, or feature requests, please open an issue on GitHub: <https://github.com/ECRIPPEL/swatplus-hms-calibration/issues>.

Author: Elzon Cassio Rippel — ORCID [0000-0002-8391-4435](https://orcid.org/0000-0002-8391-4435).

---

## License

Released under the MIT License — see [`LICENSE`](LICENSE) for the full text.

This tool depends on [SWAT+](https://swat.tamu.edu/software/plus/) and [SWATrunR](https://chrisschuerz.github.io/SWATrunR/), which are distributed under their own licenses.
