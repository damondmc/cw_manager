# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

PAWS (**P**ython **A**utomation for **W**eave **S**earches) orchestrates directed continuous gravitational wave (CW) searches using the `lalpulsar_Weave` executable. It handles: HTCondor DAG/SUB workflow generation, parameter space partitioning, post-search outlier analysis, follow-up stages, and upper limit determination.

## Installation

```bash
python -m pip install .
```

Dependencies: `numpy`, `astropy`, `scipy`, `pyfstat`, `matplotlib`, `tqdm`

## Container Build (OSG deployment)

```bash
# Build Apptainer/Singularity image for Open Science Grid
bash image.sh   # runs: apptainer build ../paws_gc.sif paws.def
```

## Configuration

Two YAML files drive everything:
- **`config.yaml`** — user/paths, HTCondor accounting, search parameters (`semi_mm`, `coh_mm`, `num_top_list`, `f0_band`, `nc_min`/`nc_max`, `cluster_n_spacing`, `followup_n_spacing`), and executable paths
- **Target YAML** (e.g., `gal.yaml`) — target name, sky position (`alpha`/`delta` in radians), uncertainties (`dalpha`/`ddelta`), and age

Both are loaded with `yaml.safe_load()` and passed as `config` and `target` dicts throughout the codebase.

## Architecture

### Core Primitives (`paws/definitions.py`)

Two functions used everywhere:
- `phase_param_name(order)` → `(["freq","f1dot",...], ["df","df1dot",...])` — returns parameter and bandwidth column names up to the given frequency derivative order (0–4)
- `task_name(target_name, stage, cohDay, order)` → `"{target}_{stage}_TCoh{cohDay}_O{order}"` — the standard naming convention for all files and DAG nodes

### Path Management (`paws/filepaths.py`)

`PathManager(config, target)` is the single source of truth for all file paths:
- `dag_file(freq, taskname, stage)` → `.dag` file path
- `condor_sub_file(freq, taskname, stage)` → `.sub` file path
- `weave_output_file(freq, taskname, job_index, stage)` → raw Weave FITS result path (under OSDF)
- `outlier_file(freq, taskname, stage, cluster=False)` → analyzed outlier FITS path (under `home_dir`)
- `sft_ensemble(freq)` → list of `osdf://` SFT URLs for H1 and L1 detectors (always OSDF; path driven by `sft_dir` in config)

### Parameter Generation (`paws/params/`)

All generators take a **model** instance that defines parameter space bounds:
- `PowerLawModel(nc_min, nc_max, tau)` — physics-based bounds derived from braking index and characteristic age; bounds vary with frequency
- `UniformModel(f1_lim, f2_lim, ...)` — flat rectangular bounds

Three generators:
- **`SearchParamGenerator(model, freq_deriv_order)`** — initial search: `generate_parameters(alpha, dalpha, delta, ddelta, f0_min, f0_max, ...)` returns a dict `{int(freq): fits.BinTableHDU}` of parameter tables (one per 1 Hz band)
- **`FollowUpParamGenerator(model)`** — follow-up: expands parameter ranges centered on previous-stage outliers; `generate_parameter(...)` returns a `fits.BinTableHDU`
- **`InjectionParamGenerator(model, ref_time, f0_band)`** — signal injections: `generate_parameters(...)` returns `(searchParamDict, injParamDict)` where each is a `{str(freq): fits.BinTableHDU}` dict

### Workflow Generation (`paws/workflow/`)

`WorkflowManager(config, target)` creates HTCondor DAG and SUB files:
- **`make_search_dag(..., tasks_per_job=10)`** — batches multiple Weave calls into a single Condor node via a bash wrapper script (`run_weave_batch_{stage}.sh`). Each node gets a task file listing raw Weave CLI arguments. Supports OSG file transfer (`use_osg=True`) and OSDF SFT access (`use_osdf=True`).
- **`make_upperlimit_dag(...)`** — single-node job that runs `upperlimit.py` (the Python upper limit script) via Python on OSG.

`writer.py` provides the low-level functions: `write_search_subfile(...)` and `write_search_dagfile(...)`.

### Pipeline (Local Execution, `paws/pipeline.py`)

For running Weave directly (not via HTCondor):
- `search_job(...)` and `injection_job(...)` — worker functions called via `multiprocessing.Pool`
- `determine_efficiency(...)` — runs parallel injection jobs, calls `ResultAnalysisManager.make_injection_outlier()`, and computes detection efficiency

### Result Analysis (`paws/analysis/`)

**`ResultAnalysisManager(config, target)`** (in `outlier.py`) is the unified result collector:
- `make_outlier(taskname, freq, mean2f_th, n_jobs, ...)` — reads all Weave FITS output files in parallel (via `ThreadPoolExecutor`), filters by `mean2F` threshold, optionally separates saturated bands, writes a multi-extension FITS outlier file, and optionally clusters results
- FITS output extensions: `{stage}_outlier`, `injection` (if injection run), `{stage}_sat_outlier` (if saturated), `info`, `non_sat_band` (search stage only)

**`clustering(data, spacing, cluster_n_spacing)`** (in `clustering.py`) — greedy clustering in phase-parameter space, sorting by `mean2F` loudness; returns `(centers_idx, cluster_size, cluster_member)`

**`SigmoidFitter`** (in `sigmoid.py`) — fits a sigmoid to h0 vs. detection efficiency data; `get_h_percentile(0.95)` returns h95 and its uncertainty. Used in `upperlimit.py` for iterative refinement until uncertainty ≤ 5%.

**`paws/analysis/astro.py`** — standalone astrophysical strain formulas (age strain limit, h0 from ellipticity/braking index, etc.)

**`paws/analysis/tools.py`** — `detection_stat_threshold(nTemp, nSeg)` computes the mean2F detection threshold from chi-squared statistics

### Upper Limit Script (`paws/upperlimit.py`)

A standalone script (with `argparse`) intended to run as a Condor job on OSG. Implements the iterative injection loop: starts with initial h0 estimates at 4 scale factors, fits a sigmoid, refines around h95 until convergence, then does a final confirmation run.

### I/O Utilities (`paws/io.py`)

- `make_dir(filenames)` — creates parent directories for file paths
- `get_spacing(data_file_path, freq_deriv_order)` — reads FITS header to extract per-dimension grid spacings (used for follow-up parameter generation and clustering)
- `get_bin_table(...)` / `get_header(...)` / `read_outlier_data(...)` — convenience wrappers for reading outlier FITS files
- Log parsers: `read_memory_usage`, `read_template_count`, `read_run_time`, `is_mismatch_exist`

## Execution Modes

The `work_in_local_dir` flag (bool, passed through most functions) controls whether file paths use absolute paths or just filenames — required when Condor transfers files to a worker node's local scratch space. The `cluster` flag on `outlier_file()` switches between unclustered (`_outlier.fts`) and clustered (`_outlier_clustered.fts`) output paths.

## Key Data Conventions

- All parameter tables use `astropy.Table` / `fits.BinTableHDU`
- Search parameter columns: `freq`, `df`, `f1dot`, `df1dot`, ... (from `phase_param_name`)
- Sky columns always present: `alpha`, `dalpha`, `delta`, `ddelta`
- Injection parameter column names capitalize `Alpha`, `Delta` and include `refTime`, `aPlus`, `aCross`, `psi`
- Frequency derivative order max is 4 (f0 through f4dot)
- Frequency bands are indexed by integer Hz (e.g., `params[97]` for the 97–98 Hz band)
