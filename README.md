# multi-method-removal-model

Code and analysis products for the manuscript **"Improving a removal model for evaluating density changes in a widespread invasive species"**.

The manuscript presents an integrated population/removal model for estimating wild pig (*Sus scrofa x domesticus*) density through time from management data collected with multiple removal methods. The repository contains the NIMBLE model, data-processing functions, empirical model-fitting workflow, simulation workflow, posterior summarizations, and figure scripts used for the paper.

## Citation

Foster, J. R., Pepin, K. M., Joseph, M. B., Tabak, M. A., & Miller, R. S. (2025). Improving a removal model for evaluating density changes in a widespread invasive species. bioRxiv. https://doi.org/10.1101/2025.10.16.680799

Currently in review at *Ecology.*

Code archive: https://doi.org/10.5281/zenodo.15784366

## Repository Contents

- [workflow_fit_to_data.R](workflow_fit_to_data.R): End-to-end workflow for preparing management data and fitting the property-level Bayesian model to the empirical MIS data.

- [workflow_simulations.R](workflow_simulations.R): HPC-oriented workflow for simulating removal datasets across starting densities, fitting the model, and writing simulation outputs.

- [config.yml](config.yml): Configuration profiles for local development and production/HPC runs.

- [R/nimble_removal_model.R](R/nimble_removal_model.R): Core NIMBLE model definition.

- [R/prep_nimble_data.R](R/prep_nimble_data.R) and [R/prep_nimble_simulation.R](R/prep_nimble_simulation.R): Data preparation (i.e. lists for NIMBLE) for data and simulated model fits.

- [R/fit_mcmc_data.R](R/fit_mcmc_data.R), [R/fit_mcmc_simulation.R](R/fit_mcmc_simulation.R), and [R/mcmc_parallel.R](R/mcmc_parallel.R): MCMC setup and parallel model fitting.

- [R/run_simulation.R](R/run_simulation.R): Simulation engine for generating properties, ecological dynamics, removals, model fits, and posterior summaries.

- [R/collate_abundance.R](R/collate_abundance.R), [R/collate_take.R](R/collate_take.R), [R/collate_parameters.R](R/collate_parameters.R), and [R/collate_scores.R](R/collate_scores.R): Collation and scoring of simulation output that was run on an HPC.

- [R/figure_*.R](R): Scripts for manuscript and supplemental figures.

## Software Requirements

This is an R project. The code uses parallel MCMC with NIMBLE and was developed for local and HPC execution.

Primary R package dependencies include:

- `nimble`
- `coda`
- `config`
- `dplyr`
- `tidyr`
- `readr`
- `purrr`
- `ggplot2`
- `ggpubr`
- `stringr`
- `boot`
- `assertthat`
- `targets`
- `lubridate`
- `testthat`
- `httpgd` for interactive plotting in some analysis scripts
- `boaR` for empirical figure preparation

## Reproducing the Analysis

Open [multi-method-removal-model.Rproj](multi-method-removal-model.Rproj) or start R from the repository root so that relative paths resolve correctly.

### 1. Configure paths and run settings

Review [config.yml](config.yml). The `default` profile is used by the empirical-data workflow and sets paths such as `data_repo`, `data_dir`, `out_dir`, and MCMC settings. The simulation workflow uses the `hpc_production` profile and expects a SLURM array environment for large production runs.

For local testing, reduce `n_iter`, `n_mcmc`, and possibly `n_pp` in [config.yml](config.yml). Production simulation fits use long MCMC chains and are not expected to finish quickly on a laptop.

### 2. Fit the empirical model

Run:

```r
source("workflow_fit_to_data.R")
```

This workflow:

1. Reads MIS removal data and land-cover covariates.
2. Filters and masks the empirical data used for model fitting.
3. Prepares constants, data, and initial values for the NIMBLE model.
4. Runs MCMC in parallel.
5. Checks MCMC output and writes posterior summaries.

### 3. Run simulations

Run:

```r
source("workflow_simulations.R")
```

The simulation workflow creates combinations of five starting densities (`0.3`, `1.475`, `2.65`, `3.825`, and `5`) and 200 simulation replicates per density. On an HPC system, each SLURM array task selects one row of this design through `SLURM_ARRAY_TASK_ID`; without that environment variable, the script defaults to the first task for local testing.

Each simulation replicate:

1. Bootstraps one-method and multi-method property histories from the empirical data.
2. Simulates ecological and removal dynamics.
3. Fits the NIMBLE model to the simulated data.
4. Checks MCMC convergence and effective sample sizes.
5. Writes known values, posterior samples, predictions, and diagnostics.

### 4. Collate simulation output

After simulation runs finish, use the collation scripts in [R](R):

```r
source("R/collate_abundance.R")
source("R/collate_take.R")
source("R/collate_parameters.R")
source("R/collate_scores.R")
```

These scripts summarize abundance recovery, removal predictions, parameter recovery, residuals, and error metrics. Collated products used by the figure scripts are stored in the `data/density_*` directories.

### 5. Make figures

The manuscript figure scripts are in [R](R):

```r
source("R/figure_simulationTimeseries.R")
source("R/figure_biasRMSLE.R")
source("R/figure_dataTimeseries.R")
source("R/figure_heatMaps.R")
source("R/figure_search.R")
source("R/figure_supplement.R")
```

Figures are written to [plots](plots). The scripts generally use the derived posterior and collated simulation files already present in [data](data).

## Model Overview

The model extends a removal/catch-effort framework to accommodate multiple management methods, including ground-shooting, fixed-wing aircraft, helicopter, snare, and trap removals. The ecological process describes latent abundance through primary periods, while the observation/removal process links method-specific effort, potential area searched, detection/capture probability, land-cover covariates, and removals. Simulations evaluate performance across a range of starting densities and population trajectories before applying the model to empirical wild pig management data.

## Notes for Reuse

- Run scripts from the repository root.
- Large MCMC runs are computationally intensive and should be run on a machine or cluster with sufficient CPU and memory.
- The HPC workflow assumes SLURM-style array jobs and parallel chains.
- Some paths in figure scripts point to `../data-store`; update these paths if your input data are stored elsewhere.
- The [data](data) directory contains masked and derived products intended to support reproducible downstream analyses without exposing original identifiers where masking was used.

## License

This repository is released under the MIT License. See [LICENSE](LICENSE) for details.
