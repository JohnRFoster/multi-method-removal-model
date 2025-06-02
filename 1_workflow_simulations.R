# --------------------------------------------------------------------
#
# Workflow script for invasive species removal simulations
#
# Author:
# John Foster
#
# General workflow (conducted with run_simulation):
#
# 1. Load config and MIS data
#
# 2. run_simulation
#   - Simulate/bootstrap data
#     - 1-method properties
#     - 2- to 5-method properties
#   - Simulate eco/take dynamics
#   - Fit MCMC
#   - Check MCMC
# 3. Summarize output
#
# --------------------------------------------------------------------

config_name <- "dev"
config <- config::get(config = config_name)

library(nimble)
library(parallel)
library(coda)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)


# -----------------------------------------------------------------
# Load MIS data ----
# -----------------------------------------------------------------
message("MIS data intake")
data_dir <- config$data_dir
df <- read_rds(file.path(data_dir, "MIS_4weekPP.rds"))
df <- df |>
  filter(property.size >= 1.8) |> # median home range size from Kay et al. (2017)
  select(-method) |>
  rename(property = agrp_prp_id,
         method = Method)

# -----------------------------------------------------------------
# Run simulation ----
# -----------------------------------------------------------------

# meant to be run in parallel on an hpc
# the array number gets assigned to task_id

args <- commandArgs(trailingOnly = TRUE)
task_id <- args[1]
if(is.na(task_id)) task_id <- 1 # for testing locally
message("Task ID: ", task_id)

source("R/run_simulation.R")
out_list <- run_simulation(config, df, task_id)

# -----------------------------------------------------------------
# Summarize output ----
# -----------------------------------------------------------------

# source("R/collate_mcmc_output.R")
# collate_mcmc_output(config, sim)

message("\nRun time: ")
print(Sys.time() - start)

message("\n\nDONE!")

