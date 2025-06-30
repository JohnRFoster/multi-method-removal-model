#---------
#
# Workflow for fitting property-level Bayes model to MIS data
#
#---------

library(dplyr)
library(tidyr)
library(readr)
library(parallel)

# config_name <- "hpc_dev"
config_name <- "default"
config <- config::get(config = config_name)
interval <- config$interval

source("R/functions_data.R")
source("R/fit_mcmc_data.R")

# ===================================================
# Data ingest ----
# ===================================================

data_repo <- config$data_repo

# the below workflow is for completeness on how we winnowed the data
# starting from the entire MIS data
## -----------------------------------------------------------------------------
## MIS data workflow ----
# file <- file.path(data_repo, config$file_mis)
# dev <- config$dev

# data_mis <- get_data(file, interval, data_repo)

# ## observation covariates ----
# file <- file.path(data_repo, config$file_land)
# data_obs <- get_obs_covars(file)

# ## join MIS with observation covariates ----
# data_join <- left_join(data_mis, data_obs, by = join_by(county_code))

# ## filter missing states ----
# data_join2 <- data_join |>
# 	filter(!st_name %in% c("CALIFORNIA", "ALABAMA", "ARIZONA", "ARKANSAS"))

# targets::tar_assert_true(!any(is.na(data_join2$c_road_den)))
# targets::tar_assert_true(!any(is.na(data_join2$c_rugged)))
# targets::tar_assert_true(!any(is.na(data_join2$c_canopy)))

# ## join with farm bill properties ----
# data_farm_bill <- read_csv(file.path(
# 	data_repo,
# 	"All_FB_Agreements_long_2024-05-30.csv"
# ))
# farm_bill_properties <- data_farm_bill |>
# 	rename(alws_agrprop_id = propertyID) |>
# 	select(-agreement_name, -property_name) |>
# 	mutate(farm_bill = 1)

# data_final <- left_join(data_join2, farm_bill_properties) |>
# 	mutate(
# 		property = as.numeric(as.factor(propertyID)),
# 		county = as.numeric(as.factor(county_code))
# 	)

# data_for_nimble <- subset_data_for_development(data_final) |>
# 	mutate(primary_period = primary_period - min(primary_period) + 1)
# -----------------------------------------------------------------------------

# the data with the 105 properties used in the analysis is here
data_for_nimble <- read_csv("data/masked_mis_data.csv") |>
	mutate(property = propertyID, county = county_code)

params_check <- config$params_check
out_dir <- config$out_dir
files_in_out_dir <- list.files(out_dir)

monitors_add <- "N"
custom_samplers <- NULL

# run first fit
informed <- FALSE

finished <- prep_and_run_mcmc(
	informed = informed,
	post_path = NULL,
	data_repo = data_repo,
	dest_mcmc = dest_mcmc,
	dest_posterior = dest_posterior,
	df = data_for_nimble,
	monitors_add = monitors_add,
	custom_samplers = custom_samplers
)

source("R/check_mcmc.R")
