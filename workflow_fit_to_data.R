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
file <- file.path(data_repo, config$file_mis)
dev <- config$dev

data_mis <- get_data(file, interval, data_repo)

# ## observation covariates ----
file <- file.path(data_repo, config$file_land)
data_obs <- get_obs_covars(file)

# ## join MIS with observation covariates ----
data_join <- left_join(data_mis, data_obs, by = join_by(county_code))

# ## filter missing states ----
data_join2 <- data_join |>
	filter(!st_name %in% c("CALIFORNIA", "ALABAMA", "ARIZONA", "ARKANSAS"))

scales <- data_join2 |>
  select(st_name, county_code, rural.road.density, mean.ruggedness, mean.canopy.density) |>
  distinct() |>
  pivot_longer(cols = -c(st_name, county_code), names_to = "land", values_to = "value") |>
  group_by(land) |>
  reframe(mu = mean(value),
          sd = sd(value))

my_scale <- function(l, v){
  scales |> filter(land == l) |> pull(all_of(v))
}

center_scale <- function(x, mu, sd){
  (x - mu) / sd
}

c_road_den <- center_scale(
  data_join2$rural.road.density,
  my_scale("rural.road.density", "mu"),
  my_scale("rural.road.density", "sd"))

c_rugged <- center_scale(
  data_join2$mean.ruggedness,
  my_scale("mean.ruggedness", "mu"),
  my_scale("mean.ruggedness", "sd"))

c_canopy <- center_scale(
  data_join2$mean.canopy.density,
  my_scale("mean.canopy.density", "mu"),
  my_scale("mean.canopy.density", "sd"))

data_c <- data_join2 |>
  mutate(
    c_road_den = c_road_den,
    c_rugged = c_rugged,
    c_canopy = c_canopy
  )

# ## join with farm bill properties ----
data_farm_bill <- read_csv(file.path(
	data_repo,
	"All_FB_Agreements_long_2024-05-30.csv"
))
farm_bill_properties <- data_farm_bill |>
	rename(alws_agrprop_id = propertyID) |>
	select(-agreement_name, -property_name) |>
	mutate(farm_bill = 1)

data_final <- left_join(data_c, farm_bill_properties)

data_for_nimble <- subset_data_for_development(data_final) |>
	mutate(primary_period = primary_period - min(primary_period) + 1) |>
  mutate(
    property = as.numeric(as.factor(propertyID)),
    county = as.numeric(as.factor(county_code))
  )

cc <- c(
	"propertyID",
	"start_dates",
	"end_dates",
	"st_name",
	"county_code",
	"method",
	"trap_count",
	"take",
	"property_area_km2",
	"effort",
	"effort_per",
	"primary_period",
	"order",
	"n_survey",
	"observed_timestep",
	"c_road_den",
	"c_rugged",
	"c_canopy",
	"property",
	"county"
)

masked_data <- data_for_nimble |>
	select(all_of(cc)) |>
	mutate(
		propertyID = as.numeric(as.factor(propertyID)),
		property = as.numeric(as.factor(property)),
		county_code = as.numeric(as.factor(county_code))
	)

write_csv(masked_data, "data/masked_mis_data.csv")
write_csv(data_for_nimble, "data/originalFitDataWithIDs.csv")
write_csv(scales, "data/originalScaleMeanSD.csv")

# -----------------------------------------------------------------------------

# the data with the 105 properties used in the analysis is here
data_for_nimble <- read_csv(file.path(data_repo, "masked_mis_data.csv")) |>
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
