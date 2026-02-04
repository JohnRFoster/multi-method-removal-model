library(dplyr)
library(tidyr)
library(readr)
library(nimble)
library(ggplot2)
library(ggpubr)

set.seed(1256)

config_name <- "default"
config <- config::get(config = config_name)
interval <- config$interval

source("R/prep_nimble_data.R")
source("R/nimble_removal_model.R")
source("R/functions_nimble.R")

data_for_nimble <- read_csv("data/masked_mis_data.csv") |>
	mutate(
		property = propertyID,
		county = county_code,
		method_vec = as.numeric(as.factor(method))
	)

posterior_parameters <- read_rds("data/1_posterior/posteriorSamples.rds")
posterior_states <- read_rds("data/1_posterior/stateSamples.rds")

posterior_samples <- bind_cols(
	posterior_parameter_samples,
	posterior_state_samples
)

constants <- nimble_constants(
	df = data_for_nimble,
	interval = interval
)

data <- nimble_data(data_for_nimble)

inits <- nimble_inits(constants, data)

model <- nimbleModel(
	code = modelCode,
	constants = constants,
	data = data,
	inits = inits,
	calculate = TRUE
)

Cmodel <- compileNimble(Rmodel)

## Ensure we have the nodes needed to simulate new datasets
dataNodes <- Rmodel$getNodeNames(dataOnly = TRUE)
parentNodes <- Rmodel$getParents(dataNodes, stochOnly = TRUE)

## Ensure we have both data nodes and deterministic intermediates (e.g., lifted nodes)
simNodes <- Rmodel$getDependencies(parentNodes, self = FALSE)

# default MCMC configuration
mcmcConf <- configureMCMC(Rmodel, useConjugacy = TRUE)

control_rw <- list(
	adaptInterval = 100,
	adaptFactorExponent = 0.6
)

mcmcConf$removeSamplers("beta1")
mcmcConf$addSampler(
	target = "beta1",
	type = "RW_block",
	control = control_rw
)

mcmcConf$removeSamplers("log_nu")
mcmcConf$addSampler(
	target = "log_nu",
	type = "slice"
)

mcmcConf$addMonitors("N")

Rmcmc <- buildMCMC(mcmcConf)
Cmcmc <- compileNimble(Rmcmc)

# don't need these samples but run to get chains moving
samples <- runMCMC(Cmcmc, niter = 50000, nburnin = 5000)

nodes <- colnames(posterior_samples)

n_samp <- nrow(posterior_samples)
pp_samples <- matrix(NA, n_samp, nrow(model_data))

pb <- txtProgressBar(max = n_samp, style = 1)
for (i in 1:n_samp) {
	for (j in seq_along(nodes)) {
		Cmodel[[nodes[j]]] <- posterior_samples[i, nodes[j]]
	}
	Cmodel$simulate(simNodes, includeData = TRUE)
	pp_samples[i, ] <- Cmodel[["y"]]
	setTxtProgressBar(pb, i)
}
close(pb)

write_rds(
	pp_samples,
	file.path("data/1_posterior/posterior_predictive_samples.rds")
)

pp_samples <- read_rds("data/1_posterior/posterior_predictive_samples.rds")

colnames(pp_samples) <- paste0("y[", seq_len(nrow(data_for_nimble)), "]")

pps <- pp_samples |>
	as_tibble() |>
	mutate(sim = seq_len(n())) |>
	pivot_longer(cols = -sim, names_to = "node", values_to = "value")

info <- data_for_nimble |>
	select(propertyID, primary_period, method) |>
	mutate(node = paste0("y[", seq_len(nrow(data_for_nimble)), "]"))

pps_info <- pps |>
	left_join(info, by = "node")

pp_by_property <- pps_info |>
	group_by(propertyID, sim) |>
	summarise(predicted_take = sum(value), .groups = "drop")

obs_by_property <- data_for_nimble |>
	group_by(propertyID) |>
	summarise(observed_take = sum(take), .groups = "drop")

# The estimated p-value is just the proportion of these S simulations for
# which the test quantity equals or exceeds its realized value
# From Bayesian Data Analysis (3rd edition) Gelman et al. 2013

# minimum take at the property level
min_take_by_property <- obs_by_property |>
	filter(observed_take == min(observed_take))

min_properties <- min_take_by_property$propertyID
min_value <- unique(min_take_by_property$observed_take)

# posterior predictive disribution of the number of PP with zero take
prop_min_take_sim <- pp_by_property |>
	filter(propertyID %in% min_properties) |>
	mutate(above_min = ifelse(predicted_take >= min_value, 1, 0))

min_p_value <- sum(prop_min_take_sim$above_min) / nrow(prop_min_take_sim)

get_p_value <- function(df_obs, df_pred, func) {
	take_by_property <- df_obs |>
		filter(observed_take == func(observed_take))

	properties <- take_by_property$propertyID
	value <- unique(take_by_property$observed_take)

	# posterior predictive disribution of the number of PP with zero take
	prop_take_sim <- df_pred |>
		filter(propertyID %in% properties) |>
		mutate(above = ifelse(predicted_take >= value, 1, 0))

	p_val <- sum(prop_take_sim$above) / nrow(prop_take_sim)
	list(simulations = prop_take_sim, p_value = p_val, func_value = value)
}

p_min <- get_p_value(
	obs_by_property,
	pp_by_property,
	min
)

p_max <- get_p_value(
	obs_by_property,
	pp_by_property,
	max
)

p_med <- get_p_value(
	obs_by_property,
	pp_by_property,
	median
)

take_by_property_sd <- sd(obs_by_property$observed_take)

# posterior predictive disribution of the number of PP with zero take
prop_take_sim <- pp_by_property |>
	group_by(sim) |>
	summarise(sd_pred = sd(predicted_take)) |>
	mutate(above = ifelse(sd_pred >= take_by_property_sd, 1, 0))

p_val <- sum(prop_take_sim$above) / nrow(prop_take_sim)
p_sd <- list(
	simulations = prop_take_sim,
	p_value = p_val,
	func_value = take_by_property_sd
)

my_hist <- function(ls, f) {
	if (f == "min") {
		x_lab <- "T(y) = Property min(y)"
	} else if (f == "max") {
		x_lab <- "T(y) = Property max(y)"
	} else if (f == "median") {
		x_lab <- "T(y) = Property median(y)"
	} else if (f == "sd") {
		x_lab <- "T(y) = sd(y)"
		ls$simulations <- ls$simulations |>
			rename(predicted_take = sd_pred)
	}

	with(ls, {
		simulations |>
			ggplot() +
			aes(x = predicted_take) +
			geom_histogram() +
			geom_vline(xintercept = func_value, color = "red") +
			labs(
				title = paste0("p-value = ", round(p_value, 3)),
				x = x_lab,
				y = "Frequency"
			) +
			theme_bw()
	})
}

g <- list()
g[[1]] <- my_hist(p_min, "min")
g[[2]] <- my_hist(p_max, "max")
g[[3]] <- my_hist(p_med, "median")
g[[4]] <- my_hist(p_sd, "sd")


ggarrange(
	plotlist = g,
	ncol = 2,
	nrow = 2,
	labels = "AUTO"
)

ggsave(
	"plots/posterior_predictive_checks.jpg",
	width = 6,
	height = 6,
	units = "in",
	dpi = "retina"
)
