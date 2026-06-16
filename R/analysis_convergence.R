library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(ggplot2)
library(httpgd)
library(ggpubr)

modelData <- read_rds("data/1_posterior/posterior_predictive_samples.rds")

map_files2 <- function(file_name) {
	get_files <- function(density_dir, file_name) {
		sim_results <- file.path(path, density_dir)
		ls <- read_rds(file.path(sim_results, file_name))
	}

	paste0("density_", c(0.3, 1.475, 2.65, 3.825, 5)) |>
		map(\(x) get_files(x, file_name)) |>
		list_rbind() |>
		mutate(start_density = as.factor(start_density))
}

path <- "../pigs-simulation/analysis/betaSurvival_uniqueAreaTrapSnare"

fname <- "all_sim_converged.rds"
all_sim_converged <- map_files2(fname)

sim_info <- all_sim_converged |>
	select(simulation, start_density, sim_died, bad_mcmc)

all_sim_converged |>
	filter(!sim_died) |>
	group_by(start_density) |>
	summarise(
		n_converged = sum(converged),
		n_not_converged = sum(!converged),
		n_sims = n(),
		prop_converged = n_converged / n_sims
	)


all_sim_converged |>
	filter(!sim_died) |>
	summarise(
		n_converged = sum(converged),
		n_not_converged = sum(!converged),
		n_sims = n(),
		prop_converged = n_converged / n_sims
	)

fname <- "all_psrf.rds"
all_psrf <- map_files2(fname) |>
	left_join(sim_info, by = c("simulation", "start_density"))


all_psrf |>
	filter(!node_converged, node_names != "psi_phi") |>
	group_by(node_names, start_density) |>
	count(node_converged) |>
	arrange(desc(n))

all_psrf |>
	filter(!node_converged, node_names != "psi_phi") |>
	group_by(node_names) |>
	count(node_converged) |>
	arrange(desc(n))

all_psrf |>
	filter(node_names == "log_rho[4]", !node_converged)


fname <- "all_known_parameters.rds"
all_params <- map_files2(fname) |>
	left_join(sim_info, by = c("simulation", "start_density")) |>
	left_join(all_psrf)

all_params |>
	filter(
		!bad_mcmc,
		node_names %in%
			c("log_rho[4]", "log_rho[5]", "log_gamma[1]", "log_gamma[2]")
	) |>
	ggplot() +
	aes(x = actual, y = `Point est.`, color = node_converged) +
	geom_point() +
	facet_wrap(~node_names, scales = "free") +
	theme_bw()


ts_wide <- all_params |>
	filter(
		!bad_mcmc,
		node_names %in%
			c("log_rho[4]", "log_rho[5]", "log_gamma[1]", "log_gamma[2]")
	) |>
	select(
		-position,
		-`Point est.`,
		-`Upper C.I.`,
		-node_converged,
		-method_idx,
		-method
	) |>
	pivot_wider(names_from = node_names, values_from = actual) |>
	left_join(all_sim_converged)

plot_pair <- function(df, x, y) {
	df |>
		ggplot() +
		aes(x = .data[[x]], y = .data[[y]], color = converged) +
		geom_point() +
		theme_bw()
}

ts_wide |> plot_pair("log_rho[4]", "log_rho[5]")
ts_wide |> plot_pair("log_rho[4]", "log_gamma[1]")
ts_wide |> plot_pair("log_rho[4]", "log_gamma[2]")

ts_wide |> plot_pair("log_rho[5]", "log_gamma[1]")
ts_wide |> plot_pair("log_rho[5]", "log_gamma[2]")

ts_wide |> plot_pair("log_gamma[1]", "log_gamma[2]")
