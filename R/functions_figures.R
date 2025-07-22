out_path <- "plots"

method_vector <- c("Firearms", "Fixed wing", "Helicopter", "Snare", "Trap")

method_lookup <- tibble(
	method_idx = 1:5,
	method = method_vector
)

land_lookup <- tibble(
	land_idx = 1:3,
	land = c("Road density", "Ruggedness", "Canopy cover")
)


map_files2 <- function(dirs_vec, file_name) {
	get_files <- function(density_dir, file_name, node) {
		sim_results <- file.path(path, density_dir)
		ls <- read_rds(file.path(sim_results, file_name))
	}

	dirs_vec |>
		map(\(x) get_files(x, file_name)) |>
		list_rbind() |>
		mutate(start_density = as.factor(start_density))
}

obs_cols <- c("Observed" = "#1b9e77", "Not observed" = "#d95f02")

make_cats <- function(df) {
	metrics <- c(
		"mbias_density",
		"norm_bias_density",
		"rmse_density",
		"nm_rmse_density",
		"rmsle_density"
	)

	delta <- scores_rds |>
		filter(obs_flag == 1) |>
		group_by(property_id) |>
		mutate(delta = c(0, diff(PPNum))) |>
		ungroup() |>
		filter(delta > 0) |>
		select(property_id, PPNum, delta)

	cats <- scores_rds |>
		ungroup() |>
		select(
			-contains("abundance"),
			-mpe_density,
			-nm_rmse_density,
			-delta,
			all_of(metrics)
		) |>
		mutate(
			density_category = if_else(density <= 2, "Low", "Medium"),
			density_category = if_else(density >= 6, "High", density_category),
			density_category = if_else(density == 0, "Extinct", density_category),
			density_category = factor(
				density_category,
				levels = c("Extinct", "Low", "Medium", "High")
			),
			obs = if_else(obs_flag == 1, "Observed", "Not observed"),
			obs = factor(obs, levels = c("Observed", "Not observed"))
		) |>
		pivot_longer(cols = all_of(metrics), names_to = "metric") |>
		left_join(delta)
}

my_summary <- function(df) {
	df |>
		summarise(
			low = quantile(value, 0.05),
			q1 = quantile(value, 0.25),
			med = quantile(value, 0.5),
			q3 = quantile(value, 0.75),
			high = quantile(value, 0.95),
			n = n()
		)
}


potential_area <- function(params, m, metric) {
	if (m == "Firearms") {
		id <- 1
	}
	if (m == "Fixed wing") {
		id <- 2
	}
	if (m == "Helicopter") {
		id <- 3
	}
	if (m == "Snare") {
		id <- 4
	}
	if (m == "Trap") {
		id <- 5
	}

	log_rho <- params |>
		pull(paste0("log_rho[", id, "]"))

	beta1 <- params |>
		pull(paste0("beta1[", id, "]"))

	min_e <- effort_range |> filter(method == m) |> pull(min_e)
	max_e <- effort_range |> filter(method == m) |> pull(max_e)
	log_effort <- seq(log(min_e), log(max_e), length.out = 50)

	log_potential_area <- tibble()

	if (id <= 3) {
		for (i in seq_along(log_effort)) {
			log_e <- log_effort[i]

			lpa <- tibble(
				log_rho = log_rho,
				beta1 = beta1,
				log_effort_per = log_e,
				lpa = log_rho + log_e
			)

			log_potential_area <- bind_rows(log_potential_area, lpa)
		}
	} else {
		log_gamma <- params |>
			pull(paste0("log_gamma[", id - 3, "]"))
		p_unique <- params |>
			pull(paste0("p_mu[", id - 3, "]")) |>
			boot::inv.logit()
		n_trap_m1 <- if_else(id == 4, 9, 1)

		for (i in seq_along(log_effort)) {
			log_e <- log_effort[i]

			lpa1 <- log(pi) +
				(2 * (log_rho + log_e - log(exp(log_gamma) + exp(log_e)))) +
				log(1 + (p_unique * n_trap_m1))

			lpa <- tibble(
				log_effort_per = log_e,
				beta1 = beta1,
				lpa = lpa1
			)

			log_potential_area <- bind_rows(log_potential_area, lpa)
		}
	}

	if (metric == "potential area") {
		log_potential_area |>
			mutate(pa = exp(lpa), effort_per = exp(log_effort_per)) |>
			group_by(effort_per) |>
			summarise(
				`5%` = quantile(pa, 0.05),
				`25%` = quantile(pa, 0.25),
				`50%` = quantile(pa, 0.5),
				`75%` = quantile(pa, 0.75),
				`95%` = quantile(pa, 0.95)
			) |>
			ungroup() |>
			mutate(method = m, value = "Potential area")
	} else if (metric == "capture") {
		log_potential_area |>
			mutate(
				pa = exp(lpa),
				effort_per = exp(log_effort_per),
				beta1 = boot::inv.logit(beta1),
				pa = beta1 * pa
			) |>
			group_by(effort_per) |>
			summarise(
				`5%` = quantile(pa, 0.05),
				`25%` = quantile(pa, 0.25),
				`50%` = quantile(pa, 0.5),
				`75%` = quantile(pa, 0.75),
				`95%` = quantile(pa, 0.95)
			) |>
			ungroup() |>
			mutate(method = m, value = "Effective area")
	}
}

get_posterior <- function(df, y) {
	dfp <- df |>
		select(starts_with(y)) |>
		pivot_longer(cols = everything(), names_to = "node", values_to = "y")

	if (grepl("log", y)) {
		dfp <- dfp |> mutate(y = exp(y))
	}
	if (y == "log_nu") {
		dfp <- dfp |> mutate(y = y / 2)
	}
	if (grepl("p_mu", y) || grepl("beta1", y)) {
		dfp <- dfp |> mutate(y = boot::inv.logit(y))
	}

	dfp
}

join_summarise_methods <- function(df, df_method_names) {
	df |>
		left_join(df_method_names) |>
		group_by(method) |>
		summarise(
			`5%` = quantile(y, 0.05),
			`25%` = quantile(y, 0.25),
			`50%` = quantile(y, 0.5),
			`75%` = quantile(y, 0.75),
			`95%` = quantile(y, 0.95)
		) |>
		ungroup() |>
		mutate(method = factor(method, levels = c("Prior", method_vector)))
}

plot_post <- function(dfg, xlab, title) {
	dfg |>
		ggplot() +
		aes(x = `50%`, xmin = `5%`, xmax = `95%`, y = method) +
		geom_linerange(position = position_dodge(width = 0.5)) +
		geom_point(position = position_dodge(width = 0.5), size = 4) +
		labs(y = "Method", x = xlab, title = title, color = element_blank()) +
		theme_bw() +
		my_theme()
}
