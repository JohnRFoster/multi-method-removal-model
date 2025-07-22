library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(ggplot2)
library(ggpubr)
library(stringr)

source("R/functions_figures.R")

top_dir <- "data"
if_dir <- "1_posterior"
posterior_path <- file.path(top_dir, if_dir, "posteriorSamples.rds")
params <- read_rds(posterior_path)

method_names <- tibble(
	node = colnames(params),
	name = NA,
	method_idx = NA,
	land_idx = NA
) |>
	mutate(
		method_idx = as.numeric(str_extract(node, "(?<=\\[)\\d")),
		method_idx = if_else(
			grepl("log_gamma", node) | grepl("p_mu", node),
			method_idx + 3,
			method_idx
		),
		land_idx = as.numeric(str_extract(node, "(?<=\\,\\s)\\d"))
	) |>
	left_join(method_lookup) |>
	left_join(land_lookup)

posterior_path <- file.path(top_dir, if_dir, "modelData.rds")
data <- read_rds(posterior_path) |>
	mutate(
		method = if_else(method == "FIREARMS", "Firearms", method),
		method = if_else(method == "FIXED WING", "Fixed wing", method),
		method = if_else(method == "HELICOPTER", "Helicopter", method),
		method = if_else(method == "SNARE", "Snare", method),
		method = if_else(method == "TRAPS", "Trap", method)
	)

effort_range <- data |>
	group_by(method) |>
	summarise(
		min_e = min((effort_per)),
		max_e = max((effort_per)),
		min_t = min(trap_count - 1),
		mean_t = mean(trap_count - 1),
		med_t = median(trap_count - 1),
		max_t = max(trap_count - 1)
	) |>
	ungroup()

pa_firearms_post1 <- potential_area(params, "Firearms", "potential area") |>
	mutate(distribution = "Posterior")

pa_fixedwing_post1 <- potential_area(params, "Fixed wing", "potential area") |>
	mutate(distribution = "Posterior")

pa_helicopter_post1 <- potential_area(params, "Helicopter", "potential area") |>
	mutate(distribution = "Posterior")

pa_snare_post1 <- potential_area(params, "Snare", "potential area") |>
	mutate(distribution = "Posterior")

pa_traps_post1 <- potential_area(params, "Trap", "potential area") |>
	mutate(distribution = "Posterior")


my_theme <- function(s = 14) {
	theme(
		title = element_text(size = s + 4),
		axis.title = element_text(size = s),
		axis.text = element_text(size = s - 2),
		strip.text = element_text(size = s - 2),
		legend.text = element_text(size = s)
	)
}

df <- bind_rows(
	pa_firearms_post1,
	pa_fixedwing_post1,
	pa_helicopter_post1,
	pa_snare_post1,
	pa_traps_post1
)

dist_colors <- c("Posterior" = "gray")

g1 <- df |>
	filter(value == "Potential area", distribution == "Posterior") |>
	mutate(method = if_else(method == "Firearms", "Sharpshooting", method)) |>
	ggplot() +
	aes(
		x = effort_per,
		y = `50%`,
		ymin = `5%`,
		ymax = `95%`,
		fill = distribution,
		linetype = distribution
	) +
	geom_ribbon(alpha = 0.5) +
	geom_ribbon(aes(ymin = `25%`, ymax = `75%`), alpha = 1) +
	geom_line() +
	facet_wrap(~method, scales = "free") +
	labs(
		x = "Effort per unit deployed",
		y = bquote("Search area (" ~ km^2 ~ ")")
	) +
	scale_fill_manual(values = dist_colors) +
	theme_bw() +
	theme(legend.position = "none")

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

beta1_1 <- get_posterior(params, "beta1")
beta1_summary <- join_summarise_methods(beta1_1, method_names)

g2 <- beta1_summary |>
	mutate(method = if_else(method == "Firearms", "Sharpshooting", method)) |>
	ggplot() +
	aes(x = `50%`, y = method) +
	geom_linerange(aes(xmin = `5%`, xmax = `95%`), linewidth = 0.5) +
	geom_linerange(aes(xmin = `25%`, xmax = `75%`), linewidth = 2) +
	geom_point(size = 3.5, color = "black") +
	geom_point(size = 2, color = "white") +
	labs(x = "Base removal rate", y = element_blank()) +
	coord_cartesian(xlim = c(0, 1)) +
	# coord_flip() +
	theme_bw()

ggarrange(g1, g2, nrow = 2, labels = c("A", "B"), heights = c(2, 1))

ggsave(
	file.path(out_path, "methodEfficiency"),
	dpi = "retina",
	device = "jpeg",
	units = "in",
	width = 6,
	height = 4
)
