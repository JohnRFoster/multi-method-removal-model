# code for manuscript figures

library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(stringr)
library(ggplot2)
library(ggpubr)

source("R/functions_collate.R")

out_path <- "plots"
density_dirs <- paste0("density_", c(0.3, 1.475, 2.65, 3.825, 5))

map_files2 <- function(dirs_vec, file_name) {
	get_files <- function(density_dir, file_name, node) {
		sim_results <- file.path("data", density_dir)
		ls <- read_rds(file.path(sim_results, file_name))
	}

	dirs_vec |>
		map(\(x) get_files(x, file_name)) |>
		list_rbind() |>
		mutate(start_density = as.factor(start_density))
}

f_name <- "all_param_samples.rds"
samples <- map_files2(density_dirs, f_name)

f_name <- "method_parameter_lookup.rds"
known_params <- map_files2(density_dirs, f_name)

f_name <- "all_beta_p.rds"
beta <- map_files2(density_dirs, f_name)

# figure s1 beta residuals----

my_theme <- function() {
	theme(
		axis.title = element_text(size = 14),
		axis.text = element_text(size = 12),
		strip.text = element_text(size = 12)
	)
}

my_linerange <- function(df) {
	ggplot(df) +
		aes(x = med, xmin = low, xmax = high, y = method) +
		geom_linerange(linewidth = 2) +
		geom_linerange(aes(xmin = q1, xmax = q3), linewidth = 4) +
		geom_point(size = 7) +
		geom_point(size = 5, color = "white") +
		geom_vline(xintercept = 0, linetype = "dashed") +
		labs(x = "Residual", y = "Method") +
		theme_bw() +
		my_theme()
}

my_summary <- function(df) {
	df |>
		reframe(
			low = quantile(value, 0.025),
			q1 = quantile(value, 0.25),
			med = quantile(value, 0.5),
			q3 = quantile(value, 0.75),
			high = quantile(value, 0.975),
			mu = mean(value),
			sd = sd(value)
		)
}

method_table <- known_params |>
	select(idx, method) |>
	distinct() |>
	rename(method_idx = idx)

land_cover_table <- tibble(
	position = 1:4,
	landCover = c("Intercept", "Road density", "Ruggedness", "Canopy cover")
)

beta_p_long <- samples |>
	select_pivot_longer("beta_p") |>
	mutate(
		method_idx = as.numeric(str_extract(node, "(?<=\\[)\\d")),
		position = as.numeric(str_extract(node, "(?<=\\, )\\d")) + 1
	)

beta_p_known <- left_join(beta_p_long, beta) |>
	mutate(value = value - actual)

beta1_long <- samples |>
	select_pivot_longer("beta1") |>
	mutate(
		method_idx = as.numeric(str_extract(node, "(?<=\\[)\\d")),
		position = 1
	)

beta_1_known <- left_join(beta1_long, beta) |>
	mutate(value = value - actual)

b1_summary <- beta_1_known |>
	group_by(method_idx) |>
	my_summary()

beta_p_known |>
	group_by(method_idx, position) |>
	my_summary() |>
	ungroup() |>
	bind_rows(mutate(b1_summary, position = 1)) |>
	left_join(method_table) |>
	left_join(land_cover_table) |>
	mutate(
		landCover = factor(
			landCover,
			levels = c("Intercept", "Road density", "Ruggedness", "Canopy cover")
		)
	) |>
	my_linerange() +
	facet_wrap(~landCover)

ggsave(
	file.path(out_path, "simulationBetas-v2.jpeg"),
	dpi = "retina",
	device = "jpeg",
	units = "cm",
	width = 16,
	height = 12
)

# figure s2 gamma residuals and median vs known ----

gamma_long <- samples |>
	select_pivot_longer("log_gamma[") |>
	mutate(idx = as.numeric(str_extract(node, "(?<=\\[)\\d")), idx = idx + 3)

gamma_known <- left_join(gamma_long, known_params) |>
	select(-rho, -p_unique) |>
	mutate(log_gamma = log(gamma), residual = value - log_gamma)

gH <- gamma_known |>
	select(simulation, start_density, method, gamma, log_gamma) |>
	distinct()

g1 <- gamma_known |>
	group_by(method) |>
	mutate(value = residual) |>
	my_summary() |>
	ungroup() |>
	my_linerange() +
	theme(
		axis.title = element_text(size = 10),
		axis.text = element_text(size = 8),
		strip.text = element_text(size = 12)
	)

g2 <- gamma_known |>
	group_by(simulation, start_density, method) |>
	mutate(value = exp(value)) |>
	my_summary() |>
	ungroup() |>
	left_join(gH) |>
	ggplot() +
	aes(x = gamma, y = med) +
	geom_point() +
	geom_smooth(method = "lm") +
	geom_abline(intercept = 0, slope = 1) +
	facet_grid(method ~ ., scales = "free_y") +
	labs(x = "Known parameter value", y = "Posterior median") +
	theme_bw() +
	theme(
		axis.title = element_text(size = 10),
		axis.text = element_text(size = 8),
		strip.text = element_text(size = 12)
	)

ggarrange(g1, g2, nrow = 1, labels = "AUTO")

ggsave(
	file.path(out_path, "simulationGammas-v2.jpeg"),
	dpi = "retina",
	device = "jpeg",
	units = "cm",
	width = 16,
	height = 12
)

# figure s3 rho residuals and median vs known ----
rho_long <- samples |>
	select_pivot_longer("log_rho[") |>
	mutate(idx = as.numeric(str_extract(node, "(?<=\\[)\\d")))

rho_known <- left_join(rho_long, known_params) |>
	select(-gamma, -p_unique) |>
	mutate(log_rho = log(rho), residual = value - log_rho)

r1 <- rho_known |>
	group_by(method) |>
	mutate(value = residual) |>
	my_summary() |>
	ungroup() |>
	my_linerange()

rH <- rho_known |>
	select(simulation, start_density, method, rho, log_rho) |>
	distinct()

r2 <- rho_known |>
	group_by(simulation, start_density, method) |>
	mutate(value = exp(value)) |>
	my_summary() |>
	ungroup() |>
	left_join(rH) |>
	ggplot() +
	aes(x = rho, y = med) +
	geom_point(size = 0.5) +
	geom_smooth(method = "lm") +
	geom_abline(intercept = 0, slope = 1) +
	facet_wrap(~method, scales = "free") +
	labs(x = "Known parameter value", y = "Posterior median") +
	theme_bw()

ggarrange(r1, r2, nrow = 1, labels = "AUTO")

ggsave(
	file.path(out_path, "simulationRhos-v2.jpeg"),
	dpi = "retina",
	device = "jpeg",
	units = "cm",
	width = 16,
	height = 12
)


# figure s4 omega residuals and median vs known ----
p_long <- samples |>
	select_pivot_longer("p_mu[") |>
	mutate(idx = as.numeric(str_extract(node, "(?<=\\[)\\d")), idx = idx + 3)

p_known <- left_join(p_long, known_params) |>
	select(-gamma, -rho) |>
	mutate(logit_p = boot::logit(p_unique), residual = value - logit_p)

p1 <- p_known |>
	group_by(method) |>
	mutate(value = residual) |>
	my_summary() |>
	ungroup() |>
	my_linerange()

pH <- p_known |>
	select(simulation, start_density, method, p_unique, logit_p) |>
	distinct()

p2 <- p_known |>
	group_by(simulation, start_density, method) |>
	mutate(value = boot::inv.logit(value)) |>
	my_summary() |>
	ungroup() |>
	left_join(pH) |>
	ggplot() +
	aes(x = p_unique, y = med) +
	geom_point(size = 0.25) +
	geom_smooth(method = "lm") +
	geom_abline(intercept = 0, slope = 1) +
	facet_grid(method ~ start_density, scales = "free") +
	labs(x = "Known parameter value", y = "Posterior median") +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggarrange(p1, p2, nrow = 1, labels = "AUTO", widths = c(1, 2))

ggsave(
	file.path(out_path, "simulationOmegas-v2.jpeg"),
	dpi = "retina",
	device = "jpeg",
	units = "cm",
	width = 16,
	height = 12
)

# figure s5 vital rates ----

phi_long <- samples |>
	select_pivot_longer("phi_mu") |>
	mutate(actual = 0.78, value = value - actual)

g1 <- phi_long |>
	my_summary() |>
	mutate(method = "Survival") |>
	my_linerange() +
	labs(y = "", title = "Survival") +
	theme(axis.text.y = element_blank())

psi_long <- samples |>
	select_pivot_longer("psi_phi") |>
	mutate(actual = 5, value = value - actual)

samples |>
	select_pivot_longer("psi_phi") |>
	group_by(start_density, simulation) |>
	my_summary() |>
	ungroup() |>
	filter(low <= 5 & high >= 5)

samples |>
	select_pivot_longer("phi_mu") |>
	group_by(start_density, simulation) |>
	my_summary() |>
	ungroup() |>
	filter(low <= 0.78 & high >= 0.78)

g2 <- psi_long |>
	group_by(start_density) |>
	my_summary() |>
	ungroup() |>
	mutate(method = start_density) |>
	my_linerange() +
	labs(y = "Starting density", title = "Shrinkage")

nu_long <- samples |>
	select_pivot_longer("log_nu") |>
	mutate(actual = log(5.290323), value = value - actual)

g3 <- nu_long |>
	my_summary() |>
	mutate(method = "Fecundity") |>
	my_linerange() +
	labs(y = "", title = "Fecundity") +
	theme(axis.text.y = element_blank())

gg1 <- ggarrange(g1, g3, nrow = 1, ncol = 2, labels = c("A", "B"))
g4 <- blank <- ggplot() + geom_blank() + theme_void()
gg2 <- ggarrange(g2, g4, nrow = 1, ncol = 2, widths = c(100, 1), labels = "C")

ggarrange(gg1, gg2, nrow = 2)

ggsave(
	file.path(out_path, "demographicParams-v2.jpeg"),
	dpi = "retina",
	device = "jpeg",
	units = "cm",
	width = 16,
	height = 12
)


# figure s6 effect of land cover ----
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

plot_post(beta_summary, xlab, title) +
	facet_wrap(~land, scales = "free") +
	geom_vline(xintercept = 0, linetype = "dashed") +
	theme(legend.position = "none")
