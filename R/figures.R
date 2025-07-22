# code for manuscript figures

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(stringr)

# figure s1 beta residuals----

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
		summarise(
			low = quantile(value, 0.05),
			q1 = quantile(value, 0.25),
			med = quantile(value, 0.5),
			q3 = quantile(value, 0.75),
			high = quantile(value, 0.95),
			mu = mean(value),
			sd = sd(value)
		)
}

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

# figure s2 gamma residuals and median vs known ----
gamma_known |>
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

# figure s3 rho residuals and median vs known ----
rho_known |>
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

# figure s4 omega residuals and median vs known ----
p_known |>
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

# figure s5 vital rates ----
g1 <- phi_long |>
	my_summary() |>
	mutate(method = "Survival") |>
	my_linerange() +
	labs(y = "", title = "Survival") +
	theme(axis.text.y = element_blank())

psi_long <- samples |>
	select_pivot_longer("psi_phi") |>
	mutate(actual = 5, value = value - actual)

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
