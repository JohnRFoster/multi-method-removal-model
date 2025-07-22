library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(ggplot2)
library(ggpubr)
library(stringr)

source("R/functions_figures.R")

single <- c(
	"Firearms",
	"Fixed wing",
	"Helicopter",
	"Snares",
	"Traps"
)

scores_rds <- read_rds(file.path("data", "abundanceScoresByPrimaryPeriod.rds"))

cats <- make_cats(scores_rds)

quintiles <- cats |>
	filter(
		extinct == "Extant",
		methods_used %in% single
	) |>
	group_by(methods_used) |>
	summarise(
		`20%` = quantile(sum_effort_per_unit, 0.2),
		`40%` = quantile(sum_effort_per_unit, 0.4),
		`60%` = quantile(sum_effort_per_unit, 0.6),
		`80%` = quantile(sum_effort_per_unit, 0.8)
	)

get_quintile <- function(method, q) {
	quintiles |>
		filter(methods_used == method) |>
		pull(q)
}

create_effort_category <- function(m) {
	cats |>
		filter(obs_flag == 1, extinct == "Extant", methods_used == m) |>
		mutate(
			effort_category = if_else(
				sum_effort_per_unit <= get_quintile(m, "20%"),
				"0-20%",
				"80-100%"
			),
			effort_category = if_else(
				sum_effort_per_unit > get_quintile(m, "20%") &
					sum_effort_per_unit <= get_quintile(m, "40%"),
				"20-40%",
				effort_category
			),
			effort_category = if_else(
				sum_effort_per_unit > get_quintile(m, "40%") &
					sum_effort_per_unit <= get_quintile(m, "60%"),
				"40-60%",
				effort_category
			),
			effort_category = if_else(
				sum_effort_per_unit > get_quintile(m, "60%") &
					sum_effort_per_unit <= get_quintile(m, "80%"),
				"60-80%",
				effort_category
			),
		)
}

firearms <- create_effort_category("Firearms")
fixed <- create_effort_category("Fixed wing")
helicopter <- create_effort_category("Helicopter")
snares <- create_effort_category("Snares")
traps <- create_effort_category("Traps")

cats <- bind_rows(firearms, fixed, helicopter, snares, traps) |>
	mutate(
		methods_used = factor(methods_used, levels = single),
		effort_category = factor(
			effort_category,
			levels = c("0-20%", "20-40%", "40-60%", "60-80%", "80-100%")
		),
		methods_used = if_else(
			methods_used == "Firearms",
			"Sharpshooting",
			methods_used
		)
	)

my_tile <- function(df, m, x, y) {
	tmp <- df |>
		filter(metric == m) |>
		group_by(methods_used, .data[[x]], .data[[y]]) |>
		summarise(v = mean(value), n = n()) |>
		ungroup() |>
		filter(n >= 30)

	l1 <- abs(min(tmp$v))
	l2 <- abs(max(tmp$v))
	ll <- max(c(l1, l2)) * 1.05

	if (grepl("bias", m)) {
		limits <- round(c(-1 * ll, ll), 2)
		breaks <- c(limits[1], 0, limits[2])
		colors = c("#e9a3c9", "#f5f5f5", "#a1d76a")
	} else {
		limits <- round(c(0, ll), 2)
		breaks <- round(seq(0, ll, length.out = 3), 2)
		colors = c("navyblue", "darkmagenta", "darkorange1")
	}

	ggplot(tmp) +
		aes(x = .data[[x]], y = .data[[y]], fill = v) +
		geom_tile() +
		facet_grid(methods_used ~ .) +
		scale_fill_gradientn(colors = colors, limits = limits, breaks = breaks) +
		theme_bw() +
		theme(
			legend.position = "bottom",
			legend.text = element_text(size = 10, vjust = -1.5)
		)
}

x <- "density_category"
xlab <- "Density category"

y <- "effort_category"
ylab <- "Binned effort per unit"

m <- "norm_bias_density"
nb_lab <- "nBias"
g1 <- my_tile(cats, m, x, y) +
	labs(x = xlab, y = ylab, fill = nb_lab) +
	theme(strip.background.y = element_blank(), strip.text = element_blank())

m <- "nm_rmse_density"
nrmse_lab <- "nRMSE"
g2 <- my_tile(cats, m, x, y) +
	labs(x = xlab, y = "", fill = nrmse_lab) +
	theme(axis.text.y = element_blank())

ggarrange(g1, g2, ncol = 2, labels = "AUTO")

ggsave(
	file.path(out_path, "simulationNbiasNrmseHeatmap"),
	dpi = "retina",
	device = "jpeg",
	units = "in",
	width = 6,
	height = 6
)
