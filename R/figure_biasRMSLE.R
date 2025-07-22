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

cats_extant <- cats |>
	filter(density_category != "Extinct")
cats_extinct <- cats |>
	filter(density_category == "Extinct") |>
	group_by(property_id) |>
	filter(PPNum == min(PPNum)) |>
	ungroup()

cats_summary <- bind_rows(cats_extant, cats_extinct) |>
	filter(metric %in% c("mbias_density", "rmsle_density")) |>
	group_by(density_category, metric, obs) |>
	my_summary() |>
	ungroup()

my_linerange <- function(df, m, x, w = 0.4) {
	df |>
		filter(metric == m) |>
		ggplot() +
		aes(x = .data[[x]], color = obs) +
		geom_linerange(
			aes(ymin = low, ymax = high),
			position = position_dodge(width = w)
		) +
		geom_linerange(
			aes(ymin = q1, ymax = q3),
			linewidth = 2,
			position = position_dodge(width = w)
		) +
		geom_point(aes(y = med), size = 4, position = position_dodge(width = w)) +
		scale_color_manual(values = obs_cols) +
		# scale_color_brewer(type = "qual", palette = 2) +
		labs(color = element_blank()) +
		theme_bw()
}

bias_lab <- expression("Bias (pigs / " ~ km^2 ~ ")")
rmsle_lab <- "RMSLE"

x <- "density_category"
xlab <- "Density"

obs_cols <- c("Observed" = "#1b9e77", "Not observed" = "#d95f02")

g1 <- my_linerange(cats_summary, "mbias_density", x) +
	labs(x = xlab, y = bias_lab) +
	geom_hline(yintercept = 0, linetype = "dashed")
g2 <- my_linerange(cats_summary, "rmsle_density", x) +
	labs(x = xlab, y = rmsle_lab)

ggarrange(
	g1,
	g2,
	ncol = 2,
	nrow = 1,
	labels = "AUTO",
	common.legend = TRUE,
	legend = "bottom"
)

ggsave(
	file.path(out_path, "biasRMSLEDensity"),
	dpi = "retina",
	device = "jpeg",
	units = "in",
	width = 6,
	height = 4
)
