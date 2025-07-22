library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(ggplot2)
library(stringr)

source("R/functions_figures.R")

analysis_dir <- "../pigs-simulation/analysis"
model_dir <- "betaSurvival_uniqueAreaTrapSnare"
path <- file.path(analysis_dir, model_dir)
density_dirs <- paste0("density_", c(0.3, 1.475, 2.65, 3.825, 5))

message("Get abundance summaries")
f_name <- "abundance_summaries.rds"
df <- map_files2(density_dirs, f_name)

# create property ID for easier joining
property_ids <- df |>
	select(start_density, simulation, property) |>
	distinct() |>
	mutate(property_id = paste(start_density, simulation, property, sep = "-"))

density <- left_join(df, property_ids) |>
	select(
		property_id,
		PPNum,
		low_density,
		med_density,
		high_density,
		density,
		obs_flag,
		start_density
	)

message("Get take summaries")
f_name <- "take_summaries.rds"
df <- map_files2(density_dirs, f_name) |>
	left_join(property_ids)

n_methods_pp <- df |>
	select(property_id, PPNum, method) |>
	distinct() |>
	group_by(property_id, PPNum) |>
	count() |>
	rename(n_methods_used = n)

n_return <- df |>
	left_join(property_ids) |>
	select(property_id, PPNum, method) |>
	distinct() |>
	pivot_wider(names_from = method, values_from = method) |>
	unite(method, -c(property_id, PPNum), sep = ", ", na.rm = TRUE) |>
	group_by(property_id, method) |>
	mutate(return_interval = c(0, diff(PPNum))) |>
	ungroup() |>
	rename(methods_used = method) |>
	left_join(n_methods_pp)

take_joint_return <- left_join(df, n_return)

take <- take_joint_return |>
	select(property_id, PPNum, sum_take, property_area, methods_used) |>
	mutate(sum_take = sum_take / property_area) |>
	unique()

density_take <- left_join(density, take)

extant_properties <- density_take |>
	group_by(property_id) |>
	filter(PPNum == max(PPNum)) |>
	ungroup() |>
	filter(density > 0) |>
	pull(property_id)

all_slopes <- tibble()
all_props <- extant_properties

for (j in seq_along(all_props)) {
	tmp <- density_take |> filter(property_id == all_props[j])

	if (nrow(tmp) == 1) {
		next
	}

	m <- lm(density ~ PPNum, data = tmp)
	ypred <- predict(m, newdata = tibble(PPNum = tmp$PPNum))
	slope <- summary(m)$coefficients[2]

	sp <- tibble(
		property_id = all_props[j],
		lambda = slope
	)

	all_slopes <- bind_rows(all_slopes, sp)
}

density_recovered <- density_take |>
	mutate(
		recovered = if_else(
			density >= low_density & density <= high_density,
			1,
			0
		)
	)

thresh <- 0.05
directions <- density_recovered |>
	filter(property_id %in% extant_properties) |>
	left_join(all_slopes) |>
	mutate(
		direction = if_else(lambda < (-1 * thresh), "Decreasing", "Stable"),
		direction = if_else(lambda > thresh, "Increasing", direction)
	)

table(directions$direction)
summary(directions$lambda)

directions |>
	group_by(direction) |>
	summarise(per_recovered = round(sum(recovered) / n() * 100, 2))

id1 <- "3.825-89-94" # steady decline
id2 <- "5-54-94" # steady incline
id3 <- "3.825-118-61" # steady

plot_df <- directions |>
	filter(property_id %in% c(id1, id2, id3))

method_shapes <- plot_df |>
	filter(obs_flag == 1) |>
	select(methods_used) |>
	distinct() |>
	mutate(shape = 1:n())

plot_df |>
	mutate(methods_used = as.factor(methods_used)) |>
	group_by(property_id) |>
	mutate(
		obs_type = if_else(obs_flag == 1, "Observed", "Not observed"),
		time = 1:n()
	) |>
	ungroup() |>
	ggplot() +
	aes(x = time) +

	geom_ribbon(
		aes(ymin = low_density, ymax = high_density, fill = "95% CI"),
		alpha = 0.6
	) +
	scale_fill_manual(values = "gray") +

	geom_line(aes(y = med_density, linetype = "Median")) +

	geom_point(aes(y = density, color = obs_type), size = 1) +
	scale_color_manual(values = obs_cols, drop = FALSE) +

	geom_point(aes(y = sum_take, shape = methods_used), size = 1) +
	scale_shape_manual(
		breaks = method_shapes$methods_used,
		values = method_shapes$shape,
		drop = FALSE
	) +

	labs(
		x = "Time",
		y = expression("Density (pigs / " ~ km^2 ~ ")"),
		color = "Known density",
		shape = "Method(s) used",
		linetype = element_blank(),
		shape = element_blank(),
		fill = element_blank()
	) +
	facet_wrap(~direction) +
	# facet_grid(direction ~ .) +
	theme_bw() +
	theme(
		axis.text = element_text(size = 8),
		axis.title = element_text(size = 10),
		strip.text = element_text(size = 10),
		legend.title = element_text(size = 8),
		legend.text = element_text(size = 8)
	)

ggsave(
	file.path("plots", "simulationTimeSeries"),
	dpi = "retina",
	device = "jpeg",
	units = "cm",
	width = 16,
	height = 12
)
