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

data_for_nimble <- read_csv(file.path("../data-store/masked_mis_data.csv")) |>
	mutate(
		property = propertyID,
		county = county_code,
		method = if_else(method == "FIREARMS", "Firearms", method),
		method = if_else(method == "FIXED WING", "Fixed wing", method),
		method = if_else(method == "HELICOPTER", "Helicopter", method),
		method = if_else(method == "SNARE", "Snare", method),
		method = if_else(method == "TRAPS", "Trap", method)
	)

all_pp <- boaR:::create_all_primary_periods(data_for_nimble) |>
	select(-timestep)

property_info <- data_for_nimble |>
	select(property, primary_period, property_area_km2) |>
	left_join(all_pp) |>
	distinct()

posterior_path <- file.path(top_dir, if_dir, "stateSamples.rds")
abundance <- read_rds(posterior_path)

abundance_long <- abundance |>
	pivot_longer(cols = everything(), names_to = "node", values_to = "value") |>
	mutate(n_id = as.numeric(stringr::str_extract(node, "(?<=\\[)\\d*(?=\\])")))

density_stats <- left_join(abundance_long, property_info) |>
	mutate(density = value / property_area_km2) |>
	group_by(
		node,
		n_id,
		property,
		primary_period,
		property_area_km2
	) |>
	summarise(
		mean = mean(density),
		variance = var(density),
		`0.025` = quantile(density, 0.025),
		`0.05` = quantile(density, 0.05),
		`0.1` = quantile(density, 0.1),
		`0.25` = quantile(density, 0.25),
		`0.5` = quantile(density, 0.5),
		`0.75` = quantile(density, 0.75),
		`0.9` = quantile(density, 0.9),
		`0.95` = quantile(density, 0.95),
		`0.975` = quantile(density, 0.975)
	) |>
	ungroup()

round(quantile(density_stats$`0.5`, c(0.05, 0.5, 0.95)), 1)

n_return <- data_for_nimble |>
	select(property, primary_period, method) |>
	distinct() |>
	pivot_wider(names_from = method, values_from = method) |>
	unite(method, -c(property, primary_period), sep = ", ", na.rm = TRUE) |>
	group_by(property, method) |>
	mutate(return_interval = c(0, diff(primary_period))) |>
	ungroup() |>
	rename(methods_used = method)

model_data <- data_for_nimble |>
	group_by(
		property,
		start_dates,
		end_dates,
		st_name,
		primary_period,
		property_area_km2,
		county_code
	) |>
	summarise(total_take = sum(take), total_effort_per = sum(effort_per)) |>
	ungroup() |>
	mutate(
		take_density = total_take / property_area_km2
	) |>
	left_join(density_stats) |>
	left_join(n_return)

all_slopes_data <- tibble()
all_props_data <- unique(model_data$property)

for (j in seq_along(all_props_data)) {
	tmp <- model_data |> filter(property == all_props_data[j])

	if (nrow(tmp) == 1) {
		next
	}

	m <- lm(`0.5` ~ primary_period, data = tmp)
	ypred <- predict(m, newdata = tibble(primary_period = tmp$primary_period))
	slope <- summary(m)$coefficients[2]

	sp <- tibble(
		property = all_props_data[j],
		lambda = slope
	)

	all_slopes_data <- bind_rows(all_slopes_data, sp)
}

thresh <- 0.05
directions_data <- model_data |>
	left_join(all_slopes_data) |>
	mutate(
		direction = if_else(lambda < (-1 * thresh), "Decreasing", "Stable"),
		direction = if_else(lambda > thresh, "Increasing", direction)
	)

directions_data |>
	select(property, direction) |>
	distinct() |>
	count(direction)

id1 <- 2 # Decreasing
id2 <- 21 # Stable
id3 <- 47 # Increasing

plot_data <- directions_data |>
	filter(property %in% c(id1, id2, id3))

method_shapes <- plot_data |>
	select(methods_used) |>
	mutate(
		methods_used = str_replace(methods_used, "Firearms", "Ground-shooting")
	) |>
	distinct() |>
	mutate(shape = seq_len(n()))

plot_data |>
	mutate(
		methods_used = as.factor(methods_used),
		methods_used = str_replace(methods_used, "Firearms", "Ground-shooting")
	) |>
	ggplot() +
	aes(x = end_dates) +
	geom_linerange(aes(
		ymin = `0.05`,
		ymax = `0.95`,
		linewidth = "90% CI"
	)) +
	geom_linerange(aes(
		ymin = `0.25`,
		ymax = `0.75`,
		linewidth = "50% CI"
	)) +
	geom_line(aes(y = `0.5`, linewidth = "Median"), linetype = "dashed") +
	geom_point(aes(y = take_density, shape = methods_used), size = 1.5) +
	scale_shape_manual(
		breaks = method_shapes$methods_used,
		values = method_shapes$shape,
		drop = FALSE
	) +
	scale_linewidth_manual(
		values = c(
			"90% CI" = 0.25,
			"50% CI" = 1,
			"Median" = 0.5
		)
	) +
	labs(
		y = expression("Density (pigs / " ~ km^2 ~ ")"),
		x = element_blank(),
		color = "Properties fit",
		shape = "Total pigs removed\nin primary period\nby method(s)",
		linewidth = "Statistic"
	) +
	scale_x_date(date_labels = "%b-%Y") +
	facet_wrap(~direction, scales = "free_x") +
	theme_bw() +
	theme(
		axis.text = element_text(size = 8),
		axis.title = element_text(size = 10),
		legend.title = element_text(size = 10),
		legend.text = element_text(size = 8),
		axis.text.x = element_text(angle = 90, vjust = 0.5)
	)

ggsave(
	file.path(out_path, "dataTimeSeries-v2.jpeg"),
	dpi = "retina",
	device = "jpeg",
	units = "cm",
	width = 16,
	height = 12
)
