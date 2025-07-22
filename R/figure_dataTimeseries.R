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

posterior_path <- file.path(top_dir, if_dir, "modelData.rds")
data <- read_rds(posterior_path) |>
	mutate(
		method = if_else(method == "FIREARMS", "Firearms", method),
		method = if_else(method == "FIXED WING", "Fixed wing", method),
		method = if_else(method == "HELICOPTER", "Helicopter", method),
		method = if_else(method == "SNARE", "Snare", method),
		method = if_else(method == "TRAPS", "Trap", method)
	)

posterior_path <- file.path(top_dir, if_dir, "densitySummaries.rds")
density <- read_rds(posterior_path)

n_return <- data |>
	select(propertyID, primary_period, method) |>
	distinct() |>
	pivot_wider(names_from = method, values_from = method) |>
	unite(method, -c(propertyID, primary_period), sep = ", ", na.rm = TRUE) |>
	group_by(propertyID, method) |>
	mutate(return_interval = c(0, diff(primary_period))) |>
	ungroup() |>
	rename(methods_used = method)

model_data <- data |>
	group_by(
		propertyID,
		agrp_prp_id,
		start_dates,
		end_dates,
		st_name,
		cnty_name,
		farm_bill,
		alws_agrprop_id,
		property,
		primary_period,
		property_area_km2,
		county_code
	) |>
	summarise(total_take = sum(take), total_effort_per = sum(effort_per)) |>
	ungroup() |>
	mutate(
		take_density = total_take / property_area_km2,
		farm_bill = if_else(is.na(farm_bill), 0, farm_bill)
	) |>
	left_join(density) |>
	left_join(n_return)

all_slopes_data <- tibble()
all_props_data <- unique(model_data$propertyID)

for (j in seq_along(all_props_data)) {
	tmp <- model_data |> filter(propertyID == all_props_data[j])

	if (nrow(tmp) == 1) {
		next
	}

	m <- lm(`0.5` ~ primary_period, data = tmp)
	ypred <- predict(m, newdata = tibble(primary_period = tmp$primary_period))
	slope <- summary(m)$coefficients[2]

	sp <- tibble(
		propertyID = all_props_data[j],
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
	select(propertyID, direction) |>
	distinct() |>
	count(direction)

id1 <- "352492-347989" # decreasing
id2 <- "359719-356143" # increasing
id3 <- "361285-357745" # stable

plot_data <- directions_data |>
	filter(propertyID %in% c(id1, id2, id3))

method_shapes <- plot_data |>
	select(methods_used) |>
	mutate(
		methods_used = str_replace(methods_used, "Firearms", "Sharpshooting")
	) |>
	distinct() |>
	mutate(shape = seq_len(n()))

plot_data |>
	mutate(
		methods_used = as.factor(methods_used),
		methods_used = str_replace(methods_used, "Firearms", "Sharpshooting")
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
	file.path(out_path, "dataTimeSeries"),
	dpi = "retina",
	device = "jpeg",
	units = "cm",
	width = 16,
	height = 12
)
