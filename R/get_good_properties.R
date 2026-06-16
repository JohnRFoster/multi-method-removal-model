library(dplyr)
library(tidyr)
library(readr)

source("R/functions_figures.R")

path <- "data"
density_dirs <- paste0("density_", c(0.3, 1.475, 2.65, 3.825, 5))

abundance_file <- "abundance_error_by_property.rds"
f_name <- abundance_file
df <- map_files2(density_dirs, f_name)

# create property ID for easier joining
property_ids <- df |>
	select(start_density, simulation, property) |>
	distinct() |>
	mutate(property_id = paste(start_density, simulation, property, sep = "-"))

df_property <- df |>
	left_join(property_ids) |>
	select(
		property_id,
		simulation,
		property,
		start_density,
		contains("density"),
		-density
	) |>
	distinct()

abundance_file <- "all_take.rds"
f_name <- abundance_file
df <- map_files2(density_dirs, f_name)

n_methods <- df |>
	left_join(property_ids) |>
	select(property_id, method) |>
	distinct() |>
	group_by(property_id) |>
	count() |>
	rename(n_methods_used = n)

n_return <- df |>
	left_join(property_ids) |>
	select(property_id, method) |>
	distinct() |>
	pivot_wider(names_from = method, values_from = method) |>
	unite(method, -c(property_id), sep = ", ", na.rm = TRUE) |>
	rename(methods_used = method) |>
	left_join(n_methods)

observed_pp <- df |>
	left_join(property_ids) |>
	select(property_id, PPNum) |>
	distinct() |>
	count(property_id) |>
	rename(n_observed_pp = n)

ts_length <- df |>
	left_join(property_ids) |>
	select(property_id, PPNum) |>
	distinct() |>
	group_by(property_id) |>
	filter(
		PPNum == min(PPNum) |
			PPNum == max(PPNum)
	) |>
	mutate(ts_length = c(0, diff(PPNum) + 1)) |>
	ungroup() |>
	filter(ts_length != 0) |>
	select(-PPNum)

take_by_property <- df |>
	left_join(property_ids) |>
	group_by(property_id, property_area, start_density) |>
	summarise(
		take = sum(take),
		effort = sum(effort_per),
		unit_count = sum(trap_count),
		n_total_events = n()
	) |>
	ungroup() |>
	left_join(n_return) |>
	left_join(observed_pp) |>
	left_join(ts_length) |>
	mutate(proportion_observed = n_observed_pp / ts_length)

property_error <- left_join(take_by_property, df_property)

thresh_mpe <- 25
thresh_nrmse <- 0.5
thresh_bias <- 0.5

acceptable_properties <- property_error |>
	filter(
		mpe_density <= thresh_mpe,
		nm_rmse_density <= thresh_nrmse,
		mbias_density <= thresh_bias,
		mbias_density >= -1 * thresh_bias
	)

nrow(acceptable_properties)
nrow(acceptable_properties) / nrow(property_error) * 100

create_range_df <- function(df, x) {
	tibble(
		metric = x,
		qmin = df |> pull(x) |> quantile(0.05),
		qmax = df |> pull(x) |> quantile(0.95)
	)
}

acceptable_metrics <- bind_rows(
	create_range_df(acceptable_properties, "property_area"),
	create_range_df(acceptable_properties, "take"),
	create_range_df(acceptable_properties, "n_total_events"),
	create_range_df(acceptable_properties, "n_observed_pp"),
	create_range_df(acceptable_properties, "ts_length"),
	create_range_df(acceptable_properties, "proportion_observed"),
	create_range_df(acceptable_properties, "effort"),
	create_range_df(acceptable_properties, "unit_count")
)

acceptable_metrics

write_csv(acceptable_metrics, file.path(path, "good_property_metrics.csv"))
