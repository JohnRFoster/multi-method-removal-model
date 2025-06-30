# code for manuscript figures

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(stringr)


# figure 2 simulation time series ----

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

# figure 3 heat maps ----

single <- c(
	"Firearms",
	"Fixed wing",
	"Helicopter",
	"Snares",
	"Traps"
)

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


# figure 4 bias and rmsle ----

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


# figure 5 real time series

method_shapes <- plot_data |>
	select(methods_used) |>
	distinct() |>
	mutate(shape = 1:n())

plot_data |>
	mutate(methods_used = as.factor(methods_used)) |>
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

# figure 6 potential area

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
