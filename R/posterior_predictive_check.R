library(dplyr)
library(tidyr)
library(readr)
library(nimble)
library(ggplot2)
library(ggpubr)

set.seed(1256)

config_name <- "default"
config <- config::get(config = config_name)
interval <- config$interval

source("R/prep_nimble_data.R")
source("R/nimble_removal_model.R")
source("R/functions_nimble.R")

data_for_nimble <- read_csv("../data-store/masked_mis_data.csv") |>
  mutate(
    property = propertyID,
    county = county_code,
    method_vec = as.numeric(as.factor(method))
  )

posterior_parameters <- read_rds("data/1_posterior/posteriorSamples.rds")
posterior_states <- read_rds("data/1_posterior/stateSamples.rds")

posterior_samples <- bind_cols(
  posterior_parameters,
  posterior_states
)


constants <- nimble_constants(
  df = data_for_nimble,
  interval = interval
)

data <- nimble_data(data_for_nimble)

attach(constants)
attach(data)

n_samp <- nrow(posterior_samples)
log_potential_area <- log_theta <- p <- lambda <- y_pred <- matrix(NA, n_samp, n_survey)
n_nodes <- grep("N[", colnames(posterior_samples), fixed = TRUE)


pb <- txtProgressBar(max = n_samp, style = 3)
for(i in seq_len(n_samp)){

  jp <- posterior_samples[i, ] # joint posterior
  log_rho <- as.numeric(jp[paste0("log_rho[", 1:n_method, "]")])
  log_gamma <- as.numeric(jp[paste0("log_gamma[", 1:2, "]")])
  p_mu <- as.numeric(jp[paste0("p_mu[", 1:2, "]")])
  p_unique <- boot::inv.logit(p_mu)
  beta1 <- as.numeric(jp[paste0("beta1[", 1:n_method, "]")])

  for(j in seq_len(n_survey)){
    log_potential_area[i, j] <- calc_log_potential_area(
      log_rho = log_rho,
      log_gamma = log_gamma,
      p_unique = p_unique,
      log_effort_per = log_effort_per[j],
      effort_per = effort_per[j],
      n_trap_m1 = n_trap_m1[j],
      log_pi = log_pi,
      method = method[j]
    )

    beta_p <- as.numeric(jp[paste0("beta_p[", method[j], ", ", 1:m_p, "]")])

    # probability of capture, given that an individual is in the surveyed area
    log_theta[i, j] <- log(
      ilogit(
        beta1[method[j]] +
          inprod(X_p[j, ], beta_p)
      )
    ) +
      min(0, log_potential_area[i, j] - log_survey_area_km2[j])
  }

  # the probability an individual is captured on the first survey
  for (j in seq_len(n_first_survey)) {
    p[i, first_survey[j]] <- exp(log_theta[i, first_survey[j]])
  }

  # the probability an individual is captured after the first survey
  for (j in seq_len(n_not_first_survey)) {
    p[i, not_first_survey[j]] <- exp(
      log_theta[i, not_first_survey[j]] +
        sum(log(
          1 - exp(log_theta[i, start[not_first_survey[j]]:end[not_first_survey[j]]])
        )))
  }


  N <- as.numeric(jp[paste0("N[", nH_p, "]")])
  lambda[i, ] <- p[i, ] * (N - y_sum)
  y_pred[i, ] <- rpois(length(N), lambda[i, ])
  setTxtProgressBar(pb, i)
}
close(pb)

colnames(y_pred) <- paste0("y[", seq_len(nrow(data_for_nimble)), "]")

pps <- y_pred |>
  as_tibble() |>
  mutate(sim = seq_len(n())) |>
  pivot_longer(cols = -sim, names_to = "node", values_to = "value")

info <- data_for_nimble |>
  select(propertyID, primary_period, method, take) |>
  mutate(node = paste0("y[", seq_len(nrow(data_for_nimble)), "]"))

pps_info <- pps |>
  left_join(info, by = "node")

pp_by_property <- pps_info |>
  group_by(propertyID, sim) |>
  summarise(predicted_take = sum(value), .groups = "drop")

obs_by_property <- data_for_nimble |>
  group_by(propertyID) |>
  summarise(observed_take = sum(take), .groups = "drop")

get_p_value <- function(df_obs, df_pred, func) {
  take_by_property <- df_obs |>
    filter(observed_take == round(func(observed_take)))

  properties <- take_by_property$propertyID
  value <- unique(take_by_property$observed_take)

  # posterior predictive distribution of the number of PP with zero take
  prop_take_sim <- df_pred |>
    filter(propertyID %in% properties) |>
    mutate(above = ifelse(predicted_take >= value, 1, 0))

  p_val <- sum(prop_take_sim$above) / nrow(prop_take_sim)
  list(simulations = prop_take_sim, p_value = p_val, func_value = value)
}

p_min <- get_p_value(
  obs_by_property,
  pp_by_property,
  min
)

p_max <- get_p_value(
  obs_by_property,
  pp_by_property,
  max
)

p_med <- get_p_value(
  obs_by_property,
  pp_by_property,
  median
)


take_by_property_sd <- sd(obs_by_property$observed_take)

# posterior predictive distribution of the number of PP with zero take
prop_take_sim <- pp_by_property |>
  group_by(sim) |>
  summarise(sd_pred = sd(predicted_take)) |>
  mutate(above = ifelse(sd_pred >= take_by_property_sd, 1, 0))

p_val <- sum(prop_take_sim$above) / nrow(prop_take_sim)
p_sd <- list(
  simulations = prop_take_sim,
  p_value = p_val,
  func_value = take_by_property_sd
)

my_hist <- function(ls, f) {
  if (f == "min") {
    x_lab <- "T(y) = Property min(y)"
  } else if (f == "max") {
    x_lab <- "T(y) = Property max(y)"
  } else if (f == "median") {
    x_lab <- "T(y) = Property median(y)"
  } else if (f == "sd") {
    x_lab <- "T(y) = sd(y)"
    ls$simulations <- ls$simulations |>
      rename(predicted_take = sd_pred)
  }

  with(ls, {
    simulations |>
      ggplot() +
      aes(x = predicted_take) +
      geom_histogram() +
      geom_vline(xintercept = func_value, color = "red") +
      labs(
        title = paste0("p-value = ", round(p_value, 3)),
        x = x_lab,
        y = "Frequency"
      ) +
      theme_bw()
  })
}

g <- list()
g[[1]] <- my_hist(p_min, "min")
g[[2]] <- my_hist(p_max, "max")
g[[3]] <- my_hist(p_med, "median")
g[[4]] <- my_hist(p_sd, "sd")

ggarrange(
  plotlist = g,
  ncol = 2,
  nrow = 2,
  labels = "AUTO"
)

ggsave(
  "plots/posterior_predictive_checks.jpg",
  width = 6,
  height = 6,
  units = "in",
  dpi = "retina"
)

