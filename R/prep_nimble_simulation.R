
prep_nimble <- function(N, take, X){

  require(dplyr)
  require(tidyr)

  N_timestep <- N |>
    select(property, PPNum) |>
    group_by(property) |>
    mutate(timestep = 1:n()) |>
    ungroup()

  # the number of timesteps for each property, sampled or not
  n_time_prop <- N_timestep |>
    group_by(property) |>
    filter(timestep == max(timestep)) |>
    pull(timestep)

  # so we can iterate through each PP for each property and keep track of index
  all_pp_wide <- N_timestep |>
    pivot_wider(names_from = timestep,
                values_from = PPNum) |>
    select(-property)

  nH <- N_timestep |>
    mutate(n_id = 1:n()) |>
    select(-PPNum) |>
    pivot_wider(names_from = timestep,
                values_from = n_id) |>
    select(-property)

  # Generate start and end indices for previous surveys
  take_timestep <- take |>
    select(property, PPNum) |>
    distinct() |>
    group_by(property) |>
    mutate(timestep = 1:n()) |>
    ungroup()

  take <- left_join(take, take_timestep)

  take$start <- 0
  take$end <- 0

  pb <- txtProgressBar(max = nrow(take), style = 1)
  for (i in 1:nrow(take)) {
    if (take$order[i] > 1) {
      idx <- which(take$county == take$county[i] &
                     take$property == take$property[i] &
                     take$timestep == take$timestep[i] &
                     take$order < take$order[i])
      take$start[i] <- idx[1]
      take$end[i] <- idx[length(idx)]
      assertthat::are_equal(idx, take$start[i]:take$end[i])
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)


  y_sum <- take |>
    group_by(property, PPNum) |>
    mutate(ysum = cumsum(take) - take) |>
    ungroup() |>
    select(property, PPNum, ysum, order)

  sum_take <- take |>
    select(property, PPNum, sum_take) |>
    distinct()

  removed_timestep <- left_join(N_timestep, sum_take) |>
    mutate(sum_take = if_else(is.na(sum_take), 0, sum_take)) |>
    select(-PPNum) |>
    pivot_wider(names_from = timestep,
                values_from = sum_take) |>
    select(-property)

  tH <- take |>
    select(property, PPNum)

  nH_p <- N |>
    select(property, PPNum, n_id) |>
    right_join(tH) |>
    pull(n_id)

  # mean litters per year from VerCauteren et al. 2019 pg 64
  data_litters_per_year <- c(1, 2, 0.86, 1, 2.28, 2.9, 0.49, 0.85, 1.57)

  # mean litter size year from VerCauteren et al. 2019 pg 63
  data_litter_size <- c(5.6, 6.1, 5.6, 6.1, 4.2, 5.0, 5.0, 6.5, 5.5, 6.8,
                        5.6, 5.9, 4.9, 5.1, 4.5, 4.7, 5.3, 5.7, 7.4, 8.4,
                        4.7, 4.9, 3.0, 3.0, 4.8, 4.8, 4.2, 5.4, 4.7, 5.2, 5.4)

  data_litter_size <- round(data_litter_size)


  constants <- list(
    n_survey = nrow(take),
    # n_lpy = length(data_litters_per_year),
    n_ls = length(data_litter_size),
    n_property = max(N$property),
    n_first_survey = length(which(take$order == 1)),
    n_not_first_survey = length(which(take$order != 1)),
    n_method = 5,
    n_units = nrow(N),
    n_time_prop = n_time_prop,
    all_pp = as.matrix(all_pp_wide),
    nH = as.matrix(nH),
    nH_p = nH_p,
    log_pi = log(pi),
    m_p = ncol(X),
    first_survey = which(take$order == 1),
    not_first_survey = which(take$order != 1),
    p_property_idx = take$property,
    p_pp_idx = take$PPNum,
    start = take$start,
    end = take$end,
    method = as.numeric(as.factor(take$method)),
    county = take$county,
    property_x = N$property,
    pp_x = N$PPNum,
    pp_len = 28,
    phi_mu_a = 3.23,
    phi_mu_b = 0.2,
    y_sum = y_sum$ysum,
    rem = as.matrix(removed_timestep),
    log_rho_mu = rep(0, 5),
    log_rho_tau = rep(0.1, 5),
    p_mu_mu = rep(0, 2),
    p_mu_tau = rep(1, 2),
    log_gamma_mu = rep(0, 2),
    log_gamma_tau = rep(0.1, 2),
    beta1_mu = rep(0, 5),
    beta1_tau = rep(1, 5),
    beta_p_mu = rep(0, 15),
    beta_p_tau = rep(1, 15),
    psi_shape = 1,
    psi_rate = 0.1,
    log_nu_mu = 2,
    log_nu_tau = 1

  )

  data <- list(
    y = take$take,
    # K = data_litters_per_year,
    J = data_litter_size,
    X_p = X,
    effort_per = take$effort_per,
    log_effort_per = log(take$effort_per),
    n_trap_m1 = take$trap_count - 1,
    log_survey_area_km2 = log(take$property_area)
  )

  return(
    list(
      constants = constants,
      data = data
    )
  )
}








