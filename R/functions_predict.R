#'@description Calculate the potential search area from posterior samples when traps or snares are used
#'@param log_rho vector of mcmc samples for rho (log scale)
#'@param log_gamma vector of mcmc samples for gamma (log scale)
#'@param p_unique vector of mcmc samples for p
#'@param effort_per the effort per trap/snare
#'@param n_trap_m1 the number of traps/snares used minus 1

trap_snare_lpa <- function(
  log_rho,
  log_gamma,
  p_unique,
  effort_per,
  n_trap_m1
) {
  log(pi) +
    (2 * (log_rho + log(effort_per) - log(exp(log_gamma) + effort_per))) +
    log(1 + (p_unique * n_trap_m1))
}

#'@description Calculate the potential search area from posterior samples when aerial methods (fixed wing, helicopter) are used
#'@param log_rho vector of mcmc samples for rho (log scale)
#'@param effort_per the effort per trap/snare

aerial_lpa <- function(log_rho, effort_per) {
  log_rho + log(effort_per)
}

#'@description Calculate the potential search area from posterior samples when firearms are used
#'@param log_rho vector of mcmc samples for rho (log scale)
#'@param effort_per the effort per trap/snare

firearms_lpa <- function(log_rho, effort_per) {
  log_rho + log(effort_per)
}

#'@description Calculate the potential search area from posterior samples for a vector of methods (used in forecasting code)
#'@param method method used for removal
#'@param method_params data frame of parameters specific to each method
#'@param effort_per the effort per unit
#'@param n_trap_m1 the number of units used minus 1

lpa <- function(method, method_params, effort_per, n_trap_m1) {
  n_ens <- nrow(method_params)
  log_potential_area <- numeric(n_ens)
  for (e in seq_len(n_ens)) {
    if (method == "TRAPS" || method == "SNARES") {
      log_potential_area[e] <- trap_snare_lpa(
        log_rho = method_params$log_rho[e],
        log_gamma = method_params$log_gamma[e],
        p_unique = nimble::ilogit(method_params$p_mu[e]),
        effort_per = effort_per,
        n_trap_m1 = n_trap_m1
      )
    } else if (method == "HELICOPTER" || method == "FIXED_WING") {
      log_potential_area[e] <- aerial_lpa(
        log_rho = method_params$log_rho[e],
        effort_per = effort_per
      )
    } else if (method == "FIREARMS") {
      log_potential_area[e] <- firearms_lpa(
        log_rho = method_params$log_rho[e],
        effort_per = effort_per
      )
    }
  }
  log_potential_area
}

#'@description calculate the capture probability for a removal event
#'@param X matrix of county-level land cover covariates, first column is 1's for an intercept term
#'@param beta vector of coeffiencents
#'@param log_potential_area matrix of potential search area, iterations (rows) by replicate (columns)
#'@param area_property scalar value of the area (mk^2) of the property

calc_p <- function(X, beta1, beta, log_potential_area, area_property, method) {
  n_mcmc <- nrow(beta)
  n_reps <- ncol(log_potential_area)
  p <- matrix(NA, n_mcmc, n_reps)
  for (mc in seq_len(n_mcmc)) {
    log_theta <- numeric(n_reps)
    for (j in seq_len(n_reps)) {
      # base probability of capture given an individual is the first survey
      log_theta[j] <- log(ilogit(beta1[mc] + inprod(X, beta[mc, ]))) +
        pmin(0, log_potential_area[mc, j] - log(area_property))

      # the probability an individual is captured on the first survey
      if (j == 1) {
        p[mc, j] <- exp(log_theta[j])
      } else {
        # the probability an individual is captured after the first survey
        for (k in 2:n_reps) {
          p[mc, j] <- exp(
            log_theta[1] +
              sum(log(1 - exp(log_theta[1:(j - 1)])))
          )
        }
      }
    }
  }
  return(p)
}


data_posteriors <- function(samples, constants, data) {
  require(dplyr)
  require(nimble)
  require(stringr)
  samples <- as.matrix(samples)
  D <- append(constants, data)

  post <- with(D, {
    log_potential_area <- matrix(NA, nrow(samples), n_survey)
    y_pred <- matrix(NA, nrow(samples), n_survey)
    y_lambda <- matrix(NA, nrow(samples), n_survey)

    pattern <- "(?<!beta_)p\\[\\d*\\]"
    p_detect <- str_detect(colnames(samples), pattern)
    if (any(p_detect)) {
      calc_p <- FALSE
      p <- samples[, which(p_detect)]
      log_theta <- NA
    } else {
      calc_p <- TRUE
      p <- matrix(NA, nrow(samples), n_survey)
      log_theta <- matrix(NA, nrow(samples), n_survey)
    }

    log_rho <- as.matrix(samples[, grep("log_rho", colnames(samples))])
    log_gamma <- as.matrix(samples[, grep("log_gamma", colnames(samples))])
    p_unique <- as.matrix(ilogit(samples[, grep("p_mu", colnames(samples))]))
    beta1 <- as.matrix(samples[, grep("beta1", colnames(samples))])
    beta_p <- as.matrix(samples[, grep("beta_p", colnames(samples))])
    n_ens <- nrow(samples)

    for (i in 1:n_survey) {
      if (method[i] == 1) {
        log_potential_area[, i] <- firearms_lpa(
          log_rho = log_rho[, method[i]],
          effort_per = effort_per[i]
        )
      } else if (method[i] == 2 | method[i] == 3) {
        log_potential_area[, i] <- aerial_lpa(
          log_rho = log_rho[, method[i]],
          effort_per = effort_per[i]
        )
      } else if (method[i] == 4 | method[i] == 5) {
        log_potential_area[, i] <- trap_snare_lpa(
          log_rho = log_rho[, method[i]],
          log_gamma = log_gamma[, method[i] - 3],
          p_unique = p_unique[, method[i] - 3],
          effort_per = effort_per[i],
          n_trap_m1 = n_trap_m1[i]
        )
      }
    }

    if (calc_p) {
      for (i in 1:n_survey) {
        M <- method[i]

        beta_p_nodes <- paste0("beta_p[", M, ", ", 1:m_p, "]")

        log_theta[, i] <- log(
          ilogit(
            beta1[, M] +
              X_p[county[i], 1] * beta_p[, beta_p_nodes[1]] +
              X_p[county[i], 2] * beta_p[, beta_p_nodes[2]] +
              X_p[county[i], 3] * beta_p[, beta_p_nodes[3]]
          )
        ) +
          pmin(0, log_potential_area[, i] - log_survey_area_km2[i])
      }

      # the probability an individual is captured on the first survey
      p[, first_survey] <- exp(log_theta[, first_survey])

      # the probability an individual is captured after the first survey
      for (e in 1:n_ens) {
        for (i in 1:n_not_first_survey) {
          p[e, not_first_survey[i]] <- exp(
            log_theta[e, start[not_first_survey[i]]] +
              sum(log(
                1 -
                  exp(log_theta[
                    e,
                    start[not_first_survey[i]]:end[not_first_survey[i]]
                  ])
              ))
          )
        }
      }
    }

    for (i in 1:n_survey) {
      N <- samples[, paste0("N[", nH_p[i], "]")]
      y_lambda[, i] <- (N - y_sum[i]) * p[, i]
      y_pred[, i] <- rpois(length(N), y_lambda[, i])
    }

    list(
      y_pred = y_pred,
      y_lambda = y_lambda,
      p_pred = p,
      potential_area_pred = exp(log_potential_area),
      theta_pred = exp(log_theta)
    )
  })
  return(post)
}
