

inits <- function(data, constants, dir = NULL){
  source("R/calc_log_potential_area.R")
  with(append(data, constants), {

    if(is.null(dir)){
      beta1 <- rnorm(n_method, 0, 0.25)
      beta_p <- matrix(rnorm(m_p*n_method, 0, 0.1), n_method, m_p)
      p_mu <- rnorm(2)
      log_gamma <- log(runif(2, 0.1, 2))
      log_rho <- log(
        c(runif(1, 0.1, 5), runif(1, 50, 150), runif(1, 50, 150), runif(1, 5, 15), runif(1, 5, 15))
      )
      psi_phi <- runif(1, 2, 4)
      phi_mu <- runif(1, 0.7, 0.8)
      mean_ls <- round(runif(1, 5, 8))
    } else {
      burn <- "parameters_burnin.rds"
      params <- read_rds(file.path(dir, burn))
      params <- as.matrix(params[[1]])
      init_mu <- apply(params, 2, mean)
      beta_p <- matrix(jitter(init_mu[grep("beta_p[", names(init_mu), fixed = TRUE)]), constants$n_method, constants$m_p)
      beta1 <- jitter(init_mu[grep("beta1[", names(init_mu), fixed = TRUE)])
      p_mu <- jitter(init_mu[grep("p_mu[", names(init_mu), fixed = TRUE)])
      phi_mu <- jitter(init_mu[grep("phi_mu", names(init_mu), fixed = TRUE)])
      psi_phi <- jitter(init_mu[grep("psi_phi", names(init_mu), fixed = TRUE)])
      log_nu <- jitter(init_mu[grep("log_nu", names(init_mu), fixed = TRUE)])
      mean_ls <- exp(log_nu)
      log_gamma <- jitter(init_mu[grep("log_gamma[", names(init_mu), fixed = TRUE)])
      log_rho <- jitter(init_mu[grep("log_rho[", names(init_mu), fixed = TRUE)])
    }

    a <- phi_mu * psi_phi
    b <- (1 - phi_mu) * psi_phi
    mean_lpy <- 1
    zeta <- mean_lpy / 365 * pp_len * mean_ls
    N <- phi <- rep(NA, n_units)
    n_init <- rep(NA, n_property)
    for(i in 1:n_property){
      n_init[i] <- round(exp(log_survey_area_km2[i]) * 5) + sum(rem[i, ], na.rm = TRUE) * 2
      N[nH[i, 1]] <- n_init[i]
      for(j in 2:n_time_prop[i]){
        phi[nH[i, j-1]] <- rbeta(1, a, b)
        z <- N[nH[i, j-1]] - rem[i, j-1]
        z <- max(2, z)
        lambda <- z * zeta / 2 + z * phi[nH[i, j-1]]

        N[nH[i, j]] <- rpois(1, lambda)
      }
    }

    buffer <- 250
    list(
      log_lambda_1 = log(n_init + buffer),
      beta_p = beta_p,
      beta1 = beta1,
      p_mu = p_mu,
      phi_mu = phi_mu,
      psi_phi = psi_phi,
      N = N + buffer,
      log_nu = log(mean_ls),
      log_gamma = log_gamma,
      log_rho = log_rho,
      phi = phi
    )
  })

}
