fit_mcmc <- function(cl, modelCode, data, constants, n_iter, n_chains, custom_samplers = NULL, monitors_add = NULL){

  require(nimble)
  require(parallel)

  export <- c(
    "modelCode",
    "data",
    "constants",
    "n_iter",
    "monitors_add",
    "custom_samplers"
  )
  clusterExport(cl, export, envir = environment())


  source("R/inits.R")
  for(i in seq_along(cl)){
    set.seed(i)
    init_dir <- "/lustrefs/ceah/feral-swine/simulation/output/test/betaSurvival/density_5/1"
    init <- inits(data, constants)
    clusterExport(cl[i], "init", envir = environment())
  }

  message("Compiling model and initial parallel sampling...")
  start <- Sys.time()

  out <- clusterEvalQ(cl, {

    library(nimble)
    source("R/calc_log_potential_area.R")

    Rmodel <- nimbleModel(
      code = modelCode,
      constants = constants,
      data = data,
      inits = init,
      calculate = TRUE
    )

    # Rmodel$initializeInfo()

    for(i in 1:constants$n_survey){
      N_model <- Rmodel$N[constants$nH_p[i]]
      n <- round(N_model - constants$y_sum[i])
      if(n <= 0){
        print(i)
        n <- ifelse(n == 0, 2, n)
        Rmodel$N[constants$nH_p[i]] <- N_model + n^2
      }
    }

    # Rmodel$simulate()
    # Rmodel$calculate()

    # default MCMC configuration
    mcmcConf <- configureMCMC(Rmodel, useConjugacy = TRUE)

    if(!is.null(monitors_add)){
      mcmcConf$addMonitors(monitors_add)
    }

    if(!is.null(custom_samplers)){
      for(i in seq_len(nrow(custom_samplers))){
        node <- custom_samplers$node[i]
        type <- custom_samplers$type[i]
        mcmcConf$removeSampler(node)
        mcmcConf$addSampler(node, type)
      }
    }

    for(i in 1:5){
      node <- paste0("beta_p[", i, ", ", 1:constants$m_p, "]")
      node <- c(paste0("beta1[", i, "]"), node)
      mcmcConf$removeSampler(node)
      mcmcConf$addSampler(node, "AF_slice")
    }

    mcmcConf$printSamplers(byType = TRUE)

    Rmcmc <- buildMCMC(mcmcConf)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc)

    runMCMC(
      Cmcmc,
      niter = n_iter,
      nchains = 1,
      nburnin = round(n_iter / 2),
      samplesAsCodaMCMC = TRUE
    )

  })

  message("Run time for ", n_iter, " iterations across ", n_chains, " chains in parallel:")
  print(Sys.time() - start)

  return(out)

}
