# prepare data for nimble, run a test build

prep_and_run_mcmc <- function(
	informed,
	post_path,
	data_repo,
	dest_mcmc,
	dest_posterior,
	df,
	monitors_add,
	custom_samplers
) {
	source("R/nimble_removal_model.R")
	source("R/prep_nimble_data.R")
	source("R/functions_nimble.R")
	source("R/mcmc_parallel.R")

	# ===================================================
	# Prepare data for NIMBLE ----
	# ===================================================

	constants <- nimble_constants(
		df = df,
		interval = 4,
		data_repo = data_repo,
		informed = informed
	)

	data <- nimble_data(df)

	n_chains <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
	if (is.na(n_chains)) {
		n_chains <- 3
	}

	message(n_chains, " chains")

	inits <- list()
	for (i in seq_len(n_chains)) {
		set.seed(i)
		inits[[i]] <- nimble_inits(constants, data)
	}

	test_build(modelCode, constants, data, inits = inits[[1]])

	# ===================================================
	# Run Nimble in parallel ----
	# ===================================================

	if (!dir.exists(dest_mcmc)) {
		dir.create(dest_mcmc, recursive = TRUE, showWarnings = FALSE)
	}
	if (!dir.exists(dest_posterior)) {
		dir.create(dest_posterior, recursive = TRUE, showWarnings = FALSE)
	}

	write_rds(df, file.path(dest_posterior, "modelData.rds"))

	message("Begin Parallel sampling and make cluster")
	cl <- makeCluster(n_chains)

	finished <- mcmc_parallel(
		cl = cl,
		model_code = modelCode,
		model_constants = constants,
		model_data = data,
		model_inits = inits,
		params_check = params_check,
		n_iters = config$n_iter,
		dest = dest_mcmc,
		monitors_add = monitors_add,
		custom_samplers = NULL
	)

	stopCluster(cl)
	return(finished)
}
