#' Run an ADMB model using the No-U-Turn Sampler (`NUTS`) or
#' Random Walk Metropolis (`RWM`) algorithms.
#'
#' @description
#' Draw Bayesian posterior samples from an AD Model Builder
#' (ADMB) model using an MCMC algorithm.
#'
#' @details
#' This function implements algorithm 6 of Hoffman and Gelman (2014),
#' and loosely follows package [rstan]. The step size can be
#' adapted or specified manually. The metric (i.e., mass matrix) can be
#' unit diagonal, adapted diagonal (default and recommended), a dense
#' matrix specified by the user, or an adapted dense matrix.
#' Further control of algorithms can be
#' specified with the `control` argument.
#' The adaptation scheme (step size and mass matrix) is based heavily on
#' those by the software `Stan`, and more details can be found in that
#' documentation and the vignette.
#'
#' For example, `sample_admb(algorithm = 'nuts')`
#' generates posterior samples from which inference can be made. Adaptation
#' schemes are used with NUTS so specifying tuning parameters is
#' not necessary, and parallel execution reduces overall run
#' time.
#'
#' The `RWM` algorithm provides no new functionality not available
#' from previous versions of `ADMB`. However, `sample_admb(algorithm = 'rwm')`
#' is setup for parallel execution. Note that the algorithms'
#' code lies in the ADMB source code, and [adnuts] provides a
#' wrapper for it.
#'
#' @param model Name of model (i.e., model.tpl)
#' @param path Path to model executable. Defaults to working
#' directory. Often best to have model files in a separate
#' subdirectory, particularly for parallel.
#' @param num_samples The number of samples to draw.
#' @param init A list of lists containing the initial parameter
#' vectors, one for each chain or a function. It is strongly
#' recommended to initialize multiple chains from dispersed
#' points. A of NULL signifies to use the starting values
#' present in the model (i.e., `obj$par`) for all chains.
#' @param num_chains The number of chains to run. If greater than 1,
#' use parallel execution (1 core per chain)
#' @param warmup The number of warm up iterations.
#' @param seeds A vector of seeds, one for each chain.
#' @param thin The thinning rate to apply to samples. Typically
#' not used with NUTS.
#' @param mceval Whether to run the model with `-mceval` on
#' samples from merged chains.
#' @param duration The number of minutes after which the model
#' will quit running.
#' @param control A list to control the sampler. See details for
#' further use.
#' @param algorithm The algorithm to use, one of `nuts` or `rwm`
#' @param skip_optimization Whether to run the optimizer before
#' running MCMC. This is rarely need as it is better to run it
#' once before to get the covariance matrix, or the estimates
#' are not needed with adaptive NUTS.
#' @param skip_monitor Whether to skip calculating diagnostics
#' (effective sample size, Rhat) via the `rstan::monitor`
#' function. This can be slow for models with high dimension or
#' many iterations. The result is used in plots and summaries
#' so it is recommended to turn on. If model run with
#' `skip_monitor = FALSE` you can recreate it post-hoc by
#' setting `fit$monitor=rstan::monitor(fit$samples,
#' fit$warmup, print=FALSE)`
#' @param skip_unbounded Whether to skip returning the unbounded
#' version of the posterior samples in addition to the bounded
#' ones. It may be advisable to set to FALSE for very large
#' models to save space
#' @param admb_args A character string which gets passed to the
#' command line, allowing finer control
#' @param adapt_delta The target acceptance rate
#' @param metric The mass metric to use. Options are: "unit" for a unit diagonal
#' matrix; `NULL` to estimate a diagonal matrix during warmup; a matrix
#' to be used directly (in untransformed space)
#' @param adapt_delta Whether adaptation of step size is turned on.
#' @param adapt_mass Whether adaptation of mass matrix is turned
#' on. Currently only allowed for diagonal metric
#' @param adapt_mass_dense Whether dense adaptation of mass
#' matrix is turned on
#' @param max_treedepth Maximum treedepth for the NUTS algorithm
#' @param stepsize The stepsize for the NUTS algorithm. If `NULL` it
#' will be adapted during warmup
#' @param adapt_init_buffer The initial buffer size during mass matrix
#' adaptation where sample information is not used (default 50)
#' @param adapt_term_buffer The terminal buffer size (default 75)
#' during mass matrix adpatation (final fast phase)
#' @param adapt_window The initial size of the mass matrix
#' adaptation window, which gets doubled each time thereafter
#' @param refresh The rate at which to refresh progress to the
#' console. Defaults to even 10%. A value of 0 turns off
#' progress updates
#' @param fn_logfile The name of a file to output `stdout` from `system()`
#' calls to. If `NULL`, `stdout` will be printed to the screen
#' @name samplers
NULL

#' Sampling from ADMB models
#'
#' @rdname samplers
#' @importFrom furrr future_map furrr_options
#' @importFrom future plan
#' @importFrom parallel detectCores
#' @export
sample_admb <- function(model,
                        path = NULL,
                        num_samples = 2000,
                        init = NULL,
                        num_chains = 3,
                        warmup = NULL,
                        seeds = NULL,
                        thin = 1,
                        mceval = FALSE,
                        duration = NULL,
                        control = NULL,
                        algorithm = c("nuts", "rwm"),
                        skip_optimization = TRUE,
                        skip_monitor = FALSE,
                        skip_unbounded = TRUE,
                        admb_args = NULL,
                        hess_step = FALSE,
                        fn_logfile = "model_output.log"){

  algorithm <- match.arg(algorithm)

  if(is.null(path)){
    stop("You must supply a path which contains the model files",
         call. = FALSE)
  }

  if(num_chains < 1){
    stop("The number of chains `num_chains` cannot be less than 1",
         call. = FALSE)
  }

  if(thin < 1){
    stop("The thinning value `thin` cannot be less than 1",
         call. = FALSE)
  }

  cores_avail  <- detectCores()
  if(num_chains > cores_avail) {
    num_chains <- cores_avail - 1
    warning(paste0("Specified number of chains greater than number ",
                   "of available cores, using number of available cores - 1 = ", num_chains))
  }
  if(!is.null(control) && !is.list(control)){
    stop("Control argument invalid, must be a list",
         call. = FALSE)
  }

  if(is.null(seeds)){
    seeds <- sample.int(1e7, size = num_chains)
    message("seeds is NULL, setting random seeds: ",
            paste(seeds, collapse = ", "))
  }
  if(num_samples <= 1 || !is.numeric(num_samples)){
    stop("`num_samples` must be > 1", call. = FALSE)
  }

  model <- get_model_executable(model, loc_dir = path)

  v <- check_admb_version(exe = model)

  if(as.numeric(v) <= 12.0 && !skip_unbounded){
    warning(paste0("Version ", v, " of ADMB is incompatible with ",
                   "`skip_unbounded` = `FALSE`, ignoring"))
    skip_unbounded <- TRUE
  }

  if(is.null(warmup)){
    warmup <- floor(num_samples / 2)
    message("`warmup` is `NULL`, setting to `floor(num_samples / 2)` = ", warmup)
  }

  if(algorithm == "nuts"){
    control <- update_control(control)
  }

  if(is.null(init)){
    message("Using MLE inits for each chain. Argument `init` was not ",
            "supplied.")
  }else if(is.function(init)){
    init <- map(seq_len(num_chains), ~{init()})
  }else if(!is.list(init)){
    stop("`init` must be `NULL`, a list, or a function",
         call. = FALSE)
  }
  if(!is.null(init) && length(init) != num_chains){
    stop("Length of `init` does not equal number of chains",
         call. = FALSE)
  }
  # Delete old *.psv, adaptation.csv and unbounded.csv files
  files_to_remove <- grep(".*\\.psv|adaptation\\.csv|unbounded\\.csv",
                          dir(path),
                          value = TRUE)
  if(length(files_to_remove)){
    unlink(file.path(path, files_to_remove))
  }

  map_func <- `if`(num_chains > 1, future_map, map)

  if(num_chains > 1){
    plan("multisession", workers = num_chains)
  }
  mcmc_out <- map_func(seq_len(num_chains), function(i)
    run_admb_sampler(chain_num = i,
                     path = path,
                     model = model,
                     duration = duration,
                     algorithm = algorithm,
                     num_samples = num_samples,
                     init = init[[i]],
                     warmup = warmup,
                     seed = seeds[i],
                     thin = thin,
                     control = control,
                     skip_optimization = skip_optimization,
                     admb_args = admb_args,
                     hess_step = hess_step,
                     fn_logfile = fn_logfile),
    .options = furrr_options(seed = TRUE))
  if(num_chains > 1){
    plan()
  }

  warmup <- mcmc_out[[1]]$warmup
  mle <- read_mle_fit(path)

  if(is.null(mle)){
    par_names <- dimnames(mcmc_out[[1]]$samples)[[2]]
    par_names <- par_names[-length(par_names)]
  }else{
    par_names <- mle$par_names
  }
  iters <- map_dbl(mcmc_out, ~{dim(.x$samples)[1]})
  if(any(iters != num_samples / thin)){
    chain_length <- rep(min(iters), length(iters))
    warning("Variable chain lengths, truncating to minimum = ", chain_length)
  }else{
    chain_length <- num_samples / thin
  }

  samples <- array(NA, dim = c(chain_length, num_chains, 1 + length(par_names)),
                   dimnames = list(NULL, NULL, c(par_names, "lp__")))

  if(skip_unbounded){
    samples_unbounded <- NULL
  }else{
    samples_unbounded <- samples
  }

  for(i in seq_len(num_chains)){
    samples[, i, ] <- mcmc_out[[i]]$samples[seq_len(chain_length), ]
    if(!skip_unbounded)
      samples_unbounded[, i, ] <-
        cbind(mcmc_out[[i]]$unbounded[seq_len(chain_length), ],
              mcmc_out[[i]]$samples[, 1 + length(par_names)])
  }
  sampler_params <- NULL
  if(algorithm == "NUTS"){
    sampler_params <- map(mcmc_out, ~{
      .x$sampler_params[seq_len(chain_length), ]
    })
  }
  cmd <- map_chr(mcmc_out, ~{.x$cmd})
  if(chain_length < warmup){
    warning("Duration too short to finish warmup period")
  }
  # When running multiple chains the psv files will either be overwritten
  # or in different folders if processed in parallel. Thus mceval needs to be
  # run after recombining chains
  message("Merging post-warmup chains into main folder: ", path)
  sa <- map(seq_len(num_chains), ~{
    samples[-seq_len(warmup), .x, -dim(samples)[3]]
  }) %>%
    do.call(rbind, .)
  write_psv(path = path,
            fn_psv = model,
            samples = sa)

  unbounded <- map(mcmc_out, ~{
    .x$unbounded
  }) %>%
    do.call(rbind, .)
  write.csv(unbounded,
            file.path(path, "unbounded.csv"),
            row.names = FALSE)

  if(mceval){
    message("Running -mceval on merged chains")
    mceval_cmd <- paste0("cd ", path, " && ", model, " -mceval")
    if(!is.null(fn_logfile)){
      mceval_cmd <- paste0(mceval_cmd, " > ", fn_logfile, " 2>&1")
    }
    system_(mceval_cmd, ignore.stdout = FALSE)
  }

  covar_est <- cov(unbounded)
  if(skip_monitor){
    message("Skipping ESS and Rhat statistics")
    mon <- NULL
  }else{
    message("Calculating ESS and Rhat (skip_monitor = TRUE will skip)")
    # The following line stops a duplicate rownames error from being thrown
    # in the rstan::monitor() function
    dimnames(samples)[[3]] <- NULL
    mon <- monitor(samples, warmup, print = FALSE)
  }
  result <- list(samples = samples,
                 sampler_params = sampler_params,
                 samples_unbounded = samples_unbounded,
                 algorithm = algorithm,
                 warmup = warmup,
                 model = model,
                 max_treedepth = mcmc_out[[1]]$max_treedepth,
                 cmd = cmd,
                 init = init,
                 covar_est = covar_est,
                 mle = mle,
                 monitor = mon)

  invisible(adfit(result))
}
