#' Bayesian inference of an ADMB model using the no-U-turn
#' sampler (NUTS) or random walk Metropolis (RWM) algorithms.
#'
#' @description Draw Bayesian posterior samples from an AD Model Builder
#' (ADMB) model using an MCMC algorithm. `sample_admb(algorithm = 'NUTS')` generates
#' posterior samples from which inference can be made. Adaptation
#' schemes are used with NUTS so specifying tuning parameters is
#' not necessary, and parallel execution reduces overall run
#' time.
#'
#' The RWM algorithm provides no new functionality not available
#' from previous versions of ADMB. However, `sample_admb(algorithm = 'RWM')` has an
#' improved console output, is setup for parallel execution, and
#' a smooth workflow for diagnostics. Note that the algorithms'
#' code lies in the ADMB source code, and [adnuts] provides a
#' wrapper for it. See vignette for more information.
#'
#' @details This function implements algorithm 6 of Hoffman and Gelman (2014),
#' and loosely follows package [rstan]. The step size can be
#' adapted or specified manually. The metric (i.e., mass matrix) can be
#' unit diagonal, adapted diagonal (default and recommended), a dense
#' matrix specified by the user, or an adapted dense matrix.
#' Further control of algorithms can be
#' specified with the `control` argument.
#' The adaptation scheme (step size and mass matrix) is based heavily on those by the
#' software Stan, and more details can be found in that documentation and the vignette
#'
#' @author Cole Monnahan
#' @param model Name of model (i.e., model.tpl)
#' @param path Path to model executable. Defaults to working
#'   directory. Often best to have model files in a separate
#'   subdirectory, particularly for parallel.
#' @param iter The number of samples to draw.
#' @param init A list of lists containing the initial parameter
#'   vectors, one for each chain or a function. It is strongly
#'   recommended to initialize multiple chains from dispersed
#'   points. A of NULL signifies to use the starting values
#'   present in the model (i.e., `obj$par`) for all chains.
#' @param chains The number of chains to run
#' @param parallel If `TRUE`, run all chains in parallel, 1 CPU for
#' each chain. If `FALSE` run all chains sequentially on one CPU
#' @param warmup The number of warm up iterations.
#' @param seeds A vector of seeds, one for each chain.
#' @param thin The thinning rate to apply to samples. Typically
#'   not used with NUTS.
#' @param mceval Whether to run the model with `-mceval` on
#'   samples from merged chains.
#' @param duration The number of minutes after which the model
#'   will quit running.
#' @param control A list to control the sampler. See details for
#'   further use.
#' @param skip_optimization Whether to run the optimizer before
#'   running MCMC. This is rarely need as it is better to run it
#'   once before to get the covariance matrix, or the estimates
#'   are not needed with adaptive NUTS.
#' @param skip_monitor Whether to skip calculating diagnostics
#'   (effective sample size, Rhat) via the `rstan::monitor`
#'   function. This can be slow for models with high dimension or
#'   many iterations. The result is used in plots and summaries
#'   so it is recommended to turn on. If model run with
#'   `skip_monitor = FALSE` you can recreate it post-hoc by
#'   setting `fit$monitor=rstan::monitor(fit$samples,
#'   fit$warmup, print=FALSE)`
#' @param skip_unbounded Whether to skip returning the unbounded
#'   version of the posterior samples in addition to the bounded
#'   ones. It may be advisable to set to FALSE for very large
#'   models to save space
#' @param admb_args A character string which gets passed to the
#'   command line, allowing finer control
#' @param adapt_delta The target acceptance rate
#' @param metric The mass metric to use. Options are: "unit" for a unit diagonal
#'   matrix; `NULL` to estimate a diagonal matrix during warmup; a matrix
#'   to be used directly (in untransformed space)
#' @param adapt_delta Whether adaptation of step size is turned on.
#' @param adapt_mass Whether adaptation of mass matrix is turned
#'   on. Currently only allowed for diagonal metric
#' @param adapt_mass_dense Whether dense adaptation of mass
#'   matrix is turned on
#' @param max_treedepth Maximum treedepth for the NUTS algorithm
#' @param stepsize The stepsize for the NUTS algorithm. If `NULL` it
#'   will be adapted during warmup
#' @param adapt_init_buffer The initial buffer size during mass matrix
#'   adaptation where sample information is not used (default
#'   50)
#' @param adapt_term_buffer The terminal buffer size (default 75)
#'   during mass matrix adpatation (final fast phase)
#' @param adapt_window The initial size of the mass matrix
#'   adaptation window, which gets doubled each time thereafter
#' @param refresh The rate at which to refresh progress to the
#'   console. Defaults to even 10%. A value of 0 turns off
#'   progress updates
#' @section Warning: The user is responsible for specifying the
#'   model properly (priors, starting values, desired parameters
#'   fixed, etc.), as well as assessing the convergence and
#'   validity of the resulting samples (e.g., through the
#'   [coda] package), or with function
#'   [launch_shinytmb()]` before making inference.
#'   Specifically, priors must be specified in the
#'   template file for each parameter. Unspecified priors will be
#'   implicitly uniform
#' @examples
#' \dontrun{
#' ## This is the packaged simple regression model
#' path.simple <- system.file('examples', 'simple', package='adnuts')
#' ## It is best to have your ADMB files in a separate folder and provide that
#' ## path, so make a copy of the model folder locally.
#' path <- 'simple'
#' dir.create(path)
#' trash <- file.copy(from=list.files(path.simple, full.names=TRUE), to=path)
#' ## Compile and run model
#' oldwd <- getwd()
#' setwd(path)
#' system_('admb simple.tpl')
#' system_('simple')
#' setwd('..')
#' init <- function() rnorm(2)
#' ## Run NUTS with defaults
#' fit <- sample_nuts(model='simple', init=init, path=path)
#' unlink(path, TRUE) # cleanup folder
#' setwd(oldwd)
#' }
#'
#' @name wrappers
NULL

#' Sampling from ADMB models
#'
#' @rdname wrappers
#' @param algorithm The algorithm to use, one of "NUTS" or "RWM"
#' @export
#' @importFrom furrr future_map
#' @importFrom future plan
sample_admb <- function(model,
                        path = NULL,
                        iter = 2000,
                        init = NULL,
                        chains = 3,
                        warmup = NULL,
                        seeds = NULL,
                        thin = 1,
                        mceval = FALSE,
                        duration = NULL,
                        control = NULL,
                        algorithm = "NUTS",
                        skip_optimization = TRUE,
                        skip_monitor = FALSE,
                        skip_unbounded = TRUE,
                        admb_args = NULL,
                        hess_step = FALSE,
                        parallel = TRUE){

  if(is.null(path)){
    stop("You must supply a path which contains the model files", call. = FALSE)
  }
  cores_avail  <- parallel::detectCores()
  if(parallel && chains > cores_avail) {
    chains <- cores_avail - 1
    warning(paste0("Specified number of chains greater than number ",
                   "of available cores, using number of available cores - 1 = ", chains))
  }
  if(!is.null(control) && !is.list(control)){
    stop("Control argument invalid, must be a list", call. = FALSE)
  }
  stopifnot(chains >= 1)
  stopifnot(thin >= 1)
  if(is.null(seeds)){
    seeds <- sample.int(1e7, size = chains)
    message("seeds is NULL, setting random seeds: ", paste(seeds, collapse = ","))
  }
  if(iter < 1 || !is.numeric(iter)){
    stop("iter must be > 1", call. = FALSE)
  }

  model <- get_model_executable(model, path)

  v <- .check_ADMB_version(model = model,
                           path = path,
                           warn = (algorithm == "NUTS"))
  if(v <= 12.0 && !skip_unbounded){
    warning(paste0("Version ", v, " of ADMB is incompatible with skip_unbounded = FALSE, ignoring"))
    skip_unbounded <- TRUE
  }
  # Update control with defaults
  if(is.null(warmup)){
    warmup <- floor(iter / 2)
    message("warning is NULL, setting to floor(iter / 2) = ", warmup)
  }
  if(!(algorithm %in% c("NUTS", "RWM"))){
    stop("Invalid algorithm specified",
         call. = FALSE)
  }
  if(algorithm == "NUTS")
    control <- .update_control(control)
  if(is.null(init)){
    warning("Using MLE inits for each chain -- strongly recommended to use dispersed inits")
  }else if(is.function(init)){
    init <- lapply(1:chains, function(x) init())
  }else if(!is.list(init)){
    stop("init must be NULL, a list, or a function",
         call. = FALSE)
  }
  if(!is.null(init) && length(init) != chains){
    stop("Length of init does not equal number of chains",
         call. = FALSE)
  }
  # Delete any psv files, adaptation.csv, and unbounded.csv in case something goes wrong we don't use old
  # values by accident
  files_to_trash <- grep("\\.psv|adaptation\\.csv|unbounded\\.csv", dir(path), value = TRUE)
  if(length(files_to_trash)){
    unlink(file.path(path, files_to_trash))
  }

  if(parallel){
    #plan("multisession", workers = chains)
    snowfall::sfStop()
    snowfall::sfInit(parallel=TRUE, cpus=chains)
    if(length(ls(envir=globalenv()))>0)
      snowfall::sfExportAll()
    on.exit(snowfall::sfStop())
    #mcmc.out <- future_map(seq_len(chains), function(i)
    mcmc.out <- snowfall::sfApply(seq_len(chains), function(i)
        sample_admb_parallel(parallel_number = i,
                           path = path,
                           model = model,
                           duration=duration,
                           algorithm = algorithm,
                           iter = iter,
                           init = init[[i]],
                           warmup = warmup,
                           seed = seeds[i],
                           thin = thin,
                           control = control,
                           skip_optimization = skip_optimization,
                           admb_args = admb_args,
                           hess_step = hess_step))
    plan()
  }else{
    mcmc.out <- lapply(seq_len(chains), function(i){
      sample_admb_parallel(parallel_number = i,
                           path = path,
                           model = model,
                           duration = duration,
                           algorithm = algorithm,
                           iter = iter,
                           init = init[[i]],
                           warmup = warmup,
                           seed = seeds[i],
                           thin = thin,
                           control = control,
                           skip_optimization = skip_optimization,
                           admb_args = admb_args,
                           hess_step = hess_step)})
  }
  warmup <- mcmc.out[[1]]$warmup
  mle <- read_mle_fit(model = model, path = path)

  if(is.null(mle)){
    par.names <- dimnames(mcmc.out[[1]]$samples)[[2]]
    par.names <- par.names[-length(par.names)]
  }else{
    par.names <- mle$par.names
  }
  iters <- unlist(lapply(mcmc.out, function(x) dim(x$samples)[1]))
  if(any(iters != iter / thin)){
    N <- min(iters)
    warning(paste0("Variable chain lengths, truncating to minimum = ", N))
  }else{
    N <- iter / thin
  }
  samples <- array(NA, dim = c(N, chains, 1 + length(par.names)),
                   dimnames = list(NULL, NULL, c(par.names, "lp__")))

  if(skip_unbounded){
    samples_unbounded <- NULL
  }else{
    samples_unbounded <- samples
  }
  for(i in 1:chains){
    samples[, i, ] <- mcmc.out[[i]]$samples[1:N, ]
    if(!skip_unbounded)
      samples_unbounded[, i, ] <- cbind(mcmc.out[[i]]$unbounded[1:N, ],
                                        mcmc.out[[i]]$samples[, 1 + length(par.names)])
  }
  if(algorithm == "NUTS"){
    sampler_params <- lapply(mcmc.out, function(x) x$sampler_params[1:N, ])
  }else{
    sampler_params <- NULL
  }
  time_warmup <- unlist(lapply(mcmc.out, function(x) as.numeric(x$time.warmup)))
  time_total <- unlist(lapply(mcmc.out, function(x) as.numeric(x$time.total)))
  cmd <- unlist(lapply(mcmc.out, function(x) x$cmd))
  if(N < warmup){
    warning("Duration too short to finish warmup period")
  }
  # When running multiple chains the psv files will either be overwritten
  # or in different folders (if parallel is used). Thus mceval needs to be
  # done posthoc by recombining chains AFTER thinning and warmup and
  # discarded into a single chain, written to file, then call -mceval.
  # Merge all chains together and run mceval
  message("Merging post-warmup chains into main folder: ", path)
  samples2 <- do.call(rbind, lapply(1:chains, function(i)
    samples[-(1:warmup), i, -dim(samples)[3]]))
  write_psv(fn = model, samples = samples2, path = path)
  # These already exclude warmup
  unbounded <- do.call(rbind, lapply(mcmc.out, function(x) x$unbounded))
  write.table(unbounded, file = file.path(path, "unbounded.csv"), sep = ",", col.names = FALSE, row.names = FALSE)
  if(mceval){
    message("Running -mceval on merged chains")
    system_(paste0("cd ", path, " && ", model, " -mceval"), ignore.stdout = FALSE)
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
                 time_warmup = time_warmup,
                 time_total = time_total,
                 algorithm = algorithm,
                 warmup = warmup,
                 model = model,
                 max_treedepth = mcmc.out[[1]]$max_treedepth,
                 cmd = cmd,
                 init = init,
                 covar_est = covar_est,
                 mle = mle,
                 monitor = mon)

  invisible(adfit(result))
}

