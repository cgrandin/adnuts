#' Run a single NUTS chain for an ADMB model
#'
#' A low level function to run a single chain. Unlikely to be used by a
#' user, instead prefer [sample_nuts()]
#' @rdname samplers
#' @param verbose Boolean for whether to print ADMB output to console.
#' @seealso [sample_nuts()]
sample_admb_nuts <- function(path,
                             model,
                             iter = 2000,
                             init = NULL,
                             chain = 1,
                             thin = 1,
                             warmup = NULL,
                             seed = NULL,
                             duration = NULL,
                             control = NULL,
                             skip_optimization = TRUE,
                             verbose = TRUE,
                             admb_args = admb_args,
                             hess_step = FALSE,
                             fn_logfile = "model_output.log"){

  control <- .update_control(control)
  eps <- control$stepsize
  stopifnot(iter >= 1)
  stopifnot(warmup <= iter)
  stopifnot(duration > 0)
  stopifnot(thin >= 1)
  if(is.null(warmup)){
    stop("Must provide warmup", call. = FALSE)
  }
  if(thin < 1 | thin > iter){
    stop("Thin must be >1 and < iter", call. = FALSE)
  }
  max_td <- control$max_treedepth
  adapt_delta <- control$adapt_delta

  # Build the command to run the model
  if(skip_optimization){
    cmd <- paste0("cd ", path, " && ", model, ifelse(hess_step, "", " -nohess"),
                  " -nox -maxfn 0 -phase 1000 -nuts -mcmc ", iter)
  } else {
    cmd <- paste0("cd ", path ," && ", model, " -hbf -nuts -mcmc ", iter)
  }
  cmd <- paste0(cmd, " -warmup ", warmup, " -chain ", chain)
  if(!is.null(seed)){
    cmd <- paste0(cmd, " -mcseed ", seed)
  }
  if(hess_step){
    cmd <- paste0(cmd, " -hess_step 10 -binp ss.bar -hbf")
  }
  if(!is.null(duration)){
    cmd <- paste0(cmd, " -duration ", duration)
  }
  cmd <- paste0(cmd, " -max_treedepth ", max_td, " -adapt_delta ", adapt_delta)
  if(!is.null(eps)){
    cmd <- paste0(cmd, " -hyeps ", eps)
  }
  if(!is.null(control$adapt_init_buffer)){
    cmd <- paste0(cmd, " -adapt_init_buffer ", control$adapt_init_buffer)
  }
  if(!is.null(control$adapt_term_buffer)){
    cmd <- paste0(cmd, " -adapt_term_buffer ", control$adapt_term_buffer)
  }
  if(!is.null(control$adapt_window)){
    cmd <- paste0(cmd, " -adapt_window ", control$adapt_window)
  }
  if(!is.null(control$refresh)){
    cmd <- paste0(cmd, " -refresh ", control$refresh)
  }
  if(control$adapt_mass){
    cmd <- paste0(cmd, " -adapt_mass")
  }
  if(control$adapt_mass_dense){
    cmd <- paste0(cmd, " -adapt_mass_dense")
  }

  # Four options for metric:
  # 1. 'mle' - Use the MLE estimates in admodel.cov without mass adaptation
  # 2. If a matrix is passed, this is written to file admodel.cov and no
  #    adaptation is done
  # 3. (default) - Adaptation starting with diagonal
  # 4. Diagonal without mass adaptation
  metric <- control$metric
  stopifnot(!is.null(metric))
  if(is.matrix(metric)){
    # User defined one will be writen to admodel.cov
    cor_user <- metric / sqrt(diag(metric) %o% diag(metric))
    if(!is.positive.definite(x = cor_user))
      stop("Invalid mass matrix passed: it is not positive definite.\n",
           "Check 'metric' argument or use different option.",
           call. = FALSE)
    write_admb_cov(metric, path, hbf = 1)
    warning("admodel.cov overwritten, revert admodel_original.cov if needed")
  }else if(is.character(metric) && metric == "unit") {
    # The default: Start from unit diag.
    cmd <- paste0(cmd, " -mcdiag")
  }else if(is.character(metric) && metric == "mle") {
    # ADMB default so do nothing special. No adaptation, will use
    # estimated MLE covariance matrix in unbounded space (read from
    # admodel.cov)
  }else{
    stop("Invalid metric option", call. = FALSE)
  }
  # Write the starting values to file. A NULL value means to use the MLE,
  # so need to run model
  if(!is.null(init)){
    cmd <- paste0(cmd, " -mcpin init.pin")
    write.table(file = file.path(path, "init.pin"),
                x = unlist(init),
                row.names = FALSE,
                col.names = FALSE)
  }
  if(!is.null(admb_args)){
    cmd <- paste0(cmd, " ", admb_args)
  }

  # Run it and get results
  if(!is.null(fn_logfile)){
    cmd <- paste0(cmd, " > ", fn_logfile, " 2>&1")
  }

  time <- system.time(system_(cmd, ignore.stdout = !verbose))[3]
  if(!file.exists(file.path(path, "adaptation.csv")) ||
     !file.exists(file.path(path, "unbounded.csv"))){

    stop(file.path(path, "adaptation.csv"), " exists = ",
         file.exists(file.path(path, "adaptation.csv")),
         "\n", file.path(path, "unbounded.csv"), " exists = ",
         file.exists(file.path(path, "unbounded.csv")),
         "\nNUTS failed to run in chain ", chain,
         call. = FALSE)
  }

  sampler_params <- as.matrix(read.csv(file.path(path, "adaptation.csv")))
  unbounded <- as.matrix(read.csv(file.path(path, "unbounded.csv"),
                                  header = FALSE))
  dimnames(unbounded) <- NULL

  pars <- read_psv(path, model)

  par_names <- names(pars)
  if(!"lp__" %in% dimnames(sampler_params)[[2]]){
    # Previous version had a bug where energy__ was stored as
    # the log-posterior. So energy is wrong, but log-posterior
    # is right here.
    # warning("ADMB version <= 12.0 has a bug where the energy statistic is wrong. Please consider updating")
    pars[, "log-posterior"] <- sampler_params[, "energy__"]
  }else{
    # Later versions has a 7th column containing the LP and 6 is
    # the energy. Both energy and lp are correct
    pars[, "log-posterior"] <- sampler_params[, "lp__"]
    # Drop the lp__ here since not used and may cause issues downstream
    sampler_params <- sampler_params[, -7]
  }
  pars <- as.matrix(pars)
  # Thin samples and adaptation post hoc for NUTS
  pars <- pars[seq(1, nrow(pars), by = thin),]
  unbounded <- unbounded[seq(1, nrow(unbounded), by = thin), ]
  sampler_params <- sampler_params[seq(1, nrow(sampler_params), by = thin), ]
  time_total <- time
  time_warmup <- NA
  warmup <- warmup / thin

  list(samples = pars,
       sampler_params = sampler_params,
       time_total = time_total,
       time_warmup = time_warmup,
       warmup = warmup,
       max_treedepth = max_td,
       model = model,
       par_names = par_names,
       cmd = cmd,
       unbounded = unbounded)
}

#' Run a single random walk Metropolis chain for an ADMB model
#'
#' A low level function to run a single chain. Unlikely to be used by a
#' user, instead prefer [sample_rwm()]
#' @rdname samplers
#' @seealso [sample_rwm()]
sample_admb_rwm <- function(path,
                            model,
                            iter = 2000,
                            init = NULL,
                            chain = 1,
                            thin = 1,
                            warmup = ceiling(iter / 2),
                            seed = NULL,
                            duration = NULL,
                            control = NULL,
                            skip_optimization = TRUE,
                            verbose = TRUE,
                            admb_args = NULL,
                            fn_logfile = "model_output.log",
                            ...){

  if(any(names(control) != "refresh")){
    warning("Only refresh control argument is used with RWM, ignoring: ",
            paste(names(control)[names(control) != "refresh"], collapse = ", "),
            call. = FALSE)
  }
  refresh <- control$refresh
  metric <- "mle" # only one allowed
  stopifnot(iter >= 1)
  stopifnot(warmup <= iter)
  stopifnot(duration > 0)
  stopifnot(thin >= 1)
  if(is.null(warmup)) stop("Must provide warmup")
  if(thin < 1 | thin > iter) stop("Thin must be > 1 and < iter")

  # Build the command to run the model
  if(skip_optimization){
    cmd <- paste0("cd ", path, " && ", model," -nox -nohess -maxfn 0 -phase 1000 -rwm -mcmc ", iter)
  } else {
    cmd <- paste0("cd ", path, " && ", model,"-rwm -mcmc ", iter)
  }

  cmd <- paste0(cmd, " -mcscale ", warmup, " -chain ", chain)
  if(!is.null(seed)) cmd <- paste0(cmd, " -mcseed ", seed)
  if(!is.null(duration)) cmd <- paste0(cmd, " -duration ", duration)
  cmd <- paste0(cmd, " -mcsave ", thin)

  # Three options for metric. NULL (default) is to use the MLE estimates
  # in admodel.cov.  If a matrix is passed, this is written to file and
  # no scaling is done. Option 'unit' means identity. Note: these are
  # all in unbounded space.
  if(is.matrix(metric)){
    ## User defined one will be writen to admodel.cov
    cor_user <- metric / sqrt(diag(metric) %o% diag(metric))
    if(!is.positive.definite(x=cor_user))
      stop("Invalid mass matrix, not positive definite")
    .write.admb.cov(metric)
  } else if(is.null(metric)){
    # NULL means default of MLE
  } else if(metric == "mle"){
    # also use mle (i.e., do nothing)
  } else if(metric == "unit") {
    # Identity in unbounded space
    cmd <- paste0(cmd, " -mcdiag")
  } else {
    stop("Invalid metric option")
  }
  # Write the starting values to file. A `NULL` value means to use the MLE,
  # so need to run model
  if(!is.null(init)){
    cmd <- paste0(cmd, " -mcpin init.pin")
    write.table(file = file.path(path, "init.pin"), x = unlist(init), row.names = FALSE, col.names = FALSE)
  }
  if(!is.null(admb_args)) cmd <- paste(cmd, admb_args)

  # Run it and get results
  if(!is.null(fn_logfile)){
    cmd <- paste0(cmd, " > ", fn_logfile, " 2>&1")
  }

  time <- system_(cmd, ignore.stdout = !verbose)[3]
  unbounded_fn <- file.path(path, "unbounded.csv")
  if(!file.exists(unbounded_fn)){
    stop(paste0("RWM failed to run in chain ", chain, ". Check inputs."))
  }
  unbounded <- as.matrix(read.csv(unbounded_fn, header = FALSE))
  dimnames(unbounded) <- NULL
  pars <- read_psv(path, model)
  par_names <- names(pars)

  lp <- as.vector(read.table(file.path(path, "rwm_lp.txt"), header = TRUE)[,1])
  pars[, "log-posterior"] <- lp
  pars <- as.matrix(pars)
  # Thinning is done internally for RWM (via -mcsave) so don't need to do
  # it here
  list(samples = pars,
       sampler_params = NULL,
       time.total = time,
       time.warmup = NA,
       warmup = warmup,
       model = model,
       par_names = par_names,
       cmd = cmd,
       unbounded = unbounded)
}
