#' Run a single random walk Metropolis chain for an ADMB model
#'
#' A low level function to run a single chain. Unlikely to be used by a
#' user, instead prefer [sample_rwm()]
#' @rdname samplers
#' @seealso [sample_rwm()]
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate
#' @importFrom microbenchmark microbenchmark
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
  metric <- "mle"
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

  runtime <- microbenchmark(system_(cmd, ignore.stdout = !verbose),
                            times = 1)

  unbounded_fn <- file.path(path, "unbounded.csv")
  if(!file.exists(unbounded_fn)){
    stop(paste0("RWM failed to run in chain ", chain, ". Check inputs."))
  }
  unbounded <- as.matrix(read.csv(unbounded_fn, header = FALSE))
  dimnames(unbounded) <- NULL
  pars <- read_psv(path)
  par_names <- names(pars)

  lp <- as.vector(read.table(file.path(path, "rwm_lp.txt"), header = TRUE)[,1])
  pars <- pars |>
    as_tibble() |>
    mutate(`log-posterior` = lp) |>
    as.matrix()
  # Thinning is done internally for RWM (via -mcsave) so don't need to do
  # it here
  list(samples = pars,
       sampler_params = NULL,
       runtime = runtime,
       warmup = warmup,
       model = model,
       par_names = par_names,
       cmd = cmd,
       unbounded = unbounded)
}
