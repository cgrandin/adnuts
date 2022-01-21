#' A wrapper for running ADMB models in parallel
#'
#' @param parallel_number The chain number
#' @param path The path name to append the chain number to
#' @param algorithm NUTS or RWM
#' @param ... Arguments passed to the [fit()] function
#'
#' @return Output from the [fit()] function
#' @export
sample_admb_parallel <- function(parallel_number,
                                 path,
                                 algorithm,
                                 ...){

  chain_dir <- file.path(path, paste0("chain_", parallel_number))
  if(dir.exists(chain_dir)){
    unlink(chain_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(chain_dir)
  if(!dir.exists(chain_dir)){
    stop("Could not create directory: ", chain_dir, call. = FALSE)
  }

  trash <- file.copy(from = list.files(path, full.names = TRUE),
                     to = chain_dir)

  if(algorithm == "NUTS"){
    fit <- sample_admb_nuts(path = chain_dir,
                            chain = parallel_number,
                            ...)
  }else if(algorithm == "RWM"){
    fit <- sample_admb_rwm(path = chain_dir,
                           chain = parallel_number,
                           ...)
  }else{
    stop("Algorithm must be either 'NUTS' or 'RWM'", call. = FALSE)
  }
  unlink(chain_dir,
         recursive = TRUE,
         force = TRUE)

  fit
}

# A wrapper for running TMB models in parallel
sample_tmb_parallel <-  function(parallel_number, obj, init, path,
                                 algorithm, lower, upper, seed, laplace, ...){
  ## Each node starts in a random work directory. Rebuild TMB model obj so
  ## can link it in each session.
  setwd(path)
  dyn.load(TMB::dynlib(obj$env$DLL))
  ## Use 'shape' attribute to obtain full length of 'map'ped parameters.
  map.index <- which(names(obj$env$parameters) %in% names(obj$env$map))
  new.par <- obj$env$parameters
  new.par[map.index] <- lapply(obj$env$parameters[map.index], function(x) attr(x, "shape"))
  obj <- TMB::MakeADFun(data=obj$env$data, parameters=new.par, random=obj$env$random,
                        map=obj$env$map, DLL=obj$env$DLL, silent=TRUE)
  obj$env$beSilent()

  ## Ignore parameters declared as random? Borrowed from tmbstan.
  if(laplace){
    par <- obj$env$last.par.best[-obj$env$random]
    fn0 <- obj$fn
    gr0 <- obj$gr
  } else {
    par <- obj$env$last.par.best
    fn0 <- obj$env$f
    gr0 <- function(x) obj$env$f(x, order=1)
  }
  ## Parameter constraints, if provided, require the fn and gr functions to
  ## be modified to account for differents in volume. There are four cases:
  ## no constraints, bounded below, bounded above, or both (box
  ## constraint).
  bounded <- !(is.null(lower) & is.null(upper))
  if(bounded){
    if(is.null(lower)) lower <- rep(-Inf, len=length(upper))
    if(is.null(upper)) upper <- rep(Inf, len=length(lower))
    cases <- .transform.cases(lower, upper)
    fn <- function(y){
      x <- .transform(y, lower, upper, cases)
      scales <- .transform.grad(y, lower, upper, cases)
      -fn0(x) + sum(log(scales))
    }
    gr <- function(y){
      x <- .transform(y, lower, upper, cases)
      scales <- .transform.grad(y, lower, upper, cases)
      scales2 <- .transform.grad2(y, lower, upper, cases)
      -as.vector(gr0(x))*scales + scales2
    }
    ## Don't need to adjust this b/c init is already backtransformed in
    ## sample_tmb.
    ## init <- .transform.inv(x=unlist(init), a=lower, b=upper, cases=cases)
  } else {
    fn <- function(y) -fn0(y)
    gr <- function(y) -as.vector(gr0(y))
  }
  if(algorithm=="NUTS")
    fit <- sample_tmb_nuts(chain=parallel_number, fn=fn, gr=gr,
                           init=init, seed=seed, ...)
  if(algorithm=="RWM")
    fit <- sample_tmb_rwm(chain=parallel_number, fn=fn, init=init,
                          seed=seed, ...)
  return(fit)
}
