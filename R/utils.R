#' Plot marginal distributions for a fitted model
#'
#' @param fit A fitted object returned by
#'   [sample_admb()].
#' @param pars A numeric or character vector of parameters which
#'   to plot, for plotting a subset of the total (defaults to all)
#' @param mfrow A custom grid size (vector of two) to be called
#'   as `par(mfrow)`, overriding the defaults.
#' @param add.mle Whether to add marginal normal distributions
#'   determined from the inverse Hessian file
#' @param add.monitor Whether to add ESS and Rhat information
#' @param breaks The number of breaks to use in [graphics::hist()],
#'   defaulting to 30
#' @export
#'
#' @details This function plots grid cells of all parameters
#' in a model, comparing the marginal posterior histogram vs
#' the asympotitic normal (red lines) from the inverse
#' Hessian. Its intended use is to quickly gauge differences
#' between frequentist and Bayesian inference on the same
#' model.
#'
#' If `fit$monitor` exists the effective sample size
#' (ESS) and R-hat estimates are printed in the top right
#' corner. See [https://mc-stan.org/rstan/reference/Rhat.html]
#' for more information. Generally Rhat>1.05 or ESS<100 (per chain)
#' suggest inference may be unreliable.
#'
#' This function is customized to work with multipage PDFs,
#' specifically:
#' `pdf('marginals.pdf', onefile=TRUE, width=7,height=5)`
#' produces a nice readable file.
#' @examples
#' fit <- readRDS(system.file('examples', 'fit.RDS', package='adnuts'))
#' plot_marginals(fit, pars=1:2)
plot_marginals <- function(fit, pars=NULL, mfrow=NULL,
                           add.mle=TRUE, add.monitor=TRUE,
                           breaks=30){
  if(!is_adfit(fit)) stop("fit is not a valid object")
  if(!is.null(mfrow)) stopifnot(is.vector(mfrow) && length(mfrow)==2)
  stopifnot(add.mle %in% c(TRUE,FALSE))
  if(add.mle & is.null(fit$mle)) {
    add.mle <- FALSE
    warning("No MLE information found in fit$mle so cannot add")
  }
  if(!add.mle) fit$mle <- NULL
  if(!add.monitor) fit$monitor <- NULL
  par.old <- par()
  on.exit(par(mfrow=par.old$mfrow, mar=par.old$mar,
              mgp=par.old$mgp, oma=par.old$oma, tck=par.old$tck))
  posterior <- extract_samples(fit, inc_lp=FALSE)
  par.names <- names(posterior)
  if(is.null(pars)) pars <- par.names
  if(is.character(pars[1])){
    pars.ind <- match(x=pars, table=par.names)
    if(any(is.na(pars.ind))){
      warning("Some par names did not match -- dropped")
      print(pars[is.na(pars.ind)])
      pars.ind <- pars.ind[!is.na(pars.ind)]
    }
    pars <- pars.ind
  } else if(any(pars > NCOL(posterior))){
    warning("Some par numbers too big -- dropped")
    print(pars[pars > NCOL(posterior)])
    pars <- pars[ pars <=NCOL(posterior)]
  }
  n <- length(pars)
  stopifnot(is.numeric(pars[1]))
  stopifnot(ncol(posterior)>1)
  par(mar=c(1.5,0,.1,0), mgp=c(2,.4,0),
      oma=c(.25,.25,.25,.25), tck=-.02)
  if(!is.null(mfrow)){
    par(mfrow=mfrow)
  } else if(n>12){
    par(mfrow=c(4,4))
  } else if(n>9){
    par(mfrow=c(4,3))
  } else if(n>6){
    par(mfrow=c(3,3))
  } else if(n>4){
    par(mfrow=c(3,2))
  } else if(n>3){
    par(mfrow=c(2,2))
  } else {
    par(mfrow=c(1,n))
  }
  for(ii in pars){
    par <- par.names[ii]
    if(!is.null(fit$mle)){
      mle <- fit$mle$est[ii]
      se <-  fit$mle$se[ii]
      x1 <- seq(qnorm(.001, mle, se), qnorm(.999, mle, se), len=100)
      y1 <- dnorm(x1, mle, se)
    } else{
      x1 <- y1 <- NULL
    }
    tmp <- hist(posterior[,ii], plot=FALSE, breaks=breaks)
    x2 <- tmp$mids; y2 <- tmp$density
    plot(0,0, type='n', xlim=range(c(x1,x2)), yaxs='i',
         ylim=c(0, max(c(y1,y2))*1.3), axes=FALSE, ann=FALSE)
    hist(posterior[,ii], breaks=breaks, add=TRUE, yaxs='i', freq=FALSE, col=gray(.8))
    axis(1);  box(col=gray(.5));
    if(!is.null(fit$mle)) lines(x1,y1, col='red', lwd=2)
    if(!is.null(fit$monitor)){
      mon <- fit$monitor
      ## add ESS and Rhat to top right
      tmp <- par("usr"); xy <- c(.85,.88)
      text.x <- tmp[1]+xy[1]*diff(tmp[1:2])
      text.y <- tmp[3]+xy[2]*diff(tmp[3:4])
      label <- paste0('ESS=', mon[ii,'n_eff'], "\nRhat=", round(mon[ii,'Rhat'],3))
      text(x=text.x, y=text.y, labels=label, cex=.8)
    }
    mtext(paste("",par), line=-1.6, adj=0, cex=.9)
  }
}

#' Plot adaptation metrics for a fitted model.
#'
#' @param fit A fitted object returned by
#' [sample_admb()].
#' @param plot Whether to plot the results
#' @return Prints and invisibly returns a ggplot object
#'
#' @details This utility function quickly plots the adaptation output of NUTS
#' chains.
#' @importFrom rlang .data
#' @export
#' @examples
#' fit <- readRDS(system.file('examples', 'fit.RDS', package='adnuts'))
#' plot_sampler_params(fit)
plot_sampler_params <- function(fit, plot=TRUE){
  if(!requireNamespace("ggplot2", quietly=TRUE))
    stop("ggplot2 package not found")
  sp <- adnuts::extract_sampler_params(fit, inc_warmup=TRUE)
  sp.long <-
    data.frame(iteration=sp$iteration, chain=factor(sp$chain),
               value=c(sp$accept_stat__, log(sp$stepsize__),
                       sp$n_leapfrog__, sp$divergent__, sp$energy__),
               variable=rep(c('accept_stat', 'log_stepsize',
                              'n_leapfrog', 'divergent',
                              'energy'), each=nrow(sp)))
  g <- ggplot2::ggplot(sp.long, ggplot2::aes(.data$iteration, y=.data$value, color=.data$chain)) +
    ggplot2::geom_point(alpha=.5) +
    ggplot2::facet_wrap('variable', scales='free_y', ncol=1) + ggplot2::theme_bw()
  if(plot) print(g)
  invisible(g)
}

#' Function to generate random initial values from a previous fit using
#' adnuts
#'
#' @param fit An outputted list from [sample_admb()]
#' @param chains The number of chains for the subsequent run, which
#'   determines the number to return.
#' @return A list of lists which can be passed back into
#'   [sample_admb()].
#' @export
sample_inits <- function(fit, chains){
  post <- extract_samples(fit)
  ind <- sample(seq_len(nrow(post)), size = chains)
  lapply(ind, function(i) as.numeric(post[i, ]))
}

#' Update algorithm for mass matrix.
#'
#' @param fn The current fn function.
#' @param gr The current gr function
#' @param y.cur The current parameter vector in unrotated (Y) space.
#' @param M The new mass matrix
.rotate_space <- function(fn, gr, M,  y.cur){
  ## Rotation done using choleski decomposition
  ## First case is a dense mass matrix
  if(is.matrix(M)){
    chd <- t(chol(M))               # lower triangular Cholesky decomp.
    chd.inv <- solve(chd)               # inverse
    ## Define rotated fn and gr functions
    fn2 <- function(x) fn(chd %*% x)
    gr2 <- function(x) {as.vector( gr(chd %*% x) %*% chd )}
    ## Now rotate back to "x" space using the new mass matrix M
    x.cur <- as.numeric(chd.inv %*% y.cur)
  } else if(is.vector(M)){
    chd <- sqrt(M)
    fn2 <- function(x) fn(chd * x)
    gr2 <- function(x) as.vector(gr(chd * x) ) * chd
    ## Now rotate back to "x" space using the new mass matrix M. M is a
    ## vector here. Note the big difference in efficiency without the
    ## matrix operations.
    x.cur <- (1/chd) * y.cur
  } else {
    stop("Mass matrix must be vector or matrix")
  }
  ## Redefine these functions
  ## Need to adjust the current parameters so the chain is
  ## continuous. First rotate to be in Y space.
  return(list(gr2=gr2, fn2=fn2, x.cur=x.cur, chd=chd))
}

#' Convert adnuts fit (named list) into a \code{shinystan} object.
#'
#' @details The shinystan packages provides several conversion functions
#' for objects of different types, such as stanfit classes (Stan ouput)
#' and simple arrays. For the latter, option NUTS information, such as
#' [sampler_params()] can be passed. This function essentially extends
#' the functionality of [as.shinystan()] to work specifically with
#' fits from adnuts (TMB or ADMB). The user can thus explore their model
#' with [launch_shinystan(.as.shinyadnuts(fit))] in the same way
#' that Stan models are examined.
#' @param fit Output list from [sample_admb()].
#' @seealso [launch_shinyadmb()]
#' @return An S4 object of class shinystan. Depending on the algorithm
#'   used, this list will have slight differences.
.as.shinyadnuts <- function(fit){
  if(fit$algorithm=="NUTS"){
    sso <- with(fit, shinystan::as.shinystan(samples, warmup=warmup, max_treedepth=max_treedepth,
             sampler_params=sampler_params, algorithm='NUTS', model_name=model))
  } else if(fit$algorithm=="HMC"){
    sso <- with(fit, shinystan::as.shinystan(samples, warmup=warmup,
             sampler_params=sampler_params, algorithm='HMC', model_name=model))
  } else {
    sso <- with(fit, shinystan::as.shinystan(samples, warmup=warmup,
             algorithm='RWM', model_name=model))
  }
  return(invisible(sso))
}

#' Launch shinystan for a TMB fit.
#'
#' @param fit A named list returned by [sample_tmb()]
#' @seealso [launch_shinyadmb()]
launch_shinytmb <- function(fit){
  shinystan::launch_shinystan(.as.shinyadnuts(fit))
}

#' Launch shinystan for an ADMB fit.
#'
#' @param fit A named list returned by [sample_admb()]
#' @seealso [launch_shinytmb()]
#' @export
launch_shinyadmb <- function(fit){
  shinystan::launch_shinystan(.as.shinyadnuts(fit))
}

#' Extract sampler parameters from a fit.
#'
#' Extract information about NUTS trajectories, such as acceptance ratio
#' and treedepth, from a fitted object.
#'
#' @details Each trajectory (iteration) in NUTS has associated information
#' about the trajectory: stepsize, acceptance ratio, treedepth, and number of
#' leapfrog steps. This function extracts these into a data.frame, which
#' may be useful for diagnosing issues in certain cases. In general, the
#' user should not need to examine them, or preferably should via
#' [plot_sampler_params()] or [launch_shinyadmb()].
#'
#' @param fit A list returned by [sample_admb()]
#' @param inc_warmup Whether to extract the `warmup` samples or not
#' (default). `warmup` samples should never be used for inference, but may
#'  be useful for diagnostics
#' @return An invisible data.frame containing samples (rows) of each
#' parameter (columns). If multiple chains exist they will be [rbind()]ed
#' together
#' @seealso [launch_shinyadmb()]
#' @export
#' @examples
#' fit <- readRDS(system.file('examples', 'fit.RDS', package='adnuts'))
#' sp <- extract_sampler_params(fit, inc_warmup=TRUE)
#' str(sp)
extract_sampler_params <- function(fit, inc_warmup=FALSE){
  x <- fit$sampler_params
  if(!is.list(x)) stop("fit$sampler_parameters is not a list -- valid fit object?")
  if(inc_warmup){
    ind <- 1:dim(x[[1]])[1]
    its <- 1:length(ind)
  } else{
    ind <- -(1:fit$warmup)
    its <- (1:length(ind)) + fit$warmup
  }
  y <- do.call(rbind, lapply(1:length(x), function(i)
    cbind(chain=i, iteration=its, x[[i]][ind,])))
  return(invisible(as.data.frame(y)))
}

#' Write matrix of samples to a binary .psv file.
#'
#' @details Useful to combine multiple MCMC runs together into a single
#' .psv file which can then be executed with `-mceval`.
#' @param fn Model executable name
#' @param samples A matrix or data.frame of samples, each column is a
#'  parameter, each row a sample.
write_psv <- function(fn,
                      samples,
                      path = NULL){

  if(is.null(path)){
    stop("`path` cannot be `NULL`",
         call. = FALSE)
  }
  if(!dir.exists(path)){
    stop("Directory `path` = `", path, "` does not exist",
         call. = FALSE)
  }

  samples <- as.matrix(samples)
  psv <- file.path(path, paste0(fn, ".psv"))
  con <- file(psv, "wb")
  on.exit(close(con), add = TRUE)
  writeBin(object = ncol(samples), con = con)
  writeBin(object = as.vector(t(samples)), con = con)
  invisible()
}

#' Find the model executable name on your filesystem
#'
#' @description
#' Find the model executable name in the local directory (`loc_dir`).
#' If not there, search the PATH.
#' If found more than once in the PATH, the first one will be used and a
#' warning issued.
#' If not on the PATH, throw an error and stop code execution.
#' If using Linux or OSX and the file permissions do not allow the user to
#' execute, throw an error and stop code execution.
#'
#' @param model The name of the executable (without any extension)
#' @param loc_dir The path to check for existence of the executable before
#' checking the path
#' @param path_only If `TRUE`, ignore `loc_dir` and search the PATH only
#'
#' @return The executable as it should be called in the OS. It is either
#' just as it was given if on the PATH or in the `loc_dir` directory on a
#' Windows machine, or prepended with './' if in the `loc_dir` directory on a
#' non-Windows machine. If not found anywhere, an error will be thrown and
#' code execution will stop
#' @export
get_model_executable <- function(model,
                                 loc_dir = getwd(),
                                 path_only = FALSE){

  stopifnot(is.character(loc_dir))
  stopifnot(is.character(model))

  if(!path_only && !dir.exists(loc_dir)){
    stop("Directory `", loc_dir, "` does not exist. Check ",
         "argument `loc_dir`",
         call. = FALSE)
  }

  check_dot <- grep("\\.", model)
  if(length(check_dot)){
    stop("model name cannot include a '.', Do not include extension in ",
         "model name",
         call. = FALSE)
  }

  os <- get_os()
  loc_dir_files <- dir(loc_dir)
  model_regex <- gsub("\\+", "\\\\+", model)
  model_name <- grep(paste0("^",
                            model_regex,
                            ifelse(os == "windows", ".exe", ""), "$"),
                     loc_dir_files,
                     value = TRUE)

  if(!path_only && length(model_name)){
    if(os != "windows"){
      model <- paste0("./", model)
    }
    found_model_loc <- "loc_dir"
  }else{
    # No local executable, check PATH for one by using which or where
    model_path <- tryCatch(system(paste("which", model),
                                  intern = TRUE),
                           warning = function(w){
                             NA
                           })
    if(is.na(model_path[1])){
      stop(model, " not found locally or on the PATH.", call. = FALSE)
    }
    if(length(model_path) > 1){
      warning("The model executable ", model,
              " was found in more than one place on your PATH.\n",
              "Using the first one found (", model_path[1], ").")
    }
    found_model_loc <- "path"
  }
  # Check to see if the file is executable (LInux and OSX only)
  if(os != "windows"){
    stat <- system(paste("stat", model_path[1]), intern = TRUE)
    perms <- grep("^Access:.*Uid.*$", stat, value = TRUE)
    perms <- gsub("\\).*", "", perms)
    perms <- gsub(".*/","",perms)
    first_perm <- substr(perms, 1, 1)
    user_exec_perm <- substr(perms, 4, 4)
    if(first_perm == "d"){
      stop("The model name given (", model_path[1], ") is a directory and not executable.",
           call. = FALSE)
    }else if(first_perm == "p"){
      stop("The model name given (", model_path[1], ") is a named pipe and not executable.",
           call. = FALSE)
    }else if(first_perm == "s"){
      stop("The model name given (", model_path[1], ") is a socket and not executable.",
           call. = FALSE)
    }else if(first_perm %in% c("b", "c")){
      stop("The model name given (", model_path[1], ") is a device block and not executable.",
           call. = FALSE)
    }else if(user_exec_perm != "x"){
      stop("The file found (", model_path[1], ") is not executable.\nCheck name and use chmod to ",
           "change permissions if name is correct.",
           call. = FALSE)
    }
  }
  if(found_model_loc == "loc_dir"){
    message("Using model found in the local directory given by the 'loc_dir' argument: ", model)
  }else if(found_model_loc == "path"){
    message("Using model executable found in the PATH: ", model)
  }
  model
}

#' Call [shell()] or [system()] depending on the Operating System
#'
#' @param ... Pass all arguments to the command function
#'
#' @return The output from the command function called
#' @export
system_ <- function(...){
  if(get_os() == "windows"){
    shell(...)
  }else{
    system(...)
  }
}
