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
extract_sampler_params <- function(fit, inc_warmup = FALSE){
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
    cbind(chain = i, iteration = its, x[[i]][ind,])))
  return(invisible(as.data.frame(y)))
}

#' Call [shell()] or [system()] depending on the Operating System
#'
#' @param ... Arguments to pass to either [shell()] or [system()]
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
