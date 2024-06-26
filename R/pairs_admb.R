#' Plot pairwise parameter posteriors and optionally the MLE points and
#' confidence ellipses.
#'
#' @author Cole Monnahan
#' @details This function is modified from the [pairs()]
#' code to work specifically with fits from the
#' [adnuts] package using either the NUTS or RWM MCMC
#' algorithms. If an invertible Hessian was found (in `fit$mle`)
#' then estimated covariances are available to
#' compare and added automatically (red ellipses). Likewise, a
#' `monitor` object from [rstan::monitor()] is attached as
#' `fit$monitor` and provides effective sample sizes (ESS)
#' and Rhat values. The ESS are used to potentially order the
#' parameters via argument `order`, but also printed on
#' the diagonal.
#' @param fit A list as returned by [sample_admb()]
#' @param pars A vector of parameter names or integers
#' representing which parameters to subset. Useful if the model
#' has a larger number of parameters and you just want to show
#' a few key ones.
#' @param label.cex Control size of outer and diagonal labels (default 1)
#' @param order The order to consider the parameters. Options are
#' `NULL` (default) to use the order declared in the model, or
#' `slow` and `fast` which are based on the effective sample
#' sizes ordered by slowest or fastest mixing respectively. See
#' example below  for usage.
#' @param diag What type of plot to include on the diagonal,
#' options are `acf` which plots the autocorrelation function
#' [stats::acf()], `hist` shows marginal posterior histograms ([graphics::hist()]),
#' and `trace` the trace plot
#' @param acf.ylim If using the [stats::acf()] function on the diagonal,
#' specify the y limit. The default is `c(-1, 1)`
#' @param ymult A vector of length `ncol(posterior)` specifying how
#' much room to give when using the `hist` option for the
#' diagonal. For use if the label is blocking part of the
#' plot. The default is `1.3` for all parameters
#' @param axis.col Color of axes
#' @param limits A list containing the ranges for each parameter
#' to use in plotting
#' @param add.monitor Boolean whether to print effective sample
#' @param add.mle Boolean whether to add 95% confidence ellipses
#' @param unbounded Whether to use the bounded or unbounded version
#' of the parameters `size` (ESS) and `Rhat` values on the diagonal
#' @param ... Arguments to be passed to plot call in lower
#' diagonal panels
#' @return Produces a plot, and returns nothing
#' @export
#' @examples
#' \dontrun{
#' fit <- readRDS(system.file('examples', 'fit.RDS', package='adnuts'))
#' pairs_admb(fit)
#' pairs_admb(fit, pars=1:2)
#' pairs_admb(fit, pars=c('b', 'a'))
#' pairs_admb(fit, pars=1:2, order='slow')
#' pairs_admb(fit, pars=1:2, order='fast')
#' }
pairs_admb <- function(fit,
                       order = NULL,
                       diag = c("trace", "acf", "hist"),
                       acf.ylim = c(-1, 1),
                       ymult = NULL,
                       axis.col = gray(.5),
                       pars = NULL,
                       label.cex = 0.8,
                       limits = NULL,
                       add.mle = TRUE,
                       add.monitor = TRUE,
                       unbounded = FALSE,
                       ...){

  if(!is_adfit(fit))
    stop("Argument 'fit' is not a valid object returned by 'sample_admb'")
  if(unbounded | !add.mle){
    mle <- NULL
  } else {
    mle <- fit$mle
  }
  posterior <- extract_samples(fit, inc_lp=TRUE, unbounded=unbounded)
  chains <- rep(1:dim(fit$samples)[2], each=dim(fit$samples)[1]-fit$warmup)
  divs <- if(fit$algorithm == "nuts")
            extract_sampler_params(fit)$divergent__ else NULL
  ptcex <- .2
  divcex <- .75
  chaincols <- 1:length(unique(chains))
  ## reset to old par when exiting
  old.par <- par(no.readonly=TRUE)
  on.exit(par(old.par))
  diag <- match.arg(diag)
  par.names <- names(posterior)
  ess <- fit$monitor$n_eff
  Rhat <- fit$monitor$Rhat
  if(is.null(ess))
    warning("No monitor information found in fitted object so ESS and Rhat not available. See details of help.")
  if(!is.null(order)){
    if(! order %in% c('slow', 'fast')){
      stop("Invalid 'order' argument, should be 'slow', 'fast', or NULL")
    }
    if(is.null(ess)){
      stop("No effective sample sizes found so cannot order by slow/fast.")
    }
    if(!is.numeric(pars[1])){
      warning("Ignoring 'order' argument because parameter names supplied in 'pars'")
    } else {
      # Get slowest or fastest parameter names
      ind <- order(ess, decreasing=(order=='fast'))
      par.names <- par.names[ind]
    }
  }
  # if(!(NCOL(posterior) %in% c(mle$nopar, mle$nopar + 1)))
  #   stop("Number of parameters in posterior and mle not the same")
  # pars will either be NULL, so use all parameters. OR a vector of
  # indices OR a vector of characters. Want to force lp__ to be the very
  # last one stylistically, and b/c there is no ellipse for it.
  if(is.null(pars)){
    # Use all or first 10
    if(NCOL(posterior)>10){
      warning("Only showing first 10 parameters, use 'pars' argument to adjust")
      pars <- par.names[1:10]
    } else {
      pars <- par.names[1:NCOL(posterior)]
    }
  } else if(is.numeric(pars[1])){
    # Index can be used instead of character names. Note this
    # can be sorted from above
    pars <- par.names[pars]
  }
  # Now pars is character and possibly sorted by fast/slow
  pars.bad <- match(x=pars, table=names(posterior))
  if(any(is.na(pars.bad))){
    warning("Some par names did not match -- dropped")
    print(pars.bad)
    pars <- pars[!is.na(pars.bad)]
  }
  # Converts character to index which is used throughout to
  # subset when looping
  pars.ind <- match(x=pars, table=names(posterior))
  n <- length(pars.ind)
  n.mle <- ifelse(is.null(mle), 0, nrow(mle$cor))
  if(n==1) stop("This function is only meaningful for >1 parameter")
  if(is.null(ymult)) ymult <- rep(1.3, n)
  # If no limits given, calculate the max range of the posterior samples and
  # parameter confidence interval.
  if(is.null(limits)){
    limits <- list()
    for(i in 1:n){
      if(pars.ind[i]<=n.mle){
        limit.temp <- mle$est[pars.ind[i]] +
          c(-1,1)*1.96*mle$se[pars.ind[i]]
      } else {
        limit.temp <- c(NA,NA)
      }
      # Multiplier for the ranges, adjusts the whitespace around the plots
      min.temp <- min(posterior[,pars.ind[i]], limit.temp[1], na.rm=TRUE)
      max.temp <- max(posterior[,pars.ind[i]], limit.temp[2], na.rm=TRUE)
      margin <- .15*(max.temp-min.temp)
      limits[[i]] <- c(min.temp-margin, max.temp+margin)
    }
  }
  # Change posterior point look depending on how many samples. Makes
  # easier to read.
  N <- NROW(posterior)
  mypch <- 16; mycol <- 1
  if(N>=1000){
    mycol <- rgb(0,0,0,.5)
  } else if(N>=10000){
    mycol <- rgb(0,0,0,.05)
  }
  if(is.null(divs)) divs <- rep(0, N)
  par(mfrow=c(n,n), mar=0*c(.1,.1,.1,.1), yaxs="i", xaxs="i", mgp=c(.25, .25,0),
      tck=-.02, cex.axis=.65, col.axis=axis.col, oma=c(2, 2, 2,2))
  temp.box <- function() box(col=axis.col, lwd=.5)
  # Row and col here are not the posterior, but the matrix of pairwise
  # combinations
  for(row in 1:n){
    for(col in 1:n){
      ii <- pars.ind[row]
      jj <- pars.ind[col]
      ## Diagonal, so add user choice
      if(row==col){
        if(diag=="hist"){
          h <- hist(posterior[,ii], plot=F)
          # Annoyingling you can't pass NULL to xlim in hist. So
          # have to split up for two cases depending on limits.
          if(is.null(limits)){
            hist(posterior[,ii], axes=F, freq=FALSE, ann=F,
                 ylim=c(0, ymult[row]*max(h$density)),
                 col=gray(.8), border=gray(.5))
          } else {
            # Else use the user provided limits
            hist(posterior[,ii], axes=F, freq=FALSE, ann=F,
                 ylim=c(0, ymult[row]*max(h$density)),
                 col=gray(.8), border=gray(.5), xlim=limits[[row]])
          }
          temp.box()
        } else if(diag=="acf") {
          acf(posterior[,ii], axes=F, ann=F, ylim=acf.ylim)
          temp.box()
        } else if(diag=="trace") {
          # Trace plots for each chain separately
          xlim <- c(1, length(chains[chains==1]))
          plot(x=0, y=0,  type="n", axes=FALSE,
               ann=FALSE, ylim=limits[[row]], xlim=xlim)
          for(ll in unique(chains)){
            lines(posterior[chains==ll,ii], col=chaincols[ll], lwd=.1)
          }
          temp.box()
        }
        # Add ESS and Rhat info to diagonal
        if(!is.null(ess) & !is.null(Rhat) & add.monitor)
          mtext(paste0('ESS=', round(ess[ii], 0), " Rhat=", format(round(Rhat[ii],2),nsmall=2)),
                cex=.8*label.cex, line=-1)
      }
      # If lower triangle and covariance known, add scatterplot
      if(row>col){
        par(xaxs="r", yaxs="r")
        plot(x=posterior[,jj], y=posterior[,ii], axes=FALSE, ann=FALSE,
             pch=mypch, cex=ptcex, col=mycol, xlim=limits[[col]],
             ylim=limits[[row]], ...)
        # replot divegences on top so they are always visible
        points(x=posterior[which(divs==1),jj], y=posterior[which(divs==1),ii],
               pch=mypch, cex=divcex, col='green')
        # can only add MLE stuff if not lp__ parameter which
        # doesn'th ave one
        if(ii<=n.mle & jj <=n.mle){
          # Add bivariate 95% normal levels from MLE
          points(x=mle$est[jj], y=mle$est[ii],
                 pch=16, cex=.5, col='red')
          # Get points of a bivariate normal 95% confidence contour
          if(!requireNamespace("ellipse", quietly=TRUE)){
            warning("ellipse package needs to be installed to show ellipses")
          } else {
            ellipse.temp <- ellipse(x=mle$cor[jj, ii],
                                    scale=mle$se[c(jj, ii)],
                                    centre= mle$est[c(jj, ii)], npoints=1000,
                                    level=.95)
            lines(ellipse.temp , lwd=.5, lty=1, col="red")
          }
        }
        par(xaxs="i", yaxs="i")
        temp.box()
      }
      if(row<col){
        # If upper triangle add text showing the empirical correlation
        plot(0,0,type="n", xlim=c(0,1), ylim=c(0,1), axes=F,ann=F)
        temp.cor <- round(cor(posterior[,c(ii,jj)])[1,2],2)
        # Set a minimum limit for this, so they're still
        # visible, but still a function of correlation. This
        # might need to be dynamic with n.
        legend("center", legend=NA, title=temp.cor,
               cex=(3*abs(temp.cor)+.25)*.5, bty='n')
        # text(.5,.5, labels=temp.cor, cex=(3*abs(temp.cor)+.5)*.9,
        #      col=1)
        temp.box()
      }
      # Add special cases of axes on the ends
      if(row==n) {
        par( mgp=c(.05, ifelse(col %% 2 ==0, 0, .5),0) )
        axis(1, col=axis.col, lwd=.5)
      }
      if(col==1 & row >1) {
        par( mgp=c(.05, ifelse(row %% 2 ==1, .15, .65),0) )
        axis(2, col=axis.col, lwd=.5)
      }
      if(col==1 & row ==1){
        par( mgp=c(.05, ifelse(row %% 2 ==1, .15, .65),0) )
        axis(2, col=axis.col, lwd=.5)
      }
      if(row==1) mtext(pars[col], line=ifelse(col %% 2 ==1, .1, 1.1),
                       cex=label.cex)
      if(col==n)
        mtext(pars[ii], side=4, line=ifelse(row %% 2 ==1, 0, 1), cex=label.cex)
    }
  }
}
