#' Plot marginal distributions for a fitted model
#'
#' @param fit A fitted object returned by
#'   [sample_admb()].
#' @param pars A numeric or character vector of parameters which
#'   to plot, for plotting a subset of the total (defaults to all)
#' @param mfrow A custom grid size (vector of two) to be called
#'   as `par(mfrow)`, overriding the defaults.
#' @param add_mle Whether to add marginal normal distributions
#'   determined from the inverse Hessian file
#' @param add_monitor Whether to add ESS and Rhat information
#' @param breaks The number of breaks to use in [graphics::hist()],
#'   defaulting to 30
#' @export
#'
#' @details This function plots grid cells of all parameters
#' in a model, comparing the marginal posterior histogram vs
#' the asympototic normal (red lines) from the inverse
#' Hessian. Its intended use is to quickly gauge differences
#' between frequentist and Bayesian inference on the same
#' model.
#'
#' If `fit$monitor` exists the effective sample size
#' (ESS) and R-hat estimates are printed in the top right
#' corner. See [https://mc-stan.org/rstan/reference/Rhat.html]
#' for more information. Generally Rhat > 1.05 or ESS < 100 (per chain)
#' suggest inference may be unreliable.
#'
#' This function is customized to work with multipage PDFs,
#' specifically:
#' `pdf('marginals.pdf', onefile=TRUE, width=7,height=5)`
#' produces a nice readable file.
#' @examples
#' fit <- readRDS(system.file('examples', 'fit.RDS', package='adnuts'))
#' plot_marginals(fit, pars=1:2)
plot_marginals <- function(fit,
                           pars = NULL,
                           mfrow = NULL,
                           add_mle = TRUE,
                           add_monitor = TRUE,
                           breaks = 30){

  if(!is_adfit(fit)){
    stop("fit is not a valid object",
         call. = FALSE)
  }
  if(!is.null(mfrow)){
    stopifnot(is.vector(mfrow) && length(mfrow) == 2)
  }
  stopifnot(add_mle %in% c(TRUE, FALSE))

  if(add_mle && is.null(fit$mle)) {
    add_mle <- FALSE
    warning("No MLE information found in `fit$mle` so cannot add")
  }else{
    fit$mle <- NULL
  }
  if(!add_monitor){
    fit$monitor <- NULL
  }
  par_old <- par()
  on.exit(par(mfrow = par_old$mfrow,
              mar = par_old$mar,
              mgp = par_old$mgp,
              oma = par_old$oma,
              tck = par_old$tck))

  posterior <- extract_samples(fit, inc_lp = FALSE)
  par_names <- names(posterior)
  if(is.null(pars)) pars <- par_names
  if(is.character(pars[1])){
    par_ind <- match(x = pars,
                      table = par_names)
    if(any(is.na(par_ind))){
      warning("Some par names did not match, removing them")
      print(pars[is.na(par_ind)])
      par_ind <- par_ind[!is.na(par_ind)]
    }
    pars <- par_ind
  }else if(any(pars > ncol(posterior))){
    warning("Some par numbers too big, removing them")
    print(pars[pars > ncol(posterior)])
    pars <- pars[ pars <= ncol(posterior)]
  }
  n <- length(pars)
  stopifnot(is.numeric(pars[1]))
  stopifnot(ncol(posterior) > 1)
  par(mar = c(1.5, 0, 0.1, 0),
      mgp = c(2, 0.4, 0),
      oma = c(0.25, 0.25, 0.25, 0.25),
      tck = -0.02)
  if(!is.null(mfrow)){
    par(mfrow = mfrow)
  }else if(n > 12){
    par(mfrow = c(4, 4))
  }else if(n > 9){
    par(mfrow = c(4, 3))
  }else if(n > 6){
    par(mfrow = c(3, 3))
  }else if(n > 4){
    par(mfrow = c(3, 2))
  }else if(n > 3){
    par(mfrow = c(2, 2))
  }else{
    par(mfrow = c(1, n))
  }
  for(ii in pars){
    par <- par_names[ii]
    if(!is.null(fit$mle)){
      mle <- fit$mle$est[ii]
      se <- fit$mle$se[ii]
      x1 <- seq(qnorm(0.001,
                      mle,
                      se),
                qnorm(0.999,
                      mle,
                      se),
                len = 100)
      y1 <- dnorm(x1, mle, se)
    }else{
      x1 <- y1 <- NULL
    }
    tmp <- hist(posterior[ ,ii],
                plot = FALSE,
                breaks = breaks)
    x2 <- tmp$mids
    y2 <- tmp$density
    plot(0,
         0,
         type = "n",
         xlim = range(c(x1, x2)),
         yaxs = "i",
         ylim = c(0,
                  max(c(y1, y2)) * 1.3),
         axes = FALSE,
         ann = FALSE)
    hist(posterior[ ,ii],
         breaks = breaks,
         add = TRUE,
         yaxs = "i",
         freq = FALSE,
         col = gray(0.8))
    axis(1)
    box(col = gray(0.5));
    if(!is.null(fit$mle)){
      lines(x1,
            y1,
            col = "red",
            lwd = 2)
    }
    if(!is.null(fit$monitor)){
      mon <- fit$monitor
      # Add ESS and Rhat to top right
      tmp <- par("usr")
      xy <- c(0.85, 0.88)
      text_x <- tmp[1] + xy[1] * diff(tmp[1:2])
      text_y <- tmp[3] + xy[2] * diff(tmp[3:4])
      label <- paste0("ESS=",
                      mon[ii, "n_eff"],
                      "\nRhat=",
                      round(mon[ii, "Rhat"],
                            3))
      text(x = text_x,
           y = text_y,
           labels = label,
           cex = 0.8)
    }
    mtext(paste("",
                par),
          line = -1.6,
          adj = 0,
          cex = 0.9)
  }
}
