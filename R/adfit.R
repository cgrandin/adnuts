#' Constructor for the `adfit` (Automatic differentiation AD) class
#'
#' @param x Output list from [sample_admb()]
#' @return An object of class `adfit`
#'
#' @export
adfit <- function(x){

  stopifnot(is.list(x))
  if(is.null(x$samples)){
    stop("Samples missing from fit",
         call. = FALSE)
  }
  if(is.null(x$algorithm)){
    stop("Algorithm missing from fit",
         call. = FALSE)
  }

  class(x) <- c(class(x), "adfit")

  x
}

#' Does object `x` inherit class `adfit`?
#'
#' @param x Output list from [sample_admb()]
#'
#' @return Logical
#'
#' @export
is_adfit <- function(x){

  inherits(x, "adfit")
}

#' Convert object of class `adfit` to a `data.frame`
#'
#' @details
#' This is just a 1-1 wrapper for [extract_samples()]

#' @param x Output list from [sample_admb()]
#' @return A data frame with parameters as columns and samples as
#' rows.
#' @export
as_data_frame_adfit <- function(x){

  extract_samples(x)
}

#' Plot object of class `adfit`
#'
#' @details
#' This is just a 1-1 wrapper for [plot_marginals()]

#' @param x Output list from [sample_admb()]
#'
#' @return Nothing
#' @export
plot_adfit <- function(x){

  plot_marginals(x)
}

#' Print summary of object of class `adfit`
#'
#' @details
#' This is just a 1-1 wrapper for [print_adfit()]
#'
#' @param x Output list from [sample_admb()]
#'
#' @return Nothing
#'
#' @export
summary_adfit <- function(x){

  print_adfit(x)
}

#' Print summary of an `adfit` object
#'
#' @param x Output list from [sample_admb()]
#' @param ... Absorb arguments intended for other functions
#'
#' @return Nothing

#' @export
print_adfit <- function(x, ...){

  iter <- dim(x$samples)[1]
  chains <- dim(x$samples)[2]
  pars <- dim(x$samples)[3]-1
  samples <- (iter - x$warmup) * chains

  cat(paste0("Model `", x$model,"`", " has ", pars, " pars, and was fit using `",
             x$algorithm, "`` with `", iter, "` iter and `", chains,
             "` chains\n"))

  rt <- sum(x$runtime) / chains

  ru <- "seconds"

  if(rt > 60 * 60 * 24) {
    rt <- rt / (60 * 60 *24)
    ru <- "days"
  } else if(rt > 60 * 60) {
    rt <- rt / (60 * 60)
    ru <- "hours"
  } else if(rt > 60){
    rt <- rt / 60
    ru <- "minutes"
  }

  cat("Average run time per chain was", round(rt, 2),  ru, "\n")

  if(!is.null(x$monitor)){
    min_ess <- min(x$monitor$n_eff)
    maxRhat <- round(max(x$monitor$Rhat), 3)
    cat(paste0("Minimum ESS = ", min_ess, " (",
               round(100 * min_ess / samples, 2), "%), and maximum Rhat = ",
               maxRhat, '\n'))
    if(min_ess < 200 | maxRhat > 1.1)
      warning("Signs of non-convergence found. Do not use for inference")
  }
  if(x$algorithm == "nuts"){
    extr <- extract_sampler_params(x)
    ndivs <- sum(extr[, "divergent__"])

    cat(paste0("There were ", ndivs, " divergences after warmup\n"))
  }
}

