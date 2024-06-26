% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-marginals.R
\name{plot_marginals}
\alias{plot_marginals}
\title{Plot marginal distributions for a fitted model}
\usage{
plot_marginals(
  fit,
  pars = NULL,
  mfrow = NULL,
  add_mle = TRUE,
  add_monitor = TRUE,
  breaks = 30
)
}
\arguments{
\item{fit}{A fitted object returned by
\code{\link[=sample_admb]{sample_admb()}}.}

\item{pars}{A numeric or character vector of parameters which
to plot, for plotting a subset of the total (defaults to all)}

\item{mfrow}{A custom grid size (vector of two) to be called
as \code{par(mfrow)}, overriding the defaults.}

\item{add_mle}{Whether to add marginal normal distributions
determined from the inverse Hessian file}

\item{add_monitor}{Whether to add ESS and Rhat information}

\item{breaks}{The number of breaks to use in \code{\link[graphics:hist]{graphics::hist()}},
defaulting to 30}
}
\description{
Plot marginal distributions for a fitted model
}
\details{
This function plots grid cells of all parameters
in a model, comparing the marginal posterior histogram vs
the asympototic normal (red lines) from the inverse
Hessian. Its intended use is to quickly gauge differences
between frequentist and Bayesian inference on the same
model.

If \code{fit$monitor} exists the effective sample size
(ESS) and R-hat estimates are printed in the top right
corner. See \link{https://mc-stan.org/rstan/reference/Rhat.html}
for more information. Generally Rhat > 1.05 or ESS < 100 (per chain)
suggest inference may be unreliable.

This function is customized to work with multipage PDFs,
specifically:
\code{pdf('marginals.pdf', onefile=TRUE, width=7,height=5)}
produces a nice readable file.
}
\examples{
fit <- readRDS(system.file('examples', 'fit.RDS', package='adnuts'))
plot_marginals(fit, pars=1:2)
}
