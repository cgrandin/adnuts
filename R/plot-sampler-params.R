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
plot_sampler_params <- function(fit, plot = TRUE){

  sp <- extract_sampler_params(fit, inc_warmup = TRUE)
  sp_long <- data.frame(iteration = sp$iteration,
                        chain = factor(sp$chain),
                        value = c(sp$accept_stat__,
                                  log(sp$stepsize__),
                                  sp$n_leapfrog__,
                                  sp$divergent__,
                                  sp$energy__),
                        variable = rep(c("accept_stat",
                                         "log_stepsize",
                                         "n_leapfrog",
                                         "divergent",
                                         "energy"),
                                       each = nrow(sp)))
  g <- ggplot(sp_long,
              aes(iteration, y = value, color = chain)) +
    geom_point(alpha = 0.5) +
    facet_wrap("variable",
               scales = "free_y",
               ncol = 1) +
    theme_bw()

  if(plot){
    return(g)
  }

  invisible(g)
}
