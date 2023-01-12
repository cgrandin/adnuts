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
  if(!is.list(x)){
    stop("`fit$sampler_parameters` is not a list",
         call. = FALSE)
  }
  if(inc_warmup){
    ind <- seq_len(dim(x[[1]])[1])
    its <- seq_along(ind)
  }else{
    ind <- -seq_len(fit$warmup)
    its <- seq_len(ind) + fit$warmup
  }
  y <- do.call(rbind,
               lapply(seq_len(x), function(i){
                 cbind(chain = i, iteration = its, x[[i]][ind, ])
               }))

  invisible(as.data.frame(y))
}
