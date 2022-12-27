#' Extract posterior samples from a model fit.
#'
#' A helper function to extract posterior samples across multiple chains
#' into a single data.frame.
#'
#' @details This function is loosely based on the \pkg{rstan} function
#'   \code{extract}. Merging samples across chains should only be used for
#'   inference after appropriate diagnostic checks. Do not calculate
#'   diagnostics like Rhat or effective sample size after using this
#'   function, instead, use \code{\link[rstan]{monitor}}. Likewise, warmup
#'   samples are not valid and should never be used for inference, but may
#'   be useful in some cases for diagnosing issues.
#'
#' @param fit A list returned by \code{sample_admb}.
#' @param inc_warmup Whether to extract the warmup samples or not
#'   (default). Warmup samples should never be used for inference, but may
#'   be useful for diagnostics.
#' @param inc_lp Whether to include a column for the log posterior density
#'   (last column). For diagnostics it can be useful.
#' @param as_list Whether to return the samples as a list (one element per
#'   chain). This could then be converted to a CODA mcmc object.
#' @param unbounded Boolean flag whether to return samples in
#'   unbounded (untransformed) space. Will only be differences
#'   when init_bounded types are used in the ADMB template. This
#'   can be useful for model debugging.
#' @return If as.list is FALSE, an invisible data.frame containing samples
#'   (rows) of each parameter (columns). If multiple chains exist they will
#'   be rbinded together, maintaining order within each chain. If as.list
#'   is TRUE, samples are returned as a list of matrices.
#' @importFrom purrr map
#' @export
#' @examples
#' ## A previously run fitted ADMB model
#' fit <- readRDS(system.file('examples', 'fit.RDS', package='adnuts'))
#' post <- extract_samples(fit)
#' tail(apply(post, 2, median))
extract_samples <- function(fit,
                            inc_warmup = FALSE,
                            inc_lp = FALSE,
                            as_list = FALSE,
                            unbounded = FALSE){

  if(!is_adfit(fit)){
    stop("`fit` is not a valid `adfit` object",
         call. = FALSE)
  }
  if(unbounded){
    x <- fit$samples_unbounded
    if(is.null(x))
      stop("No unbounded parameters in this fit",
           call. = FALSE)
  }else{
    x <- fit$samples
    if(is.null(x))
      stop("No posterior samples found",
           call. = FALSE)
  }
  if(!is.array(x))
    stop("`fit$samples` is not an array",
         call. = FALSE)

  # Set up ind to be indices to remove the first `warmup` samples
  ind <- -seq_len(fit$warmup)
  if(inc_warmup){
    ind <- seq_len(dim(x)[1])
  }

  if(inc_lp){
    y <- map(seq_len(dim(x)[2]), ~{x[ind, .x, ]})
  } else {
    y <- map(seq_len(dim(x)[2]), ~{x[ind, .x, -dim(x)[3]]})
  }

  if(as_list){
    return(invisible(y))
  }

  invisible(as.data.frame(do.call(rbind, y)))
}
