#' Read in the ADMB covariance file.
#
#' @param path Path to model
#'
#' @return A list of named values
read_admb_cov <- function(path = NULL){

  if(!dir.exists(path)){
    stop("`path` does not exist",
         call. = FALSE)
  }

  fn <- file.path(path, "admodel.cov")
  fn_bin <- file(fn, "rb")
  on.exit(close(fn_bin), add = TRUE)

  num_pars <- readBin(fn_bin, "integer", 1)
  cov_vec <- readBin(fn_bin, "numeric", num_pars ^ 2)
  cov_unbounded <- matrix(cov_vec, ncol = num_pars, nrow = num_pars)
  hybrid_bounded_flag <- readBin(fn_bin, "integer", 1)
  scale <- readBin(fn_bin, "numeric", num_pars)
  cov_bounded <- cov_unbounded * (scale %o% scale)

  list(num_pars = num_pars,
       cov_bounded = cov_bounded,
       cov_unbounded = cov_unbounded,
       hybrid_bounded_flag = hybrid_bounded_flag,
       scale = scale)
}
