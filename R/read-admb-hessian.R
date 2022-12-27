#' Read in admodel.hes file
#'
#' @details
#' This function reads in all of the information contained in the
#' admodel.hes file. Some of this is needed for relaxing the
#' covariance matrix, and others just need to be recorded and
#' rewritten to file so ADMB "sees" what it's expecting.
#'
#' @param path Path to directory containing the admodel.HES file
#'
#' @return The Hessian matrix
read_admb_hessian <- function(path){

  if(!dir.exists(path)){
    stop("`path` does not exist",
         call. = FALSE)
  }
  fn <- file.path(path, "admodel.hes")
  if(!file.exists(fn)){
    stop("The file `", fn, "` does not exist",
         call. = FALSE)
  }

  f <- file(fn, "rb")
  on.exit(close(f))

  num_pars <- readBin(f, "integer", 1)
  hes_vec <- readBin(f, "numeric", num_pars ^ 2)
  hes <- matrix(hes_vec, ncol = num_pars, nrow = num_pars)

  hybrid_bounded_flag <- readBin(f, "integer", 1)
  scale <- readBin(f, "numeric", num_pars)

  hes
}
