#' Write a covariance matrix to `admodel.cov`.
#'
#' @param path Path to model.
#' @param cov_unbounded The covariance matrix in unbounded space.
#' @param hbf The hybrid_bounded_flag value. Use `hbf = 1` for HMC.
#' @param fn_cov Name of the ADMB covariance file
write_admb_cov <- function(path = NULL,
                           cov_unbounded,
                           hbf = NULL,
                           fn_cov = "admodel.cov"){

  if(!dir.exists(path)){
    stop("`path` does not exist",
         call. = FALSE)
  }
  fn <- file.path(path, fn_cov)
  if(!file.exists(fn)){
    stop("The file `", fn, "` does not exist",
         call. = FALSE)
  }

  dest_fn <- file.path(path, "admodel_original.cov")
  tmp <- file.copy(from = fn,
                   to = dest_fn)

  # Read in the output files
  results <- read_admb_cov(path)
  if(is.null(hbf)){
    hbf = results$hybrid_bounded_flag
  }
  scale <- results$scale
  num_pars <- results$num_pars
  if(nrow(cov_unbounded) != num_pars)
    stop("Invalid size of covariance matrix, has ", nrow(cov_unbounded),
         " rows but should have ", num_pars,
         call. = FALSE)

  # Write it to file using original scales, although these are ignored.
  file_new <- file(fn, "wb")
  on.exit(close(file_new))
  writeBin(as.integer(num_pars), file_new)
  writeBin(as.vector(as.numeric(cov_unbounded)), file_new)
  writeBin(as.integer(hbf), file_new)
  writeBin(as.vector(scale), file_new)
  invisible()
}
