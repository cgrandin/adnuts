#' Write an ADMB PSV file
#'
#' @details Useful to combine multiple MCMC runs together into a single
#' `.psv file` which can then be executed with '-mceval'.
#'
#' @param path Directory to write the PSV file in
#' @param fn_psv Model PSV file name with or without the PSV extension
#' @param samples A matrix or data.frame of samples, each column is a
#' parameter, each row a sample.
write_psv <- function(path,
                      fn_psv,
                      samples){

  if(!dir.exists(path)){
    stop("Directory `path` = `", path, "` does not exist",
         call. = FALSE)
  }
  if(!length(grep("\\.(psv)|(PSV)$", fn_psv))){
    fn_psv <- paste0(fn_psv, ".psv")
  }
  fn <- file.path(path, fn_psv)

  samples <- as.matrix(samples)
  con <- file(fn, "wb")
  on.exit(close(con), add = TRUE)
  writeBin(object = ncol(samples), con)
  writeBin(object = as.vector(t(samples)), con)
}
