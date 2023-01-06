#' Trim a logfile down in size by removing all the MCMC step outputs
#'
#' @param fn The filename to modify
#'
#' @return Nothing
#' @export
trim_logfile <- function(fn){

  if(!file.exists(fn)){
    stop("File `", fn, "` does not exist",
         call. = FALSE.)
  }
  d <- readLines(fn)
  d <- d[!grepl("^ *MCMC:", d)]
  writeLines(d, fn)

  invisible()
}

#' Trim logfiles down in size by removing all the MCMC step outputs
#'
#' @details
#' Searches a `path` for all files matching `chain_XX_model_output.log`
#' where XX are digits representing chain numbers
#' @param path The path where the logfiles reside
#'
#' @return Nothing
#' @importFrom purrr walk
#' @export
trim_logfiles <- function(path){

  if(!dir.exists(path)){
    stop("Directory `", path, "` does not exist",
         call. = FALSE.)
  }
  pat <- "chain_[0-9]+_model_output\\.log$"
  fns <- list.files(path,
                    full.names = TRUE,
                    pattern = pat)
  if(!length(fns)){
    stop("No files in directory `", path, "` match pattern `", pat, "`",
         call. = FALSE)
  }

  walk(fns, ~{
    trim_logfile(.x)
  })

  invisible()
}