#' Extract runtime, final acceptance ratio, and final stepsize
#' information from logffiles after running adnuts
#'
#' @details
#' The information is extracted from the logfiles, 1 for each
#' chain
#' @param path
#'
#' @return A list of numeric vectors including warmup and sampling runtimes,
#' final acceptance ratio, and final stepsize
#' @export
get_runtimes <- function(path){

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
  fns <- sort(fns)
  log_contents <- map(fns, ~{
    tmp <- readLines(.x)
    tail(tmp, 10)
  })

  get_time <- function(str, chain_num){
    # Is it in hours, minutes, or seconds
    secs <- length(grep("seconds", str))
    mins <- length(grep("minutes", str))
    hrs <- length(grep("hours", str))
    days <- length(grep("days", str))

    time_type <- 1000 * days + 100 * hrs + 10 * mins + secs
    if(!time_type){
      stop("No time information found in the logfile for chain ", chain_num,
           call. = FALSE)
    }

    if(!time_type %in% c(1, 10, 100, 1000)){
      stop("Problem with the time information in the logfile for ",
           "chain ", chain_num,
           call. = FALSE)
    }
    time_type <- as.character(time_type)
    switch(time_type,
           `1` = {
             tm <- as.numeric(gsub(".* ([0-9]+\\.[0-9]+) +seconds.*$",
                                   "\\1",
                                   str))
           },
           `10` = {
             tm <- as.numeric(gsub(".* ([0-9]+\\.[0-9]+) +minutes.*$",
                                   "\\1",
                                   str))
             tm <- tm * 60
           },
           `100` = {
             tm <- as.numeric(gsub(".* ([0-9]+\\.[0-9]+) +hours.*$",
                                   "\\1",
                                   str))
             tm <- tm * 3600
           },
           `1000` = {
             tm <- as.numeric(gsub(".* ([0-9]+\\.[0-9]+) +days.*$",
                                   "\\1",
                                   str))
             tm <- tm * 86400
           })
    # Time is now in seconds no matter how it was reported
    tm
  }

  out <- list()
  out$final_step_size <- imap_dbl(log_contents, ~{
    ind <- grep("Final step size", .x)
    if(!length(ind)){
      return(NA_real_)
    }
    as.numeric(gsub(".*Final step size *= *([0-9]+\\.[0-9]+);.*$", "\\1",
                    .x[ind]))
  })
  out$final_acceptance_ratio <- imap_dbl(log_contents, ~{
    ind <- grep("Final acceptance ratio", .x)
    if(!length(ind)){
      return(NA_real_)
    }
    as.numeric(gsub(".*Final acceptance ratio *= *([0-9]+\\.[0-9]+),.*$", "\\1",
                    .x[ind]))
  })
  out$warmup_time <- imap_dbl(log_contents, ~{
    ind <- grep("Warmup ", .x)
    if(!length(ind)){
      return(NA_real_)
    }
    get_time(.x[ind], .y)
  })
  out$sampling_time <- imap_dbl(log_contents, ~{
    ind <- grep("Sampling ", .x)
    if(!length(ind)){
      return(NA_real_)
    }
    get_time(.x[ind], .y)
  })
  out$total_time <- map_dbl(log_contents, ~{
    ind <- grep("Total ", .x)
    if(!length(ind)){
      return(NA_real_)
    }
    get_time(.x[ind], .y)
  })

  attr(out, "units") <- "seconds"
  out
}
