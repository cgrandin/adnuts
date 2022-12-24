#' A wrapper for running ADMB models in parallel
#'
#' @param parallel_number The chain number
#' @param path The path name to append the chain number to
#' @param algorithm NUTS or RWM
#' @param ... Arguments passed to the [fit()] function
#'
#' @return Output from the [fit()] function
#' @export
sample_admb_parallel <- function(parallel_number,
                                 path,
                                 algorithm,
                                 ...){

  # pad the chain number to two digits
  chain_num_char <- as.character(parallel_number)
  nc <- nchar(chain_num_char)
  chain_num_char <- if(nc == 1,
                       paste0("0", chain_num_char),
                       chain_num_char)
  chain_dir <- file.path(path, paste0("chain_", chain_num_char))
  if(dir.exists(chain_dir)){
    unlink(chain_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(chain_dir)
  if(!dir.exists(chain_dir)){
    stop("Could not create directory: ", chain_dir, call. = FALSE)
  }

  trash <- file.copy(from = list.files(path, full.names = TRUE),
                     to = chain_dir)

  if(algorithm == "NUTS"){
    fit <- sample_admb_nuts(path = chain_dir,
                            chain = parallel_number,
                            ...)
  }else if(algorithm == "RWM"){
    fit <- sample_admb_rwm(path = chain_dir,
                           chain = parallel_number,
                           ...)
  }else{
    stop("Algorithm must be either 'NUTS' or 'RWM'", call. = FALSE)
  }
  unlink(chain_dir,
         recursive = TRUE,
         force = TRUE)

  fit
}
