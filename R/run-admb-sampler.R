#' Setup directories and run a sampler, either [sample_admb_rwm()]
#' or [sample_admb_nuts()]
#'
#' @param chain_num The chain number
#' @param path The path name to append the chain number to
#' @param algorithm The algorithm to use, one of `nuts` or `rwm`
#' @param ... Arguments passed to the [fit()] function
#'
#' @return Output from the [fit()] function
#' @export
run_admb_sampler <- function(chain_num,
                             path,
                             algorithm = c("nuts", "rwm"),
                             ...){

  algorithm <- match.arg(algorithm)

  # pad the chain number to two digits
  chain_num_char <- as.character(chain_num)
  nc <- nchar(chain_num_char)
  chain_num_char <- ifelse(nc == 1,
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

  if(algorithm == "nuts"){
    fit <- sample_admb_nuts(path = chain_dir,
                            chain = chain_num,
                            ...)
  }else if(algorithm == "rwm"){
    fit <- sample_admb_rwm(path = chain_dir,
                           chain = chain_num,
                           ...)
  }else{
    stop("Algorithm must be either `nuts` or `rwm`",
         call. = FALSE)
  }
  unlink(chain_dir,
         recursive = TRUE,
         force = TRUE)

  fit
}
