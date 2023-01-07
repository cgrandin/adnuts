#' Setup directories and run a sampler, either [sample_admb_rwm()]
#' or [sample_admb_nuts()]
#'
#' @rdname samplers

#' @param chain_num The chain number
#' @param ... Arguments passed to the [fit()] function
#'
#' @return Output from the [sample_admb_nuts()] and [sample_admb_rwm()]
#' functions
#'
#' @export
run_admb_sampler <- function(chain_num,
                             path,
                             algorithm = c("nuts", "rwm"),
                             fn_logfile = "model_output.log",
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
                            fn_logfile = fn_logfile,
                            ...)
  }else if(algorithm == "rwm"){
    fit <- sample_admb_rwm(path = chain_dir,
                           chain = chain_num,
                           fn_logfile = fn_logfile,
                           ...)
  }else{
    stop("Algorithm must be either `nuts` or `rwm`",
         call. = FALSE)
  }

  # Copy log files to the parent directory
  if(!is.null(fn_logfile)){
    fn_logfile_from <- file.path(chain_dir, fn_logfile)
    chain_num_char <- as.character(chain_num)
    if(nchar(chain_num_char) == 1){
      chain_num_char <- paste0("0", chain_num_char)
    }
    fn_logfile_to <- file.path(path,
                               paste0("chain_",
                                      chain_num_char,
                                      "_",
                                      fn_logfile))
    if(file.exists(fn_logfile_from)){
      file.copy(fn_logfile_from, fn_logfile_to, overwrite = TRUE)
    }
  }

  unlink(chain_dir,
         recursive = TRUE,
         force = TRUE)

  fit
}
