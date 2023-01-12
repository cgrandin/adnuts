#' Function to generate random initial values from a previous fit using
#' adnuts
#'
#' @param fit An outputted list from [sample_admb()]
#' @param chains The number of chains for the subsequent run, which
#' determines the number to return.
#' @return A list of lists which can be passed back into
#' [sample_admb()].
#' @export
sample_inits <- function(fit, chains){
  post <- extract_samples(fit)
  ind <- sample(seq_len(nrow(post)), size = chains)
  lapply(ind, function(i) as.numeric(post[i, ]))
}

#' Launch shinystan for an ADMB fit.
#'
#' @param fit A named list returned by [sample_admb()]
#' @seealso [launch_shinytmb()]
#' @export
launch_shinyadmb <- function(fit){
  launch_shinystan(as_shiny_nuts(fit))
}

#' Call [shell()] or [system()] depending on the Operating System
#'
#' @param ... Arguments to pass to either [shell()] or [system()]
#'
#' @return The output from the command function called
#' @export
system_ <- function(...){

  if(get_os() == "windows"){
    shell(...)
  }else{
    system(...)
  }
}
