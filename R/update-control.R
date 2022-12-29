#' Update the control list.
#'
#' @param control A list passed from a sampling function
#' @return A list with default control elements updated by those supplied
#' in `control`
update_control <- function(control){

  lst <- list(adapt_delta = 0.8,
              metric = "unit",
              stepsize = NULL,
              adapt_mass = TRUE,
              adapt_mass_dense = FALSE,
              max_treedepth = 12)

  # Special case if user is doing mle they probably don't want
  # mass adaptation turned on. They have to override it by
  # setting TRUE for either adaptation option
  if(is.character(control$metric) | is.matrix(control$metric)){
    if(is.null(control$adapt_mass) & is.null(control$adapt_mass_dense)){
      lst$adapt_mass <- lst$adapt_mass_dense <- FALSE
    }
  }

  if(!is.null(control)){
    for(i in names(control)){
      lst[[i]] <- control[[i]]
    }
  }
  if(lst$adapt_mass_dense & lst$adapt_mass){
    lst$adapt_mass <- FALSE
  }

  lst
}
