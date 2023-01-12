#' Convert adnuts fit (named list) into a `shinystan` object.
#'
#' @details The `shinystan` package provides several conversion functions
#' for objects of different types, such as stanfit classes (Stan ouput)
#' and simple arrays. For the latter, option NUTS information, such as
#' [sampler_params()] can be passed. This function essentially extends
#' the functionality of [shinystan::as.shinystan()] to work specifically with
#' fits from adnuts (TMB or ADMB). The user can thus explore their model
#' with [launch_shinystan(as_shiny_nuts(fit))] in the same way
#' that Stan models are examined.
#' @param fit Output list from [sample_admb()].
#' @seealso [launch_shinyadmb()]
#' @importFrom shinystan as.shinystan
#' @return An S4 object of class `shinystan`. Depending on the algorithm
#' used, this list will have slight differences.
as_shiny_nuts <- function(fit){

  fit$algorithm <- toupper(fit$algorithm)
  if(fit$algorithm == "NUTS"){
    sso <- with(fit,
                as.shinystan(samples,
                             warmup = warmup,
                             max_treedepth = max_treedepth,
                             sampler_params=sampler_params,
                             algorithm = "NUTS",
                             model_name = model))
  }else if(fit$algorithm == "HMC"){
    sso <- with(fit,
                as.shinystan(samples,
                             warmup = warmup,
                             sampler_params = sampler_params,
                             algorithm = "HMC",
                             model_name = model))
  }else{
    sso <- with(fit,
                as.shinystan(samples,
                             warmup = warmup,
                             algorithm = "RWM",
                             model_name = model))
  }

  invisible(sso)
}
