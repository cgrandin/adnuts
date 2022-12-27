#' Check identifiability from model Hessian
#'
#' @param path Path to model folder, that contains `admodel.hes`
#' @param model Model name without file extension
#' @details Read in the admodel.hes file and check the eigenvalues to
#' determine which parameters are not identifiable and thus cause the
#' Hessian to be non-invertible. Use this to identify which parameters
#' are problematic. This function was converted from a version in the
#' `FishStatsUtils` package.
#' @return Prints output of bad parameters and invisibly returns it.
#' @export
check_identifiable <- function(model, path){

  # Check eigen decomposition
  fit <- read_mle_fit(model, path)
  hes <- read_admb_hessian(path)
  ev <-eigen(hes)
  which_bad <- which(ev$values < sqrt(.Machine$double.eps))
  if(length(which_bad) == 0){
    message( "All parameters are identifiable" )
  } else {
    # Check for parameters
    if(length(which_bad == 1)){
      row_max <- abs(ev$vectors[, which_bad])
    } else {
      row_max  <-  apply(ev$vectors[, which_bad],
                         MARGIN = 1,
                         FUN = function(vec){max(abs(vec))} )
    }
    bad <- data.frame(ParNum = 1:nrow(hes),
                      Param = fit$par.names,
                      MLE = fit$est[1:nrow(hes)],
                      Param_check = ifelse(row_max > 0.1, "Bad", "OK"))
    row.names(bad) <- NULL
    bad <- bad[bad$Param_check == "Bad", ]
    print(bad)
    invisible(bad)
  }
}
