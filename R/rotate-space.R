#' Update algorithm for mass matrix.
#'
#' @details
#' Rotation done using choleski decomposition
#' First case is a dense mass matrix
#'
#' @param fn The current fn function.
#' @param gr The current gr function
#' @param y_cur The current parameter vector in unrotated (Y) space.
#' @param M The new mass matrix
rotate_space <- function(fn, gr, M, y_cur){

  if(is.matrix(M)){
    # Lower triangular Cholesky decomp.
    chd <- t(chol(M))
    # Inverse
    chd_inv <- solve(chd)
    # Define rotated fn and gr functions
    fn2 <- function(x){
      fn(chd %*% x)
    }
    gr2 <- function(x){
      as.vector(gr(chd %*% x) %*% chd )
    }
    # Now rotate back to "x" space using the new mass matrix M
    x_cur <- as.numeric(chd_inv %*% y_cur)
  }else if(is.vector(M)){
    chd <- sqrt(M)
    fn2 <- function(x){
      fn(chd * x)
    }
    gr2 <- function(x){
      as.vector(gr(chd * x) ) * chd
    }
    # Now rotate back to "x" space using the new mass matrix M. M is a
    # vector here. Note the big difference in efficiency without the
    # matrix operations.
    x_cur <- 1 / chd * y_cur
  }else{
    stop("Mass matrix must be vector or matrix",
         call. = FALSE)
  }
  # Redefine these functions
  # Need to adjust the current parameters so the chain is
  # continuous. First rotate to be in Y space.
  list(gr2 = gr2,
       fn2 = fn2,
       x_cur = x_cur,
       chd = chd)
}
