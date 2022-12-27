#' Return the Operating System name in lower case
#'
#' @return A character string describing the OS, e.g. 'windows',
#' 'linux', or 'osx' (Macbook)
#' @export
get_os <- function(){

  sysinf <- as.list(tolower(Sys.info()))
  names(sysinf) <- tolower(names(sysinf))

  if(is.null(sysinf)){
    os <- .Platform$OS.type
    if(grepl("^darwin", R.version$os)){
      os <- "osx"
    }
    if(grepl("linux-gnu", R.version$os)){
      os <- "linux"
    }
  }else{
    os <- sysinf$sysname
    if(os == "Darwin")
      os <- "osx"
  }

  tolower(os)
}
