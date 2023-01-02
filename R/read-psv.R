#' Read in a binary ADMB PSV file
#'
#' @param path The path to the location of the PSV file
#' @param nms Column names for the output. If `NULL`, they will be sequential
#' starting with the letter 'V'
#'
#' @return A data frame with the PSV file output
read_psv <- function(path,
                     nms = NULL){

  if(!dir.exists(path)){
    stop("Directory `path` = `", path, "` does not exist",
         call. = FALSE)
  }

  fns <- list.files(path)
  psv_fn_ind <- grep("\\.(psv)|(PSV)$", fns)
  if(!length(psv_fn_ind)){
    stop("No `psv` file found in directory:\n",
         path,
         call. = FALSE)
  }
  if(length(psv_fn_ind) > 1){
    stop("More than one `psv` file found in directory:\n",
         path,
         call. = FALSE)
  }

  fn <- file.path(path, fns[psv_fn_ind])

  # Read in the binary file
  psv_conn <- file(fn, open = "rb")
  psv_ncol <- readBin(psv_conn, "int")
  fs <- file.info(fn)$size
  isize <- 4
  dsize <- 8
  f_contents <- matrix(readBin(psv_conn,
                               "double",
                               n = (fs - isize) / dsize),
                       byrow = TRUE,
                       ncol = psv_ncol)
  close(psv_conn)

  if(is.null(nms)){
    nms <- paste0("V", seq(ncol(f_contents)))
  }
  if(length(nms) < ncol(f_contents)){
    warning("`nms` too short, more columns than names. Generating dummy `nms`")
    nms <- c(nms, paste0("V", seq(length(nms) + 1, ncol(f_contents))))
  }else if(length(nms) < ncol(f_contents)){
    warning("`nms` too long, more names than columns. Truncating `nms`")
    nms <- nms[seq(ncol(f_contents))]
  }
  colnames(f_contents) <- nms
  as.data.frame(f_contents)
}
