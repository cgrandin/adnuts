#' Read in an ADMB PSV file
#'
#' @param path The path to the location of the PSV file
#' @param fn_psv The name of the psv file, with or without the extension
#' psv or .PSV
#' @param nms Column names for the output. If `NULL`, they will be sequential
#' starting with the letter 'V'
#'
#' @return A data frame with the PSV file output
read_psv <- function(path,
                     fn_psv = NULL,
                     nms = NULL){

  if(!dir.exists(path)){
    stop("Directory `path` = `", path, "` does not exist",
         call. = FALSE)
  }
  if(is.null(fn_psv)){
    stop("PSV file name `fn_psv` is `NULL`",
         call. = FALSE)
  }
  if(!length(grep("\\.(psv)|(PSV)$", fn_psv))){
    fn_psv <- paste0(fn_psv, ".psv")
  }
  fn <- file.path(path, fn_psv)

  if(!file.exists(fn)){
    # Sometimes ADMB will shorten the name of the psv file,
    # so search for it here
    fns <- list.files(path)
    mtch <- grep("\\.(psv)|(PSV)$", fns)
    if(!length(mtch)){
      stop("No .psv file found", call. = FALSE)
    }
    if(length(mtch) > 1){
      stop("More than one .psv file found", call. = FALSE)
    }

    ff <- fns[mtch]
    warning("PSV file `",
            fn,
            "` not found, using only PSV file found: ", ff)
    fn <- file.path(path, ff)
  }

  # Read in the binary file
  f <- file(fn, open = "rb")
  nv <- readBin(f, "int")
  fs <- file.info(fn)$size
  isize <- 4
  dsize <- 8
  f_contents <- matrix(readBin(f, "double", n = (fs - isize) / dsize),
                       byrow = TRUE,
                       ncol = nv)
  close(f)

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
