#' Check that the  model is compiled with the right version of ADMB which is
#' 12.0 or later by default
#'
#' @details
#' Some functionality of packages [adnuts] is embedded in the ADMB
#' source code so that when a model is compiled it is contained in the model
#' executable. If this code does not exist, `adnuts` will fail. The solution
#' is to update ADMB and recompile the model.
#'
#' @param exe Executable filename, compiled using ADMB
#' @param check_path_env Logical. If `TRUE`, extract the full path for the
#' simple executable name given by `exe` from the $PATH environment
#' variable and use that for `exe`. If `FALSE`, use `exe` as the full
#' path executable name.
#' @param min_version Minimum valid version (numeric). Defaults to 12.0.
#' @param warn Logical. If `TRUE`, show warnings
#'
#' @return The version number as a character string, invisibly
#' @export
check_admb_version <- function(exe = NULL,
                               check_path_env = TRUE,
                               min_version = 12,
                               warn = TRUE){

  if(is.null(exe)){
    if(!check_path_env){
      stop("If `exe is `NULL`, `check_path_env`cannot be `FALSE`",
           call. = FALSE)

    }
  }

  fn_exe <- exe
  if(check_path_env){
    fn_exe <- suppressWarnings(system_(paste("which", exe),
                                       wait = TRUE,
                                       intern = TRUE))
    if(!is.null(attributes(fn_exe)) && attr(fn_exe, "status") == 1){
      stop("Operating system did not find the `exe` = `", exe,
           "` on the $PATH.\nCommand run was `", paste("which", exe),
           "`",
           call. = FALSE)
    }
  }

  if(!file.exists(fn_exe))
    stop("File `fn` = `", fn_exe, "` does not exist",
         call. = FALSE)

  # If executable in current dir and on linux, need to prepend dot slash
  if(dirname(fn_exe) == "." && get_os() != "windows"){
    fn_exe <- paste0("./", fn_exe)
  }

  # Get version number from the exe
  ver <- suppressWarnings(system_(paste(fn_exe, "-version"),
                                  wait = TRUE,
                                  intern = TRUE))
  if(!is.null(attributes(ver)) && attr(ver, "status") %in% 1:2){
    stop("There was an error running the version command. Check the executable ",
         "from a terminal. The command was:\n", c,
         call. = FALSE)
  }


  ver_line <- grep("^ADMB-.*compiled", ver)
  if(!length(ver_line) || length(ver_line) > 1){
    stop("Could not find the unique line of output with the version number.\n",
         "The command run was:\n", paste(fn_exe, "-version"), "\n\n",
         "Which output:\n", paste(ver, collapse = "\n"), "\n\n",
         "The regular expression used was:\n",
         "grep('^ADMB-.*compiled', ver)\n\n",
         "Which mathced the following lines of the output:\n",
         paste(ver_line, collapse = ", "),
         call. = FALSE)
  }

  line <- ver[ver_line]
  v <- gsub("^ADMB-([0-9]+\\.[0-9]+).*$", "\\1", line)
  v_num <- as.numeric(v)

  if(is.na(v_num) | !is.numeric(v_num)){
    stop("Could not verify the version number. Attempted to extract it ",
         "from the following line:\n\n",
         line, "\n\nIf it is present in that line, the regular ",
         "expression needs updating, cnotact package maintainer",
         call. = FALSE)
  }
  if(v_num < min_version)
    stop("`", exe, "` was compiled with old version of ADMB. ",
         "Version 12.0 or greater is required, found: `", v, "`\n",
         "`adnuts` is incompatible with this version. Update ADMB ",
         "and try again",
         call. = FALSE)

  invisible(v)
}

