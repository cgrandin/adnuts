#' Find the model executable name on your filesystem
#'
#' @description
#' Find the model executable name in the local directory (`loc_dir`).
#' If not there, search the `$PATH`.
#' If found more than once in the `$PATH`, the first one will be used and a
#' warning issued.
#' If not on the `$PATH`, throw an error and stop code execution.
#' If using Linux or OSX and the file permissions do not allow the user to
#' execute, throw an error and stop code execution.
#'
#' @param model The name of the executable (without any extension)
#' @param loc_dir The path to check for existence of the executable before
#' checking the path
#' @param path_only If `TRUE`, ignore `loc_dir` and search the PATH only
#'
#' @return The executable as it should be called in the OS. It is either
#' just as it was given if on the `$PATH` or in the `loc_dir` directory on a
#' Windows machine, or pre-pended with `./` if in the `loc_dir` directory on a
#' non-Windows machine. If not found anywhere, an error will be thrown and
#' code execution will stop
#' @export
get_model_executable <- function(model,
                                 loc_dir = getwd(),
                                 path_only = FALSE){

  stopifnot(is.character(loc_dir))
  stopifnot(is.character(model))

  if(!path_only && !dir.exists(loc_dir)){
    stop("Directory `", loc_dir, "` does not exist. Check ",
         "argument `loc_dir`",
         call. = FALSE)
  }

  check_dot <- grep("\\.", model)
  if(length(check_dot)){
    stop("model name cannot include a '.', Do not include extension in ",
         "model name",
         call. = FALSE)
  }

  os <- get_os()
  loc_dir_files <- dir(loc_dir)
  model_regex <- gsub("\\+", "\\\\+", model)
  model_name <- grep(paste0("^",
                            model_regex,
                            ifelse(os == "windows", ".exe", ""), "$"),
                     loc_dir_files,
                     value = TRUE)

  if(!path_only && length(model_name)){
    if(os != "windows"){
      model <- paste0("./", model)
    }
    found_model_loc <- "loc_dir"
  }else{
    # No local executable, check PATH for one by using which or where
    model_path <- tryCatch(system(paste("which", model),
                                  intern = TRUE),
                           warning = function(w){
                             NA
                           })
    if(is.na(model_path[1])){
      stop(model, " not found locally or on the PATH.", call. = FALSE)
    }
    if(length(model_path) > 1){
      warning("The model executable ", model,
              " was found in more than one place on your PATH.\n",
              "Using the first one found (", model_path[1], ").")
    }
    found_model_loc <- "path"
  }
  # Check to see if the file is executable (LInux and OSX only)
  if(os != "windows"){
    stat <- system(paste("stat", model_path[1]), intern = TRUE)
    perms <- grep("^Access:.*Uid.*$", stat, value = TRUE)
    perms <- gsub("\\).*", "", perms)
    perms <- gsub(".*/","",perms)
    first_perm <- substr(perms, 1, 1)
    user_exec_perm <- substr(perms, 4, 4)
    if(first_perm == "d"){
      stop("The model name given (", model_path[1], ") is a directory and not executable.",
           call. = FALSE)
    }else if(first_perm == "p"){
      stop("The model name given (", model_path[1], ") is a named pipe and not executable.",
           call. = FALSE)
    }else if(first_perm == "s"){
      stop("The model name given (", model_path[1], ") is a socket and not executable.",
           call. = FALSE)
    }else if(first_perm %in% c("b", "c")){
      stop("The model name given (", model_path[1], ") is a device block and not executable.",
           call. = FALSE)
    }else if(user_exec_perm != "x"){
      stop("The file found (", model_path[1], ") is not executable.\nCheck name and use chmod to ",
           "change permissions if name is correct.",
           call. = FALSE)
    }
  }
  if(found_model_loc == "loc_dir"){
    message("Using model found in the local directory given by the 'loc_dir' argument: ", model)
  }else if(found_model_loc == "path"){
    message("Using model executable found in the PATH: ", model)
  }
  model
}
