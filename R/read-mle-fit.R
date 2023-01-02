#' Read maximum likelihood fit for ADMB model
#'
#' @details
#' This is based loosely on [r4ss::read.admbFit()]
#' Sequentially read .par file which contains model size, minimum NLL, and
#' maxgrad at the top.
#'
#' @param path The path containing the MLE output you wish to load (contains
#' a `par` and `cor` file)
#'
#' @return A list containing, MLE estimates, standard errors, covariance
#' and correlation matrices, and other output from ADMB
#'
#' @export
read_mle_fit <- function(path = NULL){
  # if(model == "ss3"){
  #   model <- "ss"
  # }

  fns <- list.files(path)
  par_fn_ind <- grep("\\.par$", fns)
  if(!length(par_fn_ind)){
    stop("No `par` file found in directory:\n",
         path,
         call. = FALSE)
  }
  if(length(par_fn_ind) > 1){
    stop("More than one `par` file found in directory:\n",
         path,
         call. = FALSE)
  }
  par_fn <- file.path(path, fns[par_fn_ind])
  par_header <- scan(par_fn,
                     what = "",
                     n = 16,
                     quiet = TRUE)
  par <- par_header[which(par_header == "=") + 1]
  if(length(par) != 3){
    stop("The `par` file did not contain exactly three equals signs in ",
         "the header row which means one of `number of parameters`, ",
         "`objective value`, or `maximum gradient` were missing:\n",
         par_fn,
         call. = FALSE)
  }
  nopar <- par[1]
  nll <- par[2]
  maxgrad <- par[3]

  # The .cor file contains parameter (and derived quantity) names,
  # estimates, and se's. This is more convenient to read in than the .par
  # file.
  cor_fn_ind <- grep("\\.(cor)|(COR)$", fns)
  if(!length(cor_fn_ind)){
    stop("No `cor` file found in directory:\n",
         path,
         call. = FALSE)
  }
  if(length(cor_fn_ind) > 1){
    stop("More than one `cor` file found in directory:\n",
         path,
         call. = FALSE)
  }

  cor_fn <- file.path(path, fns[cor_fn_ind])

  cor_content <- readLines(cor_fn)
  # Remove top teo header lines
  tot_non_header <- cor_content[3:length(cor_content)]
  tot_par <- length(tot_non_header)
  if(tot_par < nopar){
    warning("`cor` file contains `", tot_par, "` parameters and `par` file ",
            "contiains `", nopar, "` parameters. MLE output not available (`NULL`)",
            call. = FALSE)
    return(NULL)
  }

  lines <- map(strsplit(tot_non_header, " "), ~{.x[.x != ""]})
  names_all <- map_chr(lines, ~{.x[2]})
  par_names <- names_all[seq_len(nopar)]

  est <- map_dbl(lines, ~{as.numeric(.x[3])})
  se <- map_dbl(lines, ~{as.numeric(.x[4])})

  # The correlation in the bounded space.
  cor <- matrix(NA, tot_par, tot_par)
  cor_vec <- imap(lines, ~{as.numeric(.x[5:(.y + 4)])}) |>
    unlist()
  cor[upper.tri(cor, diag = TRUE)] <- as.numeric(cor_vec)
  cor[lower.tri(cor)] <- t(cor)[lower.tri(cor)]

  # Covariance matrix
  # cov <- cor*(std %o% std)
  list(nopar = nopar,
       nll = nll,
       maxgrad = maxgrad,
       par_names = par_names,
       names_all = names_all,
       est = est,
       est = est,
       se = se,
       cor = cor[1:nopar, 1:nopar])
}
