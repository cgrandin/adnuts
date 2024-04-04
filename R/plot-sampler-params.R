#' Plot adaptation metrics for a fitted model.
#'
#' @details This utility function quickly plots the adaptation output of NUTS
#'  chains.
#'
#' @param fit A fitted object returned by [sample_admb()].
#' @param ... Arguments passed to [extract_sampler_params()]
#'
#' @return Prints and invisibly returns a ggplot object
#'
#' @importFrom rlang .data
#' @export
#' @examples
#' fit <- readRDS(system.file('examples', 'fit.RDS', package='adnuts'))
#' plot_sampler_params(fit)
plot_sampler_params <- function(fit,
                                label_size = 12,
                                label_face = c("plain",
                                               "bold",
                                               "italic"),
                                ...){

  label_face <- match.arg(label_face)

  # `inc_warmup` may be passed to `extract_sampler_params()`` using `...``
  sp <- extract_sampler_params(fit, ...)
  # Remove trailing double-underscores from labels
  names(sp) <- gsub("__", "", names(sp))
  # Replace lower case letters at beginning of labels with capitals
  names(sp) <- gsub("(?<=^|_)([a-z])",
                    "\\U\\1",
                    names(sp),
                    perl = TRUE)
  # Remove underscores from labels
  names(sp) <- gsub("_", " ", names(sp))

  sp_long <- sp |>
    mutate(Stepsize = log(Stepsize),
           Chain = factor(Chain)) |>
    rename(`Log Stepsize` = Stepsize) |>
    pivot_longer(cols = !c("Iteration", "Chain"))

  g <- sp_long |>
    ggplot(aes(Iteration,
               y = value,
               color = Chain)) +
    geom_point(alpha = 0.5) +
    facet_wrap("name",
               scales = "free_y",
               ncol = 1) +
    scale_color_viridis_d() +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = label_size,
                                    face = label_face))

  g
}
