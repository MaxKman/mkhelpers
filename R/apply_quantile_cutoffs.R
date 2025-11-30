#' Applying quantile cutoffs to a numeric vector
#'
#' @description Elements in the vector <= the lower cutoff will be set to the lower cutoff.
#' @description Elements in the vector >= the upper cutoff will be set to the upper cutoff.
#' @description All other elements will be left unchanged
#' @param x A numeric vector
#' @param quantiles Lower and upper quantiles to apply
#'
#' @return A numeric vector
#' @export
#'
#' @examples
#' library(mkhelpers)
#' mtcars$disp %>% sort
#' apply_quantile_cutoffs(mtcars$disp) %>% sort
apply_quantile_cutoffs <- function(x, quantiles = c(0.1, 0.9)) {
  quant_lower <- quantile(x, quantiles[[1]])
  quant_higher <- quantile(x, quantiles[[2]])
  x <- if_else(x <= quant_lower,
          quant_lower,
          if_else(x >= quant_higher,
                  quant_higher,
                  x))
  unname(x)
}
