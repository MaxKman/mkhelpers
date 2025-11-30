#' Wrapper for plot_layout with mm as default
#'
#' @param width Plot width in mm
#' @param height Plot height in mm
#' @param ... Parameters passed on to patchwork::plot_layout
#'
#' @return A patchwork
#' @export
#'
#' @examples
#' library(mkhelpers)
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(mpg, cyl)) + geom_point()
#' p + plot_layout_mm(30, 30)
plot_layout_mm <- function(width, height, ...) {
  patchwork::plot_layout(..., widths = grid::unit(width, "mm"), heights = grid::unit(height, "mm"))
}
