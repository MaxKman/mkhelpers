#' Export ggplot as png and include in knitted document
#'
#' @description Wrapper for ggsave with unit set to "mm" and automatic inclusion of the png image in the knitted document
#' @param plot A ggplot object
#' @param file The output path. Filename has to end in png.
#' @param show_plot Whether to include the plot in the knitted document.
#' @param ... Other parameters to be passed to ggsave
#'
#' @export
#'
#' @examples
#' library(mkhelpers)
#' p <- ggplot(mtcars, aes(mpg, cyl)) +
#'   geom_point()
#'   png_save_show(p, "~/mtcars.png", width = 89, height = 89)

png_save_show <- function(plot, file, show_plot = TRUE, ...) {
  ggplot2::ggsave(plot = plot,
                  filename = file,
                  units = "mm",
                  limitsize = FALSE,
                  ...)
  if(show_plot) {
    knitr::include_graphics(file)
  }
}

