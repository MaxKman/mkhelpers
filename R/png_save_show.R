#' Export ggplot as png and include in knitted document
#'
#' @description Wrapper for bro::bro_ggsave_paged with unit set to "mm" and automatic inclusion of the png image in the knitted document
#' @param plot A ggplot object
#' @param file The output path. Filename has to end in png.
#' @param show_plot Whether to include the plot in the knitted document.
#' @param ... Other parameters to be passed to bro_ggsave_paged
#'
#' @export
#'
#' @examples
png_save_show <- function(plot, file, show_plot = TRUE, ...) {
  bro::bro_ggsave_paged(gg = plot,
                        filename = file,
                        units = "mm",
                        ...)
  if(show_plot) {
    knitr::include_graphics(file)
  }
}
