#' Previews a vector of hexadecimal colors
#'
#' @param color_vector A hexadecimal color vector to plot
#'
#' @return
#' @export
#'
#' @examples
#' library(mkhelpers)
#' color_preview(c("#E04B4B", "#EC9393", "#3B537F", "#8998B2", "#969696", "#F2BF44"))
color_preview <- function(color_vector) {
  names(color_vector) <- color_vector
  tbl_x <- tibble(colors = factor(color_vector, levels = color_vector), value = 1)
  ggplot(tbl_x, aes(x = value, y = colors, fill = colors)) +
    geom_bar(position="stack", stat="identity") +
    theme_classic() +
    theme(legend.position = "none") +
    ggeasy::easy_remove_x_axis() +
    ggeasy::easy_remove_y_axis(what = c("ticks", "title", "line")) +
    scale_fill_manual(values = color_vector) +
    coord_fixed()
}
