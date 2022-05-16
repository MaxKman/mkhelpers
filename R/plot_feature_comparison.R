#' Title
#'
#' @param tbl_x A tibble containing feature expression values and group designations
#' @param title The plot title (e.g. the gene that is plotted)
#' @param expr_col The column, which contains the expression information
#' @param group_col The column, which designates the group (e.g. single cell clusters or treatment groups)
#' @param colors Here one can set the group colors manually by providing a named character vector (e.g. c(group_1 = "#666666", group_2 = "#FFC966"))
#' @param output Return a ggplot object ("plot") or save as a png image ("image").
#' @param output_path Output path.
#' @param dpi Resolution (dpi) for saved png.
#' @param point_size Point size passed to geom_jitter()
#' @param alpha Alpha passed to geom_jitter().
#' @param plot_width Width of the coordinate system in mm. Set to NA to leave undetermined.
#' @param plot_height Height of the coordinate system in mm. Set to NA to leave undetermined.
#' @param out_width Width in mm of the png output.
#' @param out_height Height in mm of the png output.
#'
#' @return A ggplot object (when output is set to "plot")
#' @export
#'
#' @examples
#'# See ?group_markers
plot_feature_comparison <- function(tbl_x, title, expr_col, group_col, colors = "auto", output = c("image", "plot"), output_path, dpi = 600, point_size = 1.2, alpha = 0.8, plot_width = 25, plot_height = 15, out_width = 40, out_height = 30) {
  p <- ggplot(tbl_x, aes(!!sym(group_col), !!sym(expr_col), color = !!sym(group_col))) +
    geom_boxplot(outlier.shape = NA, size = 0.2, color = "#7A7A7A", width = 0.8) +
    geom_jitter(width = 0.3, height = 0, stroke = 0, size = point_size, alpha = alpha) +
    theme_mk +
    theme(legend.position = "none") +
    labs(title = title)
  if(!is.na(plot_width) & !is.na(plot_height)) {
    p <- p + plot_layout_mm(width = plot_width, height = plot_height)
  }
  if(colors[[1]] != "auto") {
    p <- p + scale_color_manual(values = colors)
  }
  if(output[[1]] == "image") {
    png_save_show(p, output_path, dpi = dpi, width = out_width, height = out_width)
  } else {
    return(p)
  }
}
