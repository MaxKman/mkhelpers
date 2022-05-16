#' Title
#'
#' @param tbl_x A tibble containing feature expression values and group designations
#' @param feature_list A list of features to be plotted from feature_col
#' @param feature_col The tibble column containing the feature names
#' @param n_cols Number of columns to plot
#' @param out_width Width of image output in mm
#' @param output_path Output path and file name (should end in .png)
#' @param ... Other parameters passed on to plot_umap. Note, that expr_cola and group_col need to be provided.
#'
#' @return
#' @export
#'
#' @examples
plot_feature_comparison_grid <- function(tbl_x, feature_list, feature_col = "gene", n_cols = 3, out_width = 180, output_path, ...) {
  n_plots <- length(feature_list)
  col_width <- out_width / n_cols
  out_height <- ceiling(n_plots / n_cols) * col_width * 0.9
  plot_l <- map(feature_list, function(feature_x) {
    tbl_temp = tbl_x %>% filter(!!sym(feature_col) == feature_x)
    plot_feature_comparison(tbl_x = tbl_temp, title = feature_x, output = "plot", output_path = output_path, ...)
  })
  p <- patchwork::wrap_plots(plot_l, ncol = n_cols, widths = col_width)
  png_save_show(p, output_path, width = out_width, height = out_height)
}
