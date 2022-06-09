#' Project a list of discrete or continous features on a umap plot
#'
#' @param tbl_x A tibble containing umap coordinates and feature columns
#' @param feature_list List of feature columns to be plotted
#' @param quantile_limits Quantile limits applied to continuous feature. E.g. c(0.1, 0.9) applies the dynamic range of the color scale only to to values above the 10th percentile and below the 90th percentile. Set to a one item list, e.g. list(c(0.1, 0.9)) to affect all plots or a provide separate settings for each feature as a list.
#' @param n_cols Number of columns to plot
#' @param out_width Width of image output in mm
#' @param output_path Output path and file name (should end in .png)
#' @param show_labels Whether to label the mean coordinates of discrete feature values. Set to TRUE or FALSE to affect all plots or give a list of booleans for all features
#' @param show_legend Whether to show the ggplot legend. Set to TRUE or FALSE to affect all plots or give a list of booleans for all features.
#' @param order_values Whether to plot higher values in front of lower values ("sorted") or in random order ("random"). Provide a list to set individually for each plot.
#' @param invert_sort_direction Whether to invert sort direction (plot lower values in front). Provide a list to set individually for each plot.
#' @param ... Other parameters passed on to plot_umap. Note, that umap_dim_col_1 and umap_dim_col_2 need to be provided.
#'
#' @return
#' @export
#'
#' @examples
#'\dontrun{
#'  library(Seurat)
#'  library(SeuratData)
#'  InstallData("pbmc3k")
#'  seu_pbmc <- LoadData("pbmc3k", "pbmc3k.final")
#'
#'  tbl_x <- seu_extract_tbl(seu_x = seu_pbmc,
#'    reduction = "umap",
#'    metadata_cols = "seurat_clusters",
#'    extract_expr = TRUE,
#'    genes_extract = c("CD3D", "CD19"),
#'    assay_extract = "RNA",
#'    slot_extract = "data",
#'    expr_format = "wide")
#' plot_umap_grid(tbl_x = tbl_x, umap_dim_col_1 = UMAP_1, umap_dim_col_2 = UMAP_2, output_path = "~/test_plot_umap_grid.png", feature_list = c("seurat_clusters", "CD3D", "CD19"), show_labels = TRUE, point_size = 0.5)
#' }
plot_umap_grid <- function(tbl_x, feature_list, quantile_limits = list(c(0.3, 0.99)), n_cols = 3, out_width = 180, output_path, show_labels = FALSE, show_legend = TRUE, order_values = list(c("sorted", "random")), invert_sort_direction = FALSE, ...) {
  n_plots <- length(feature_list)
  col_width <- out_width / n_cols
  out_height <- ceiling(n_plots / n_cols) * col_width * 0.9

  plot_l <- pmap(list(feature_list, quantile_limits, show_labels, show_legend, order_values, invert_sort_direction), function(feature_x, quantile_limits_x, show_labels_x, show_legend_x, order_values_x, invert_sort_direction_x) {
    plot_umap(tbl_x = tbl_x, feature_x = !!sym(feature_x), plot_width = 0.6 * col_width, plot_height = 0.6 * col_width, output = "plot", title = feature_x, quantile_limits = quantile_limits_x, show_labels = show_labels_x, show_legend = show_legend_x, order_values = order_values_x, invert_sort_direction = invert_sort_direction_x, ...)
  })

  p <- patchwork::wrap_plots(plot_l, ncol = n_cols, widths = col_width)
  png_save_show(p, output_path, width = out_width, height = out_height)
}
