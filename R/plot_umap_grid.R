#' Project a list of discrete or continous features on a umap plot
#'
#' @param tbl_x A tibble containing umap coordinates and feature columns
#' @param feature_list List of feature columns to be plotted
#' @param quantile_limits_list List of quantile limits for each plot (e.g. list(c(0.1,0.9), c(0.3,0.99))). Will only be applied to features with continuous scale
#' @param n_cols Number of columns to plot
#' @param out_width Width of image output in mm
#' @param output_path Output path and file name (should end in .png)
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
plot_umap_grid <- function(tbl_x, feature_list, quantile_limits_list = NULL, n_cols = 3, out_width = 180, output_path, ...) {
  n_plots <- length(feature_list)
  col_width <- out_width / n_cols
  out_height <- ceiling(n_plots / n_cols) * col_width * 0.8
  if(is.null(quantile_limits_list)) {
    quantile_limits <- list(c(0.1, 0.9))
    quantile_limits_list <- rep(quantile_limits, length(feature_list))
  }
  plot_l <- map2(feature_list, quantile_limits_list, function(feature_x, quantile_limits_x) {
    plot_umap(tbl_x = tbl_x, feature_x = !!sym(feature_x), plot_width = 0.6 * col_width, plot_height = 0.6 * col_width, output = "plot", title = feature_x, quantile_limits = quantile_limits_x, ...)
  })
  p <- patchwork::wrap_plots(plot_l, ncol = n_cols, widths = col_width)
  png_save_show(p, output_path, width = out_width, height = out_height)
}
