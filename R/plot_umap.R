#' Project discrete or continous features on a umap plot
#'
#' @param tbl_x A tibble containing umap coordinates and feature columns
#' @param umap_dim_col_1 Name of column containing umap dimension 1
#' @param umap_dim_col_2 Name of column containing umap dimension 2
#' @param feature_x Name of column containing feature data
#' @param quantile_limits Quantile limits applied to continuous feature. E.g. c(0.1, 0.9) applies the dynamic range of the color scale only to to values above the 10th percentile and below the 90th percentile.
#' @param feature_colors Colors to be used for discrete features. Ideally provide as named vector, where the vector names are the feature levels and the values the colors.
#' @param title Plot title. Also used for naming the output file.
#' @param output Return a ggplot object ("plot") or save as a png image ("image").
#' @param output_path Output path
#' @param dpi Resolution (dpi) for saved png.
#' @param order_values Whether to plot higher values in front of lower values ("sorted") or in random order ("random").
#' @param point_size Point size passed to geom_point()
#' @param alpha Alpha passed to geom_point()
#' @param plot_width Width of the coordinate system in mm. Set to NA to leave undetermined.
#' @param plot_height Height of the coordinate system in mm. Set to NA to leave undetermined.
#' @param out_width Width in mm of the png output.
#' @param out_height Height in mm of the png output.
#' @param show_legend Whether to show the ggplot legend.
#' @param show_labels Whether to label the mean coordinates of discrete feature values.
#' @param label_size Size of labels.
#'
#' @return A ggplot object (when output is set to "plot")
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
#'    genes_extract = "CD3D",
#'    assay_extract = "RNA",
#'    slot_extract = "data",
#'    expr_format = "wide")
#'  plot_umap(tbl_x = tbl_x, umap_dim_col_1 = UMAP_1, umap_dim_col_2 = UMAP_2, output_path = "~", title = "Seurat clusters", feature_x = seurat_clusters, order_values = "random", show_legend = FALSE, show_labels = TRUE, point_size = 0.5)
#'  plot_umap(tbl_x = tbl_x, umap_dim_col_1 = UMAP_1, umap_dim_col_2 = UMAP_2, output_path = "~", title = "CD3 expression", feature_x = CD3D, point_size = 0.5)
#' }
plot_umap <- function(tbl_x, umap_dim_col_1, umap_dim_col_2, feature_x, quantile_limits = c(0.3, 0.99), feature_colors = NULL, title, output = c("image", "plot"), output_path, dpi = 300, order_values = c("sorted", "random"), point_size = 0.3, alpha = 1, plot_width = 60, plot_height = 60, out_width = 89, out_height = 89, show_legend = TRUE, show_labels = FALSE, label_size = 2) {

  if(order_values[[1]] == "sorted") {
    tbl_x <- tbl_x %>%
      arrange({{feature_x}})
  } else if(order_values[[1]] == "random") {
    tbl_x <- tbl_x %>%
      slice_sample(n = nrow(tbl_x), replace = FALSE)
  }

  p <- ggplot2::ggplot(tbl_x, aes({{umap_dim_col_1}}, {{umap_dim_col_2}}, color = {{feature_x}})) +
    mkhelpers::theme_mk +
    labs(x = "UMAP 1", y = "UMAP 2", title = title) +
    geom_point(stroke = 0, size = point_size, alpha = alpha)

  if(tbl_x %>% pull({{feature_x}}) %>% class == "numeric") {mode <- "continuous"}
  if(tbl_x %>% pull({{feature_x}}) %>% class %in% c("character", "factor")) {mode <- "discrete"}

  if(mode == "continuous") {
    feature_x_cutoffs <- tbl_x %>% pull({{feature_x}}) %>% quantile(quantile_limits)
    p <- p + viridis::scale_color_viridis(limits = feature_x_cutoffs, oob=scales::squish, option = "cividis")
  }

  if(mode == "discrete") {
    if(!is.null(feature_colors)) {
      p <- p + scale_color_manual(values = feature_colors)
    }
    if(show_labels) {
      label_df <- tbl_x %>% select({{umap_dim_col_1}}, {{umap_dim_col_2}}, {{feature_x}}) %>%
        group_by({{feature_x}}) %>%
        summarise({{umap_dim_col_1}} := mean({{umap_dim_col_1}}), {{umap_dim_col_2}} := mean({{umap_dim_col_2}}))
      p <- p + geom_text(data = label_df, aes(label = {{feature_x}}), size = label_size, color = "black")
    }
  }

  if(!show_legend) {
    p <- p + theme(legend.position = "none")
  }

  if(!is.na(plot_width) & !is.na(plot_height)) {
    p <- p + plot_layout_mm(width = plot_width, height = plot_height)
  }

  if(output[[1]] == "image") {
    png_save_show(p, glue::glue("{output_path}/{title}.png", width = out_width, height = out_height), dpi = dpi)
  } else if (output[[1]] == "plot") {
    return(p)
  }
}
