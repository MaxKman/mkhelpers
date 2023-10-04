#' Project discrete or continous features on a umap plot
#'
#' @param tbl_x A tibble containing umap coordinates and feature columns
#' @param umap_dim_col_1 Name of column containing umap dimension 1
#' @param umap_dim_col_2 Name of column containing umap dimension 2
#' @param feature_x Name of column containing feature data
#' @param quantile_limits Quantile limits applied to continuous feature. E.g. c(0.1, 0.9) applies the dynamic range of the color scale only to to values above the 10th percentile and below the 90th percentile.
#' @param feature_colors For discrete features a color vector with one color per feature needs to be provided. A fixed assignment of colors to classes can be achieved by providing a named vector, where the vector names are the feature levels and the values the colors. For continous features provide either a vector of two or three colors to build a color gradient or leave this parameter empty to use the Cividis scale.
#' @param title Plot title.
#' @param output Return a ggplot object ("plot") or save as a png image ("image").
#' @param output_path Output path, shold include '.../imagename.png'.
#' @param dpi Resolution (dpi) for saved png.
#' @param order_values Whether to plot higher values in front of lower values ("sorted") or in random order ("random").
#' @param invert_sort_direction Whether to invert sort direction (plot lower values in front)
#' @param point_size Point size passed to geom_point().
#' @param alpha Alpha passed to geom_point().
#' @param plot_width Width of the coordinate system in mm. Set to NA to leave undetermined.
#' @param plot_height Height of the coordinate system in mm. Set to NA to leave undetermined.
#' @param out_width Width in mm of the png output.
#' @param out_height Height in mm of the png output.
#' @param show_legend Whether to show the ggplot legend.
#' @param show_labels Whether to label the mean coordinates of discrete feature values.
#' @param pretty_labels Whether to make labels pretty (takes longer)
#' @param pretty_labels_precision Tune the precision of label positioning. Higher values are more precise and take more processing time.
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
plot_umap <- function(tbl_x, umap_dim_col_1, umap_dim_col_2, feature_x, quantile_limits = c(0.3, 0.99), feature_colors = NULL, title, output = c("image", "plot"), output_path, dpi = 300, order_values = c("sorted", "random"), invert_sort_direction = FALSE, point_size = 0.3, point_size_legend = 1.5, alpha = 1, plot_width = 60, plot_height = 60, out_width = 89, out_height = 89, show_legend = TRUE, show_labels = FALSE, pretty_labels = FALSE, pretty_labels_precision = 50, label_size = 2) {

  if(order_values[[1]] == "sorted" & !invert_sort_direction) {
    tbl_x <- tbl_x %>%
      arrange({{feature_x}})
  } else if(order_values[[1]] == "sorted" & invert_sort_direction) {
    tbl_x <- tbl_x %>%
      arrange(desc({{feature_x}}))
  } else if(order_values[[1]] == "random") {
    tbl_x <- tbl_x %>%
      slice_sample(n = nrow(tbl_x), replace = FALSE)
  }

  p <- ggplot2::ggplot(tbl_x, aes({{umap_dim_col_1}}, {{umap_dim_col_2}}, color = {{feature_x}})) +
    mkhelpers::theme_mk +
    labs(x = "UMAP 1", y = "UMAP 2", title = title) +
    geom_point(stroke = 0, size = point_size, alpha = alpha)

  if(tbl_x %>% pull({{feature_x}}) %>% class %in% c("integer", "numeric")) {mode <- "continuous"}
  if(tbl_x %>% pull({{feature_x}}) %>% class %in% c("character", "factor")) {mode <- "discrete"}

  if(mode == "continuous") {
    feature_x_cutoffs <- tbl_x %>% mutate({{feature_x}} := na_if({{feature_x}}, 0)) %>% pull({{feature_x}}) %>% quantile(quantile_limits, na.rm = TRUE) #zero values are not taken into account when determining the cutoff
    p <- p +
      guides(color = guide_colourbar(barwidth = unit(0.9*plot_width, "mm"), barheight = unit(0.025*plot_width, "mm"), title.position = "top", title.hjust = 0.5))
    if(is.null(feature_colors)) {
      p <- p +
        viridis::scale_color_viridis(limits = feature_x_cutoffs, oob=scales::squish, option = "cividis")
    } else if (length(feature_colors) == 1) {
      print("Only a single color provided, defaulting to Cividis scale")
      p <- p +
        viridis::scale_color_viridis(limits = feature_x_cutoffs, oob=scales::squish, option = "cividis")
    }
    else if (length(feature_colors) == 2) {
      p <- p +
        scale_color_gradient(low = feature_colors[[1]], high = feature_colors[[2]], limits = feature_x_cutoffs, oob=scales::squish)
    } else {
      p <- p +
        scale_color_gradient2(low = feature_colors[[1]], mid = feature_colors[[2]], high = feature_colors[[3]], limits = feature_x_cutoffs, oob=scales::squish)
    }
  }

  if(mode == "discrete") {
    p <- p + guides(color = guide_legend(override.aes = list(size = point_size_legend), title.position = "top", title.hjust = 0.5)) # This fixes legend key size
    if(!is.null(feature_colors)) {
      p <- p +
        scale_color_manual(values = feature_colors)
    }
    if(show_labels & pretty_labels) {
      print("Building pretty labels, stay tuned ...")
      umap_1_min <- min(tbl_x %>% pull({{umap_dim_col_1}}))
      umap_1_max <- max(tbl_x %>% pull({{umap_dim_col_1}}))
      umap_2_min <- min(tbl_x %>% pull({{umap_dim_col_2}}))
      umap_2_max <- max(tbl_x %>% pull({{umap_dim_col_2}}))

      # Helper function for euclidean distance
      euclidean_dist <- function(x1, y1, x2, y2) {
        sqrt((x2 - x1)^2 + (y2 - y1)^2)
      }

      # Define a grid of points
      grid_points <- expand.grid(seq(umap_1_min, umap_1_max, length.out = pretty_labels_precision), seq(umap_2_min, umap_2_max, length.out = pretty_labels_precision)) %>%
        rename(label_coord_1 = Var1, label_coord_2 = Var2) %>%
        mutate(point_id = row_number())

      # Function to find nearest grid point
      find_grid_point <- function(coords) {
        grid_points %>%
          rowwise() %>%
          mutate(euclidean_dist = euclidean_dist(coords[[1]], coords[[2]], label_coord_1, label_coord_2)) %>%
          ungroup %>%
          slice_min(euclidean_dist, with_ties = FALSE) %>%
          pull(point_id)
      }

      # Find nearest grid point for a subset of each cluster and identify points with most cells nearby as reference for labelling
      label_df_tmp <- tbl_x %>%
        group_by({{feature_x}}) %>%
        sample_n(pretty_labels_precision, replace = TRUE) %>%
        rowwise() %>%
        mutate(point_id = find_grid_point(c({{umap_dim_col_1}}, {{umap_dim_col_2}}))) %>%
        group_by(point_id, {{feature_x}}) %>%
        mutate(n_cells = n()) %>%
        slice_max(n_cells, n = ceiling(nrow(grid_points) * 0.9), with_ties = TRUE) # This ensures that the labels are not trying to dodge individual cells situated between clusters

      label_df <- label_df_tmp %>%
        select(point_id, {{feature_x}}, n_cells) %>%
        distinct %>%
        group_by({{feature_x}}) %>%
        slice_max(n_cells, with_ties = FALSE) %>%
        ungroup %>%
        left_join(grid_points) %>%
        rename({{umap_dim_col_1}} := label_coord_1, {{umap_dim_col_2}} := label_coord_2) %>%
        select(-point_id, -n_cells) %>%
        ungroup %>%
        bind_rows(label_df_tmp %>% select({{umap_dim_col_1}}, {{umap_dim_col_2}}) %>% mutate({{feature_x}} := "") %>% relocate({{feature_x}}))

        # Label the clusters
        p <- p + ggrepel::geom_text_repel(data = label_df, aes(label = {{feature_x}}), size = label_size, color = "black", min.segment.length = 0, segment.size = 0.2, box.padding = 0.3, max.overlaps = nrow(tbl_x)*0.5, max.time = 10, max.iter = 1000000, force = 1)

      } else if(show_labels) {
        label_df <- tbl_x %>% select({{umap_dim_col_1}}, {{umap_dim_col_2}}, {{feature_x}}) %>%
          group_by({{feature_x}}) %>%
          summarise({{umap_dim_col_1}} := mean({{umap_dim_col_1}}), {{umap_dim_col_2}} := mean({{umap_dim_col_2}}))
        p <- p + geom_text(data = label_df, aes(label = {{feature_x}}), size = label_size, color = "black")
      }
  }

  if(!show_legend) {
    p <- p + theme(legend.position = "none")
  }

  if(show_legend) {
    p <- p + theme(legend.position = "bottom", legend.margin = margin(unit = "mm"))
  }

  if(!is.na(plot_width) & !is.na(plot_height)) {
    p <- p + plot_layout_mm(width = plot_width, height = plot_height)
  }

  if(output[[1]] == "image") {
    png_save_show(p, glue::glue("{output_path}", width = out_width, height = out_height), dpi = dpi)
  } else if (output[[1]] == "plot") {
    return(p)
  }
}
