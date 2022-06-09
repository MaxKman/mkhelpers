#' Add group membership column for each group
#' @description This function is mainly designed to support plotting of cell membership for individual groups, for example individual clusters
#' @param tbl_x A tibble containing a cell name column and a group column
#' @param cell_name_col The name of the column that contains the cell identifier
#' @param group_col The column identifying group membership
#'
#' @return
#' @export
#'
#' @examples
#' library(Seurat)
#' library(SeuratData)
#' InstallData("pbmc3k")
#' seu_pbmc <- LoadData("pbmc3k", "pbmc3k.final")
#' tbl_x <- seu_extract_tbl(seu_x = seu_pbmc,
#'                          reduction = "umap",
#'                          metadata_cols = "seurat_clusters")
#' tbl_x <- tbl_x %>% add_group_id_columns(cell_name_col = cell_name, group_col = seurat_clusters)
#' show_labels_list <- c(TRUE, rep(FALSE, length(tbl_x$seurat_clusters %>% unique)))
#' plot_umap_grid(tbl_x = tbl_x, umap_dim_col_1 = UMAP_1, umap_dim_col_2 = UMAP_2, output_path = "~/test_plot_cluster_ids.png", feature_list = str_subset(colnames(tbl_x), "^seurat_clusters"), show_labels = show_labels_list, point_size = 0.5, show_legend = FALSE)
add_group_id_columns <- function(tbl_x, cell_name_col, group_col) {
  group_col_str <- deparse(substitute(group_col))
  tbl_temp <- tbl_x %>%
    select({{cell_name_col}}, {{group_col}}) %>%
    mutate(cluster_membership = factor("1", levels = c("1", "0"))) %>%
    pivot_wider(names_from = {{group_col}}, names_prefix = str_c(group_col_str, "_"), values_from = cluster_membership) %>%
    replace(is.na(.), "0")
  tbl_x %>% left_join(tbl_temp)
}
