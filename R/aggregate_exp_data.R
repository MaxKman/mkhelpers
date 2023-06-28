#' Title Aggregate single cell expression data in pseudobulks
#' @description This function aggregates single cell data according to a predefined group (see aggr_col), separately for each donor
#' @param m A single cell matrix with genes as rows and cells as columns
#' @param md Metadata dataframe (must contain aggr_col, sample_col)
#' @param aggr_col The metadata column determining which cells to aggregate (e.g. clusters or experimental groups)
#' @param sample_col The metadata column identifying the sample to which each cell belongs (e.g. the donor)
#' @param n_cells_min The minimum number of cells for a cluster/sample combination to be included in the analysis. E.g. if set to 20, samples will be excluded for clusters, for which they contain less than 20 cells. Set to 0 to ignore this setting.
#' @param min_n_samples_aggr Minimum number of samples in each aggregation group. If less than min_n_samples_aggr aggregation will not be performed.
#' @param mode Whether to aggregate using mean() ("mean") or sum() ("sum").
#' @param n_cells_normalize The number of cells to which the pseudobulking will be standardized (only applied for mode = "sum"). E.g. if set to 10000, the counts aggregated from a cluster/sample combination will be normalized as follows: count_matrix / n_cells * n_cells_normalize.
#' @param return_matrix If set to TRUE returns a matrix with each aggregation groups as columns and genes as rows. If set to FALSE returns a dataframe
#' @param expr_format Return gene expression information in wide or long format (default: "wide")
#' @param cell_name_col The name of the metadata column that contains the cell identifier
#' @param subset_col Provide a name of metadata column here to use for subsetting the data
#' @param subset Provide identifiers, found in subset_col, on which to subset (e.g. a character vector of cell names or cluster identities)
#' @param verbose Whether to print messages. Default: TRUE
#'
#' @return A matrix with aggregation groups as columns and genes as rows (return_matrix = TRUE), or a dataframe containing the same information  (return_matrix = FALSE)
#' @export
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(SeuratData)
#' InstallData("pbmc3k")
#' seu_pbmc <- LoadData("pbmc3k", "pbmc3k.final")
#'
#' # Simulate some donors
#' add_donors <- seu_pbmc@meta.data %>%
#'   rownames_to_column("cell_name") %>%
#'   group_by(seurat_clusters) %>%
#'   mutate(simulated_donors = sample(5, n(), replace = TRUE)) %>%
#'   ungroup %>%
#'   column_to_rownames("cell_name") %>%
#'   select(simulated_donors)
#' add_donors$cell_name <- rownames(add_donors)
#' seu_pbmc <- AddMetaData(seu_pbmc, add_donors)
#'
#' # Basic example with mode = "mean"
#' aggregate_exp_data(m = seu_pbmc@assays$RNA@data, md = seu_pbmc@meta.data, aggr_col = seurat_clusters, sample_col = simulated_donors, n_cells_min = 20, min_n_samples_aggr = 3, mode = "mean")
#'
#' # Basic example with mode = "count"
#' aggregate_exp_data(m = seu_pbmc@assays$RNA@data, md = seu_pbmc@meta.data, aggr_col = seurat_clusters, sample_col = simulated_donors, n_cells_min = 20, n_cells_normalize = 10000, min_n_samples_aggr = 3, mode = "count")
#'
#' # Basic example with mode = "count", return dataframe
#' aggregate_exp_data(m = seu_pbmc@assays$RNA@data, md = seu_pbmc@meta.data, aggr_col = seurat_clusters, sample_col = simulated_donors, n_cells_min = 20, n_cells_normalize = 10000, min_n_samples_aggr = 3, mode = "count", return_matrix = FALSE)
#'
#' # Test if min_n_samples_aggr is applied correctly
#' md_modified <- seu_pbmc@meta.data %>% filter(!(seurat_clusters == 1 & simulated_donors %in% 1:3))
#' test_results <- aggregate_exp_data(m = seu_pbmc@assays$RNA@data, md = md_modified, aggr_col = seurat_clusters, sample_col = simulated_donors, n_cells_min = 20, n_cells_normalize = 10000, min_n_samples_aggr = 3, mode = "count", return_matrix = FALSE, expr_format = "long")
#' test_results$seurat_clusters %>% unique %>% sort
#'
#' # Test if n_cells_min is applied correctly
#' md_modified <- seu_pbmc@meta.data %>% filter(!(seurat_clusters == 1 & simulated_donors %in% 1:3))
#' md_sample <- seu_pbmc@meta.data %>% filter(seurat_clusters == 1 & simulated_donors %in% 1:3) %>% group_by(seurat_clusters, simulated_donors) %>% slice_sample(n = 19)
#' md_modified <- bind_rows(md_modified, md_sample)
#' test_results <- aggregate_exp_data(m = seu_pbmc@assays$RNA@data, md = md_modified, aggr_col = seurat_clusters, sample_col = simulated_donors, n_cells_min = 20, n_cells_normalize = 10000, min_n_samples_aggr = 3, mode = "count", return_matrix = FALSE, expr_format = "long")
#' test_results$seurat_clusters %>% unique %>% sort
#' test_results <- aggregate_exp_data(m = seu_pbmc@assays$RNA@data, md = md_modified, aggr_col = seurat_clusters, sample_col = simulated_donors, n_cells_min = 18, n_cells_normalize = 10000, min_n_samples_aggr = 3, mode = "count", return_matrix = FALSE, expr_format = "long")
#' test_results$seurat_clusters %>% unique %>% sort
#'  }
aggregate_exp_data <- function(m, md, aggr_col, sample_col = none, n_cells_min, n_cells_normalize = 0, min_n_samples_aggr, mode = c("mean", "sum"), return_matrix = TRUE, cell_name_col = cell_name, expr_format = c("wide", "long"), subset = "none", subset_col, verbose = TRUE) {

  if(is.null(dim(m))) {
    stop("The dimensions of the input matrix are undefined. Have you tried to subset a matrix on a single row or column? Try adding drop = FALSE to your subsetting call, e.g. m[1, , drop = FALSE]")
  }

  #Subset data
  if(subset != "none") {
    md <- md %>% filter({{subset_col}} %in% subset)
  }

  # First create a little report which samples are going to be included / excluded:
  group_stats <- md %>%
    group_by({{aggr_col}}, {{sample_col}}) %>%
    summarise(n_cells = n()) %>%
    ungroup() %>%
    filter(n_cells >= n_cells_min) %>%
    group_by({{aggr_col}}) %>%
    summarise(n_samples = n()) %>%
    filter(n_samples >= min_n_samples_aggr)
  all_groups <- md %>%
    pull({{aggr_col}}) %>%
    unique
  included_groups <- group_stats %>% pull({{aggr_col}}) %>% unique
  excluded_groups <- all_groups[!(all_groups %in% included_groups)]

  if(verbose) {
    print(glue::glue("The following aggregation groups have sufficient cells and samples to be aggregated:\n{included_groups %>% stringr::str_replace_na() %>% stringr::str_c(collapse = ', ')}\n\nThe following aggregation groups are excluded as they contain less than {min_n_samples_aggr} samples with >= {n_cells_min} cells:\n{excluded_groups %>% stringr::str_replace_na() %>% stringr::str_c(collapse = ', ')}\n\n"))
    }


  sample_col_str <- deparse(substitute(sample_col))
  cell_name_col_str <- deparse(substitute(cell_name_col))

  groups <- md %>% pull({{aggr_col}}) %>% unique

  # Remove samples with less than n_cells_min cells
  aggregated_exp_list <- map(groups, function(group_x) {
    if(sample_col_str != "none") {
      samples_cells <- md %>%
        filter({{aggr_col}} == group_x) %>%
        select({{sample_col}}, {{cell_name_col}}) %>%
        group_nest({{sample_col}}, .key = cell_name_col_str) %>%
        mutate(n_cells = map({{cell_name_col}}, nrow)) %>%
        filter(n_cells >= n_cells_min) %>%
        mutate({{cell_name_col}} := map({{cell_name_col}}, deframe))

      # Only proceed with groups that still have min_n_samples_aggr (otherwise return NA)
      if(nrow(samples_cells) >= min_n_samples_aggr) {
        m_aggr_group <- map2(samples_cells %>% pull({{sample_col}}), samples_cells %>% pull({{cell_name_col}}), function(sample_x, cells_sample_x) {
          if(verbose) {print(glue::glue("... aggregating {group_x}, {sample_x} ..."))}
          m <- m[, cells_sample_x, drop=FALSE]
          if(mode[1] == "mean") {
            m_aggr <- Matrix::rowMeans(m)
          } else {
            n_cells_x <- ncol(m)
            m_aggr <- Matrix::rowSums(m)
            if(n_cells_normalize > 0) {
              m_aggr <- m_aggr / n_cells_x * n_cells_normalize
              m_aggr <- round(m_aggr)
            }
          }
          m_aggr <- tibble(gene = names(m_aggr), exp = m_aggr, {{aggr_col}} := group_x, {{sample_col}} := sample_x)
        })
        bind_rows(m_aggr_group)
      } else {
        return(NA)
      }
    } else {
      group_cells <- md %>% filter({{aggr_col}} == group_x) %>% pull({{cell_name_col}})
      if(length(group_cells) > n_cells_min) {
        m <- m[, group_cells, drop=FALSE]
        if(mode[1] == "mean") {
          m_aggr <- Matrix::rowMeans(m)
        } else {
          n_cells_x <- ncol(m)
          m_aggr <- Matrix::rowSums(m)
          if(n_cells_normalize > 0) {
            m_aggr <- m_aggr / n_cells_x * n_cells_normalize
            m_aggr <- round(m_aggr)
          }
        }
        m_aggr <- tibble(gene = names(m_aggr), exp = m_aggr, {{aggr_col}} := group_x)
      } else {
        return(NA)
      }
    }
  })
  aggregated_exp_list <- aggregated_exp_list[!is.na(aggregated_exp_list)]
  aggregated_exp <- bind_rows(aggregated_exp_list)

  if(return_matrix & sample_col_str != "none") {
    aggregated_exp %>% unite(group_sample, c({{aggr_col}}, {{sample_col}}), sep = "_._") %>% pivot_wider(names_from = group_sample, values_from = exp) %>% column_to_rownames("gene") %>% as.matrix()
  } else if (return_matrix & sample_col_str == "none") {
    aggregated_exp %>% pivot_wider(names_from = {{aggr_col}}, values_from = exp) %>% column_to_rownames("gene") %>% as.matrix()
  } else if (expr_format[[1]] == "long") {
    aggregated_exp
  } else {
    aggregated_exp %>% pivot_wider(names_from = {{aggr_col}}, values_from = exp)
  }
}
